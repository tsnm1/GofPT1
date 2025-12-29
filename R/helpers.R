# --- 3.1 Data Splitting Helper ---
split_data <- function(X, Y) {
  X <- as.matrix(X)
  n <- nrow(X)
  idx <- sort(sample(1:n, floor(n / 2)))
  list(
    x1 = X[-idx, , drop = FALSE], y1 = Y[-idx],
    x2 = X[idx, , drop = FALSE],  y2 = Y[idx]
  )
}

# --- 3.2 Estimation Function (Crucial: Handles CPB vs Hybrid Logic) ---
#' @details
#' If \code{penalize} is \code{NULL} (Hybrid mode), the function mimics the original code's
#' logic: it attempts Lasso selection followed by a GLM refit. If the GLM coefficients contain NAs
#' (indicating collinearity), it uses \code{caret::findCorrelation} to remove highly correlated
#' predictors (cutoff 0.9) and refits.
get_residuals_and_beta <- function(x, y, family = "gaussian", penalize = TRUE) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  n <- nrow(x)
  p <- ncol(x)
  pred <- NULL
  beta_hat <- rep(0, p)
  intercept <- 0

  # --- Logic Branch A: CPB Method (Strictly GLMNET based) ---
  if (!is.null(penalize)) {
    if (penalize) {
      # 1. Penalized: Lasso Selection -> Unpenalized GLMNET Refit
      fit_cv <- cv.glmnet(x, y, family = family)
      beta_lasso <- coef(fit_cv, s = "lambda.min")[-1]
      idx <- which(beta_lasso != 0)

      if (length(idx) == 0) {
        pred <- rep(mean(y), n)
        # beta_hat remains 0
      } else {
        fit_refit <- glmnet(x[, idx, drop = FALSE], y, family = family, lambda = 0)
        pred <- predict(fit_refit, newx = x[, idx, drop = FALSE], type = "response")
        beta_sub <- coef(fit_refit)[-1]
        beta_hat[idx] <- as.numeric(beta_sub)
        intercept <- coef(fit_refit)[1]
      }
    } else {
      # 2. Unpenalized: Full GLMNET
      fit_full <- glmnet(x, y, family = family, lambda = 0)
      pred <- predict(fit_full, newx = x, type = "response")
      beta_hat <- coef(fit_full)[-1]
      intercept <- coef(fit_full)[1]
    }
  }

  # --- Logic Branch B: Hybrid Method (penalize = NULL) ---
  else {
    # Original Pcvm_Pls logic: Lasso Selection -> GLM (stats) -> Caret Cleanup
    glm_fam <- if (family == "gaussian") gaussian() else binomial(link = "logit")

    cv_fit <- cv.glmnet(x, y, family = family, intercept = TRUE)
    b_lasso <- as.numeric(coef(cv_fit, s = "lambda.min")[-1])
    intercept <- as.numeric(coef(cv_fit, s = "lambda.min")[1])
    idx <- which(b_lasso != 0)
    beta_hat <- b_lasso # Default to Lasso beta

    if (length(idx) > 0) {
      x_sub <- x[, idx, drop = FALSE]

      # Attempt standard GLM
      fit_sub <- tryCatch(glm(y ~ x_sub, family = glm_fam), error = function(e) NULL)

      b_sub_check <- if (!is.null(fit_sub)) coef(fit_sub)[-1] else NA

      # If GLM fails or produces NA (Collinearity), use Caret logic
      if (is.null(fit_sub) || any(is.na(b_sub_check))) {
        # Original logic: "if(sum(is.na(sec_beta1))>0) { library(caret) ... }"
        if (ncol(x_sub) > 1) {
          # Use pairwise complete obs for correlation, just like original code
          cor_matrix <- cor(x_sub, use = "pairwise.complete.obs")
          high_cor <- tryCatch(caret::findCorrelation(cor_matrix, cutoff = 0.9), error = function(e) NULL)

          if (!is.null(high_cor) && length(high_cor) > 0) {
            x_sub <- x_sub[, -high_cor, drop = FALSE]
            # Update index to remove dropped columns
            idx <- idx[-high_cor]
            # Refit GLM
            fit_sub <- tryCatch(glm(y ~ x_sub, family = glm_fam), error = function(e) NULL)
          }
        }
      }

      # Apply GLM results if valid
      if (!is.null(fit_sub) && !any(is.na(coef(fit_sub)[-1]))) {
        beta_hat[] <- 0 # Reset
        beta_hat[idx] <- coef(fit_sub)[-1]
        intercept <- coef(fit_sub)[1]
        pred <- predict(fit_sub, type = "response")
      } else {
        # Fallback to Lasso prediction if GLM still fails
        pred <- predict(cv_fit, newx = x, s = "lambda.min", type = "response")
      }
    } else {
      # Intercept only model
      pred <- rep(mean(y), n)
    }
  }

  # Final fallback for predictions
  if (is.null(pred)) {
    lin <- as.numeric(x %*% beta_hat + intercept)
    pred <- if (family == "gaussian") lin else plogis(lin)
  }
  pred <- as.vector(pred)
  U <- y - pred

  # Calculate beta projection (Logic from original code)
  bnorm <- sqrt(sum(beta_hat^2))
  b_pro <- if (bnorm > 1e-10) beta_hat / bnorm else rep(1, p) / sqrt(p)

  # Calculate derivative for Binomial Martingale test (Logic from original code)
  deri <- NULL
  if (family == "binomial") {
    # Original: deri_link1 <- plogis(...) * (1 - plogis(...))
    lin_pred <- x %*% beta_hat + intercept
    probs <- plogis(lin_pred)
    deri <- as.vector(probs * (1 - probs))
  }

  return(list(
    U = as.matrix(U), beta_pro = matrix(b_pro, ncol = 1),
    beta_hat = beta_hat, intercept = intercept, deri_link = deri
  ))
}

# --- 3.3 SDR Projections ---
get_sdr_projections <- function(x, y_list, family) {
  n <- nrow(x)
  p <- ncol(x)
  screen_num <- floor(n / log(n))

  run_sdr <- function(tar, is_cat) {
    tryCatch(
      {
        beta <- NULL
        if (p <= screen_num) {
          fit <- LassoSIR::LassoSIR(x, tar, categorical = is_cat, screening = FALSE, H = 10, solution.path = FALSE, nfolds = 5)
          beta <- fit$beta
        } else {
          met <- if (is_cat) "MV-SIS" else "DC-SIS"
          rank <- VariableScreening::screenIID(x, tar, method = met)$rank
          idx <- which(rank <= screen_num)
          if (length(idx) > 0) {
            fit <- LassoSIR::LassoSIR(x[, idx, drop = FALSE], tar, categorical = is_cat, screening = FALSE, H = 10, solution.path = FALSE, nfolds = 5)
            beta <- matrix(0, p, ncol(fit$beta))
            beta[idx, ] <- fit$beta
          }
        }
        beta <- as.matrix(beta)
        if (ncol(beta) == 0) {
          return(matrix(NA, p, 1))
        }
        norm <- sqrt(colSums(beta^2))
        norm[norm < 1e-12] <- 1
        return(t(t(beta) / norm))
      },
      error = function(e) matrix(NA, p, 1)
    )
  }

  return(list(
    U_pro = run_sdr(y_list$U, FALSE),
    y_pro = run_sdr(y_list$original_y, family == "binomial")
  ))
}

# --- 3.4 Martingale Test Statistics ---
calc_martingale_pvals <- function(x_source, est, pro_target, family) {
  n <- nrow(x_source)
  num_pro <- ncol(pro_target)
  pvals <- numeric(num_pro)
  U <- as.vector(est$U)

  for (q in 1:num_pro) {
    x_pro <- as.vector(x_source %*% pro_target[, q])
    # Original logic: Indictor <- ifelse(x_pro <= t(x_pro), 1, 0)
    Ind <- outer(x_pro, x_pro, "<=") * 1

    # --- Branch A: Linear (Gaussian) ---
    if (family == "gaussian") {
      hat_A <- rbind(t(x_pro), rep(1, n))

      G_inv <- array(0, c(2, 2, n))
      for (i in 1:n) {
        # Original: ginv((hat_A%*%diag(Ind[i,]))%*%t(hat_A))
        G_inv[, , i] <- MASS::ginv((hat_A %*% diag(Ind[i, ])) %*% t(hat_A))
      }

      Integ <- hat_A %*% diag(as.vector(U)) %*% t(Ind)

      sec <- matrix(0, n, n)
      for (l in 1:n) {
        sec[l, ] <- (t(hat_A[, l]) %*% G_inv[, , l] %*% Integ[, l]) * Ind[l, ]
      }

      m_sta <- (1 / sqrt(n)) * (t(U) %*% Ind - colSums(sec))

      ord <- sort(x_pro)
      t_val <- ord[floor(0.99 * n)]
      sigma2 <- mean(U^2)
      F_val <- (mean(x_pro <= t_val))^2
      stat <- (1 / (sigma2 * max(F_val, 1e-10))) * mean((x_pro <= t_val) * (as.vector(m_sta)^2))

      pvals[q] <- pvalue_integ_Brown(stat)
    }

    # --- Branch B: Logit (Binomial) ---
    else {
      deri <- est$deri_link # Already calculated as p(1-p) in helper

      # The third term in band_y: (x %*% beta_hat) * deri
      term3 <- as.vector((x_source %*% est$beta_hat + est$intercept) * deri)
      band_y <- cbind(U^2, deri, term3)

      # Original logic: bandwidth_choice specific to logit
      h <- bandwidth_choice_logit(x_pro, band_y)

      # Kernel Construction (Epanechnikov)
      dist <- outer(x_pro, x_pro, "-") / h
      K <- (3 / 4) * (1 - dist^2) * (abs(dist) <= 1)

      # Original A matrices calculation
      denom <- as.vector(K %*% U^2) + n^(-12)
      A1 <- ((K %*% deri) * x_pro) / denom
      A2 <- (K %*% term3) / denom
      A3 <- (K %*% deri) / denom
      hat_A <- rbind(t(A1), t(A2), t(A3))

      G_inv <- array(0, c(3, 3, n))
      for (i in 1:n) {
        # Original logic: ((rep(1,3)%*%t(U)^2)*hat_A) ...
        term <- ((rep(1, 3) %*% t(U)^2) * hat_A) %*% (t(hat_A) * (as.matrix(Ind[i, ]) %*% rep(1, 3))) + n^(-12)
        G_inv[, , i] <- MASS::ginv(term)
      }

      Integ <- ((rep(1, 3) %*% t(U)) * hat_A) %*% t(Ind)

      sec <- matrix(0, n, n)
      for (l in 1:n) {
        val <- (U[l]^2 * t(hat_A[, l])) %*% G_inv[, , l] %*% Integ[, l]
        sec[l, ] <- as.numeric(val) * Ind[l, ]
      }

      m_sta <- (1 / sqrt(n)) * (t(U) %*% Ind - colSums(sec))

      ord <- sort(x_pro)
      t_val <- ord[floor(0.99 * n)]
      # Original: psi_1 <- (1/n1) * sum(U1^2 * (x1_pro <= t_1))
      psi <- mean(U^2 * (x_pro <= t_val))
      # Original: (1/psi_1^2) * mean(...)
      stat <- (1 / max(psi^2, 1e-10)) * mean(U^2 * (x_pro <= t_val) * (as.vector(m_sta)^2))

      pvals[q] <- pvalue_integ_Brown(stat)
    }
  }
  return(pvals)
}

# --- 3.5 PLS Test Statistics ---
calc_pls_pvals <- function(x_source, U, pro_target) {
  n <- nrow(x_source)
  h <- n^(-2 / 9)
  num_pro <- ncol(pro_target)
  pvals <- numeric(num_pro)
  EE <- U %*% t(U)

  for (q in 1:num_pro) {
    x_pro <- as.vector(x_source %*% pro_target[, q])
    dist <- outer(x_pro, x_pro, "-") / h
    # Epanechnikov kernel: (3/4)*(1-x^2)*I(|x|<=1)
    K <- (3 / 4) * (1 - dist^2) * (abs(dist) <= 1)

    KEE <- K * EE
    numer <- sum(KEE) - sum(diag(KEE))
    denom_sq <- sum(KEE^2) - sum(diag(KEE^2))
    Tn <- numer / sqrt(2 * denom_sq + 1e-12)
    pvals[q] <- 1 - pnorm(Tn)
  }
  return(pvals)
}

# --- 3.6 P-value Combination (with Truncation) ---
combine_pvals <- function(p) {
  p <- as.vector(p[!is.na(p)])
  k <- length(p)
  if (k == 0) {
    return(c(cauchy = NA, fisher = NA, min_p = NA))
  }

  # --- Truncation Logic (Safety Belt) ---
  # Cauchy Protection: prevent tan(pi/2) = Inf
  p_c <- p
  p_c[p_c > 0.999999] <- 0.999999
  p_c[p_c < 1e-16] <- 1e-16
  cauchy <- 1 - pcauchy(mean(tan((0.5 - p_c) * pi)))

  # Fisher Protection: prevent log(0) = -Inf
  p_f <- p
  p_f[p_f < 1e-300] <- 1e-300
  fisher <- 1 - pchisq(-2 * sum(log(p_f)), df = 2 * k)

  # Min-P
  min_p <- 1 - (1 - min(p))^k

  return(c(cauchy = cauchy, fisher = fisher, min_p = min_p))
}

# --- 3.7 Helpers for Bandwidth and Brownian Table ---
bandwidth_choice_logit <- function(x, y_mat) {
  n <- length(x)
  # Original code: c_h <- seq(0.25, 1.25, 0.15)
  c_h <- seq(0.25, 1.25, 0.15)
  h_seq <- c_h * n^(-2 / 9)
  CV <- numeric(length(h_seq))

  for (i in seq_along(h_seq)) {
    h <- h_seq[i]
    dist <- outer(x, x, "-") / h
    K <- (3 / 4) * (1 - dist^2) * (abs(dist) <= 1)
    diag(K) <- 0
    denom <- rowSums(K) + n^(-12)

    cv_err <- 0
    # Cross Validation error sum across the columns of band_y
    for (k in 1:ncol(y_mat)) {
      est <- (K %*% y_mat[, k]) / denom
      cv_err <- cv_err + mean((y_mat[, k] - est)^2)
    }
    CV[i] <- cv_err
  }
  return(h_seq[which.min(CV)])
}

pvalue_integ_Brown <- function(x) {
  p0 <- c(
    1, 1, 0.9994, 0.9945, 0.9824, 0.9642, 0.9417, 0.9169, 0.891, 0.8648, 0.839, 0.8138, 0.7894, 0.7659,
    0.7434, 0.7218, 0.7012, 0.6814, 0.6626, 0.6445, 0.6273, 0.6108, 0.5949, 0.5798, 0.5652, 0.5513,
    0.5378, 0.5249, 0.5125, 0.5006, 0.489, 0.4779, 0.4672, 0.4568, 0.4468, 0.4371, 0.4278, 0.4187,
    0.4099, 0.4014, 0.3931, 0.3851, 0.3773, 0.3697, 0.3623, 0.3552, 0.3482, 0.3414, 0.3348, 0.3284,
    0.3222, 0.3161, 0.3101, 0.3043, 0.2987, 0.2932, 0.2878, 0.2825, 0.2774, 0.2724, 0.2675, 0.2627,
    0.258, 0.2534, 0.2489, 0.2446, 0.2403, 0.2361, 0.232, 0.228, 0.224, 0.2202, 0.2164, 0.2127, 0.2091,
    0.2056, 0.2021, 0.1987, 0.1953, 0.1921, 0.1889, 0.1857, 0.1826, 0.1796, 0.1767, 0.1738, 0.1709,
    0.1681, 0.1654, 0.1627, 0.16, 0.1574, 0.1549, 0.1524, 0.1499, 0.1475, 0.1451, 0.1428, 0.1405, 0.1383,
    0.1361, 0.1339, 0.1318, 0.1297, 0.1277, 0.1257, 0.1237, 0.1218, 0.1198, 0.118, 0.1161, 0.1143, 0.1125,
    0.1108, 0.1091, 0.1074, 0.1057, 0.1041, 0.1025, 0.1009, 0.0994, 0.0978, 0.0963, 0.0949, 0.0934, 0.092,
    0.0906, 0.0892, 0.0878, 0.0865, 0.0852, 0.0839, 0.0826, 0.0814, 0.0802, 0.0789, 0.0778, 0.0766, 0.0754,
    0.0743, 0.0732, 0.0721, 0.071, 0.07, 0.0689, 0.0679, 0.0669, 0.0659, 0.0649, 0.0639, 0.063, 0.0621, 0.0611,
    0.0602, 0.0593, 0.0585, 0.0576, 0.0568, 0.0559, 0.0551, 0.0543, 0.0535, 0.0527, 0.0519, 0.0512, 0.0504, 0.0497,
    0.049, 0.0482, 0.0475, 0.0469, 0.0462, 0.0455, 0.0448, 0.0442, 0.0435, 0.0429, 0.0423, 0.0417, 0.0411, 0.0405,
    0.0399, 0.0393, 0.0388, 0.0382, 0.0376, 0.0371, 0.0366, 0.036, 0.0355, 0.035, 0.0345, 0.034, 0.0335, 0.033, 0.0326,
    0.0321, 0.0317, 0.0312, 0.0308, 0.0303, 0.0299, 0.0295, 0.029, 0.0286, 0.0282, 0.0278, 0.0274, 0.027, 0.0266, 0.0263,
    0.0259, 0.0255, 0.0252, 0.0248, 0.0245, 0.0241, 0.0238, 0.0234, 0.0231, 0.0228, 0.0225, 0.0221, 0.0218, 0.0215,
    0.0212, 0.0209, 0.0206, 0.0203, 0.0201, 0.0198, 0.0195, 0.0192, 0.019, 0.0187, 0.0184, 0.0182, 0.0179, 0.0177,
    0.0174, 0.0172, 0.0169, 0.0167, 0.0165, 0.0162, 0.016, 0.0158, 0.0156, 0.0153, 0.0151, 0.0149, 0.0147, 0.0145,
    0.0143, 0.0141, 0.0139, 0.0137, 0.0135, 0.0133, 0.0132, 0.013, 0.0128, 0.0126, 0.0124, 0.0123, 0.0121, 0.0119,
    0.0118, 0.0116, 0.0114, 0.0113, 0.0111, 0.011, 0.0108, 0.0107, 0.0105, 0.0104, 0.0102, 0.0101, 0.01, 0.0098,
    0.0097, 0.0096, 0.0094, 0.0093, 0.0092, 0.009, 0.0089, 0.0088, 0.0087, 0.0086, 0.0084, 0.0083, 0.0082, 0.0081,
    0.008, 0.0079, 0.0078, 0.0077, 0.0076, 0.0075, 0.0074, 0.0073, 0.0072, 0.0071, 0.007, 0.0069, 0.0068, 0.0067,
    0.0066, 0.0065, 0.0064, 0.0063, 0.0062, 0.0062, 0.0061, 0.006, 0.0059, 0.0058, 0.0057, 0.0057, 0.0056, 0.0055,
    0.0054, 0.0054, 0.0053, 0.0052, 0.0052, 0.0051, 0.005, 0.0049, 0.0049, 0.0048, 0.0047, 0.0047, 0.0046, 0.0046,
    0.0045, 0.0044, 0.0044, 0.0043, 0.0043, 0.0042, 0.0041, 0.0041, 0.004, 0.004, 0.0039, 0.0039, 0.0038, 0.0038,
    0.0037, 0.0037, 0.0036, 0.0036, 0.0035, 0.0035, 0.0034, 0.0034, 0.0033, 0.0033, 0.0032, 0.0032, 0.0032, 0.0031,
    0.0031, 0.003, 0.003, 0.003, 0.0029, 0.0029, 0.0028, 0.0028, 0.0028, 0.0027, 0.0027, 0.0026, 0.0026, 0.0026,
    0.0025, 0.0025, 0.0025, 0.0024, 0.0024, 0.0024, 0.0023, 0.0023, 0.0023, 0.0023, 0.0022, 0.0022, 0.0022, 0.0021,
    0.0021, 0.0021, 0.0021, 0.002, 0.002, 0.002, 0.0019
  )

  x0 <- seq(0, 3.99, by = 0.01)
  if (is.na(x)) {
    return(NA)
  }
  if (x >= 8) {
    return(0)
  }
  if (x > 3.99) {
    return(-0.0019 * x / 4.01 + 8 * 0.0019 / 4.01)
  }

  i <- sum(x0 < x)
  if (i == 0) i <- 1
  return((p0[i + 1] - p0[i]) * x / 0.01 + (x0[i + 1] * p0[i] - p0[i + 1] * x0[i]) / 0.01)
}
