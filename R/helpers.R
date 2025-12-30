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
get_residuals_and_beta <- function(x, y, family = c("gaussian", "binomial"),
                                   penalize = NULL) {
  family <- match.arg(family)
  x <- as.matrix(x)
  y <- as.numeric(y)
  n <- nrow(x)
  p <- ncol(x)

  if (is.null(penalize)) {
    penalize <- (p >= n - 1)
  }

  if (isTRUE(penalize)) {
    # --- Step 1: Lasso Selection ---
    lasso_model1 <- glmnet::cv.glmnet(x, y, family = family, intercept = TRUE)
    lasso_beta1 <- coef(lasso_model1, s = "lambda.min")[-1]
    index_beta1_non0 <- which(as.numeric(lasso_beta1) != 0)

    if (length(index_beta1_non0) == 0) {
      # --- Case: No variables selected ---
      intercept1 <- mean(y)
      U1 <- y - mean(y)
      beta1_hat <- lasso_beta1
      beta1_pro <- rep(1, p) / sqrt(p)
    } else {
      # --- Case: Variables selected -> Refit ---

      if (family == "gaussian") {
        # --- Gaussian Logic: GLM + Caret Cleanup ---                                            # second estimation
        x1_sec <- x[, index_beta1_non0, drop = FALSE]
        sec_model1 <- glm(y ~ x1_sec, family = gaussian)
        sec_beta1 <- unname(sec_model1$coefficients)[-1]
        lasso_beta1[index_beta1_non0] <- sec_beta1
        beta1_hat <- lasso_beta1
        intercept1 <- sec_model1$coefficients[1]
        pred1 <- predict(sec_model1, newx = x1_sec, type = "response")
        pred1 <- matrix(unname(pred1), ncol = 1)
        U1 <- y - pred1
        beta1_pro <- beta1_hat / sqrt(sum(beta1_hat^2))


        if (any(is.na(sec_beta1))) {
          library(caret)
          lasso_beta1 <- coef(lasso_model1, s = "lambda.min")[-1]

          cor_matrix <- cor(x1_sec, use = "pairwise.complete.obs")
          high_cor_vars <- findCorrelation(cor_matrix, cutoff = 0.9, verbose = TRUE)
          index_beta1_non0 <- index_beta1_non0[-high_cor_vars]
          x1_sec <- x[, index_beta1_non0]

          sec_model1 <- glm(y ~ x1_sec, family = gaussian)
          sec_beta1 <- unname(sec_model1$coefficients)[-1]
          lasso_beta1[index_beta1_non0] <- sec_beta1
          beta1_hat <- lasso_beta1
          pred1 <- predict(sec_model1, newx = x1_sec, type = "response")
          pred1 <- matrix(unname(pred1), ncol = 1)
          U1 <- y - pred1 # residual based on x2 y2
          # U1 <- sec_model1$residuals
          beta1_pro <- beta1_hat / sqrt(sum(beta1_hat^2))
        }
      } else {
        # --- Binomial Logic: GLM + cv.glmnet fallback ---
        x1_sec <- x[, index_beta1_non0, drop = FALSE]
        sec_model1 <- tryCatch(glm(y ~ x1_sec, family = binomial(link = "logit")), error = function(e) NULL)
        lasso_beta1[index_beta1_non0] <- unname(sec_model1$coefficients)[-1]
        beta1_hat <- lasso_beta1
        intercept1 <- unname(sec_model1$coefficients)[1]
        pred1 <- predict(sec_model1, newx = x1_sec, type = "response")
        pred1 <- matrix(unname(pred1), ncol = 1)
        U1 <- y - pred1
        beta1_pro <- beta1_hat

        if (is.null(sec_model1) || any(is.na(coef(sec_model1)))) {
          lasso_model1_1 <- cv.glmnet(x1_sec, y, family = "binomial", alpha = 1)
          lasso_beta1[index_beta1_non0] <- coef(lasso_model1_1, s = "lambda.min")[-1]
          beta1_hat <- lasso_beta1
          intercept1 <- coef(lasso_model1_1, s = "lambda.min")[1]
          pred1 <- predict(lasso_model1_1, newx = x1_sec, s = "lambda.min", type = "response")
          U1 <- y - as.vector(pred1)
          beta1_pro <- beta1_hat
        }
      }
    }
  } else {
    # --- Case: Non-penalize (Full GLM) ---
    fit1 <- glm(y ~ x, family = if (family == "gaussian") gaussian() else binomial(link = "logit"))
    beta1_hat <- coef(fit1)[-1]
    intercept1 <- coef(fit1)[1]
    U1 <- y - predict(fit1, type = "response")
    beta1_pro <- beta1_hat / sqrt(sum(beta1_hat^2))
  }

  # Derivative for Martingale (Strictly original binomial formula)
  deri <- if (family == "binomial") {
    # p_val <- as.vector(y - U1)
    # p_val * (1 - p_val)
    plogis(x %*% beta1_hat + intercept1) * (1 - plogis(x %*% beta1_hat + intercept1))
  } else {
    NULL
  }

  return(list(
    U = as.matrix(U1), beta_pro = matrix(beta1_pro, ncol = 1),
    beta_hat = beta1_hat, intercept = intercept1, deri_link = deri
  ))
}

# --- 3.3 SDR Projections ---
get_sdr_projections <- function(x, y_list, family, methods = c("local", "global")) {
  n <- nrow(x)
  p <- ncol(x)
  screen_num <- floor(n / log(n))

  # --- 1. sir_Upro (Residuals): Always H=10, categorical=FALSE ---
  sir_Upro <- tryCatch(
    {
      U <- y_list$U
      if (p <= screen_num) {
        fit <- LassoSIR::LassoSIR(x, U,
          H = 10, choosing.d = "automatic",
          solution.path = FALSE, categorical = FALSE,
          nfolds = 5, screening = FALSE
        )
        beta <- fit$beta
        sir_Upro <- beta / sqrt(colSums(beta^2))
      } else {
        # Screening for U always uses DC-SIS
        rank_U <- VariableScreening::screenIID(X = x, Y = U, method = "DC-SIS")
        idx_U <- seq(1:p)[rank_U$rank <= screen_num]
        fit <- LassoSIR::LassoSIR(x[, idx_U, drop = FALSE], U,
          H = 10, choosing.d = "automatic",
          solution.path = FALSE, categorical = FALSE,
          nfolds = 5, screening = FALSE
        )
        beta_f <- matrix(0, nrow = p, ncol = ncol(fit$beta))
        beta_f[idx_U, ] <- fit$beta
        sir_Upro <- beta_f / sqrt(colSums(beta_f^2))
      }
    },
    error = function(e) {
      return(NA)
    }
  )

  # --- 2. sir_ypro (Response): H and categorical depend on family ---
  if (methods == "local") {
    sir_ypro <- tryCatch(
      {
        Y <- as.numeric(y_list$original_y)

        # Logic for categorical and parameters based on family
        if (family == "binomial") {
          categorical_Y <- TRUE
          H_val <- 2
          met <- "MV-SIS"
        } else {
          categorical_Y <- FALSE
          H_val <- 10
          met <- "DC-SIS"
        }

        if (p <= screen_num) {
          fit <- LassoSIR::LassoSIR(x, Y,
            H = H_val, choosing.d = "automatic",
            solution.path = FALSE, categorical = categorical_Y,
            nfolds = 5, screening = categorical_Y
          )
          beta <- fit$beta
          sir_ypro <- beta / sqrt(colSums(beta^2))
        } else {
          # Screening method matches the data type
          rank_y <- VariableScreening::screenIID(X = x, Y = Y, method = met)
          idx_y <- seq(1:p)[rank_y$rank <= screen_num]
          fit <- LassoSIR::LassoSIR(x[, idx_y, drop = FALSE], Y,
            H = H_val, choosing.d = "automatic",
            solution.path = FALSE, categorical = categorical_Y,
            nfolds = 5, screening = categorical_Y
          )
          beta_f <- matrix(0, nrow = p, ncol = ncol(fit$beta))
          beta_f[idx_y, ] <- fit$beta
          sir_ypro <- beta_f / sqrt(colSums(beta_f^2))
        }
      },
      error = function(e) {
        return(NULL)
      }
    )
  } else {
    sir_ypro <- NULL
  }


  return(list(U_pro = sir_Upro, y_pro = sir_ypro))
}

# --- 3.4 Pcvm Test Statistics ---
calc_martingale_pvals <- function(x_source, est, sir_pro_target, beta_pro_target, family) {
  n <- nrow(x_source)
  U <- as.vector(est$U)

  # clean <- function(m) m[, !apply(m, 2, function(v) any(is.na(v))), drop = FALSE]
  sir_pro <- sir_pro_target
  beta_pro <- beta_pro_target

  # --- Branch A: Gaussian Case (Linear) ---
  if (family == "gaussian") {
    all_pro <- cbind(sir_pro, beta_pro)
    n_pro <- ncol(all_pro)
    pvals <- numeric(n_pro)

    for (q in 1:n_pro) {
      x_pro <- x_source %*% all_pro[, q]
      x_Indictor <- ifelse(x_pro %*% rep(1, n) <= rep(1, n) %*% t(x_pro), 1, 0)

      hat_A <- rbind(t(x_pro), rep(1, n))
      Gamma_inv <- array(0, dim = c(2, 2, n))
      for (i in 1:n) {
        Gamma_inv[, , i] <- MASS::ginv((hat_A %*% diag(x_Indictor[i, ])) %*% t(hat_A))
      }
      Integral <- hat_A %*% diag(as.vector(U)) %*% t(x_Indictor)

      martingle_sec <- diag(0, n)
      for (l in 1:n) {
        martingle_sec[l, ] <- (t(hat_A[, l]) %*% Gamma_inv[, , l] %*% Integral[, l]) %*% x_Indictor[l, ]
      }

      martingle_sta <- (1 / sqrt(n)) * t(U) %*% x_Indictor - (1 / sqrt(n)) * colSums(martingle_sec)

      # Denominator logic based on original linear code
      t_val <- sort(x_pro)[floor(0.99 * n)]
      sigma_square <- mean(U^2)
      F_val <- (mean(x_pro <= t_val))^2
      PCvM <- (1 / (sigma_square * F_val)) * mean((x_pro <= t_val) * (t(martingle_sta)^2))
      pvals[q] <- pvalue_integ_Brown(PCvM)
    }
  }

  # --- Branch B: Binomial Case (Logit) ---
  else {
    n_sir <- ncol(sir_pro)
    pvals <- numeric(n_sir + 1)
    deri_link <- est$deri # Expected to be p*(1-p)

    # 1. Loop for sir_Upro (Non-parametric 3x3 Branch)
    for (q in 1:n_sir) {
      x_pro <- x_source %*% sir_pro[, q]
      band_y <- cbind(U^2, deri_link, (x_source %*% est$beta_hat) * deri_link)
      h <- bandwidth_choice(x_pro, band_y)

      Ker_inter <- (x_pro %*% rep(1, n) - rep(1, n) %*% t(x_pro)) / h
      Ker_indictor <- ifelse(abs(Ker_inter) <= 1, 1, 0)
      kernel <- (3 / 4) * (1 - Ker_inter^2) * Ker_indictor

      Indictor <- ifelse(x_pro %*% rep(1, n) <= rep(1, n) %*% t(x_pro), 1, 0)

      # Construct hat_A components
      denom <- kernel %*% U^2 + n^(-12)
      A_1 <- ((kernel %*% deri_link) * x_pro) / denom
      A_2 <- (kernel %*% ((x_source %*% est$beta_hat) * deri_link)) / denom
      A_3 <- (kernel %*% deri_link) / denom
      hat_A <- rbind(t(A_1), t(A_2), t(A_3))

      Gamma_inv <- array(0, dim = c(3, 3, n))
      for (i in 1:n) {
        Gamma_inv[, , i] <- MASS::ginv(((rep(1, 3) %*% t(U)^2) * hat_A) %*% (t(hat_A) * (as.matrix(Indictor[i, ]) %*% rep(1, 3))) + n^(-12))
      }
      Integral <- ((rep(1, 3) %*% t(U)) * hat_A) %*% t(Indictor)

      martingle_sec <- matrix(0, n, n)
      for (l in 1:n) {
        martingle_sec[l, ] <- (U[l]^2 * t(hat_A[, l]) %*% Gamma_inv[, , l] %*% Integral[, l]) %*% Indictor[l, ]
      }
      martingle_sta <- (1 / sqrt(n)) * t(U) %*% Indictor - (1 / sqrt(n)) * colSums(martingle_sec)

      # Denominator logic based on original logit code
      t_val <- sort(x_pro)[floor(0.99 * n)]
      psi <- mean(U^2 * (x_pro <= t_val))
      PCvM <- (1 / psi^2) * mean(U^2 * (x_pro <= t_val) * (t(martingle_sta)^2))

      pvals[q] <- pvalue_integ_Brown(PCvM)
    }

    # 2. Calculation for beta_pro (Parametric 2x2 Branch)
    x_pro_beta <- x_source %*% beta_pro[, 1]
    beta_Indictor <- ifelse(x_pro_beta %*% rep(1, n) <= rep(1, n) %*% t(x_pro_beta), 1, 0)

    hat_beta_A <- rbind(t(x_pro_beta), rep(1, n))
    Gamma_beta_inv <- array(0, dim = c(2, 2, n))
    for (i in 1:n) {
      Gamma_beta_inv[, , i] <- MASS::ginv(((rep(1, 2) %*% t(U)^2) * hat_beta_A) %*% (t(hat_beta_A) * (as.matrix(beta_Indictor[i, ]) %*% rep(1, 2))) + n^(-10))
    }

    Integral_beta <- ((rep(1, 2) %*% t(U)) * hat_beta_A) %*% t(beta_Indictor)

    martingle_beta_sec <- matrix(0, n, n)
    for (l in 1:n) {
      martingle_beta_sec[l, ] <- (U[l]^2 * t(hat_beta_A[, l]) %*% Gamma_beta_inv[, , l] %*% Integral_beta[, l]) %*% beta_Indictor[l, ]
    }
    martingle_beta_sta <- (1 / sqrt(n)) * t(U) %*% beta_Indictor - (1 / sqrt(n)) * colSums(martingle_beta_sec)

    t_beta <- sort(x_pro_beta)[floor(0.99 * n)]
    psi_beta <- (1 / n) * sum(U^2 * (x_pro_beta <= t_beta))
    PCvM_beta <- (1 / psi_beta^2) * mean(U^2 * (x_pro_beta <= t_beta) * (t(martingle_beta_sta)^2))

    pvals[n_sir + 1] <- pvalue_integ_Brown(PCvM_beta)
  }

  return(pvals)
}

# --- 3.5 PLS Test Statistics ---
calc_pls_pvals <- function(x_source, U, sir_pro_target, beta_pro_target, h) {
  # --- Step 1: Matrix Construction (Matching original logic) ---
  # Combine SDR directions and Beta direction into a single projection set
  pro_pls <- cbind(sir_pro_target, beta_pro_target)

  # Clean NA columns as per tryCatch logic in original snippets
  pro_pls <- pro_pls[, !apply(pro_pls, 2, function(v) any(is.na(v))), drop = FALSE]

  n <- nrow(x_source)
  pro_pls_num <- ncol(pro_pls)
  pval_matrix_PLS <- matrix(nrow = 1, ncol = pro_pls_num)

  # errormat: Residual matrix based on current split
  errormat <- U %*% t(U)
  epsilon <- 1e-12 # Matching the epsilon in your PLS snippet

  # --- Step 2: Loop through each projection ---
  for (q in 1:pro_pls_num) {
    x_pro <- x_source %*% pro_pls[, q]

    # Kernel function matrix (Epanechnikov kernel)
    x_pro_mat <- ((x_pro) %*% matrix(1, 1, n) - matrix(1, n, 1) %*% (t(x_pro))) / h
    indictor <- ifelse(abs(x_pro_mat) <= 1, 1, 0)
    kermat <- (3 / 4) * (1 - x_pro_mat^2) * indictor

    KE <- kermat * errormat
    num <- sum(KE) - psych::tr(KE)
    # den <- sqrt(2 * (sum(KE^2) - psych::tr(KE^2)) + epsilon)
    den <- sqrt(2 * (sum(KE^2) - psych::tr(KE^2)))

    Tn <- num / den

    pval_matrix_PLS[, q] <- 1 - pnorm(Tn)
  }

  # Return as a vector to be consistent with p_mart
  return(as.vector(pval_matrix_PLS))
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

# --- 3.7 The choice for the bandwidth for martingale-based test ---
bandwidth_choice <- function(x, y_mat) {
  n <- length(x)
  p <- ncol(x)
  c_h <- seq(0.25, 1.25, 0.15)
  len_h <- length(c_h)
  h <- c_h * n^(-2 / 9)
  CV <- rep(0, len_h)

  for (i in 1:len_h) {
    Ker_inter <- (x %*% rep(1, n) - rep(1, n) %*% t(x)) / h[i] # kernel function matrix
    # kernel <- exp(-0.5*Ker_inter^2)                                          # Gaussian kernel

    Ker_indictor <- ifelse(abs(Ker_inter) <= 1, 1, 0)
    kernel <- (3 / 4) * (1 - Ker_inter^2) * Ker_indictor

    diag(kernel) <- rep(0, n)

    R1 <- mean((y_mat[, 1] - (kernel %*% y_mat[, 1]) / (rowSums(kernel) + n^(-12)))^2)
    R2 <- mean((y_mat[, 2] - (kernel %*% y_mat[, 2]) / (rowSums(kernel) + n^(-12)))^2)
    R3 <- mean((y_mat[, 3] - (kernel %*% y_mat[, 3]) / (rowSums(kernel) + n^(-12)))^2)
    CV[i] <- R1 + R2 + R3
  }
  index <- which.min(CV)
  h0 <- c_h[index] * n^(-2 / 9) # The bandwidth for the martingale with beta1 and alpha as projections
  return(c(h0))
}

# --- 3.8 The p-value of int_0^1 B(t)^2 dt where B(t) is the standard brown motion ---
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
    p <- NA
  } else if (x > 3.99 & x < 8) {
    p <- -0.0019 * x / 4.01 + 8 * 0.0019 / 4.01
  } else if (x >= 8) {
    p <- 0
  } else if (x == 0) {
    p <- 1
  } else {
    i <- sum(x0 < x)
    p <- (p0[i + 1] - p0[i]) * x / (x0[i + 1] - x0[i]) + (x0[i + 1] * p0[i] - p0[i + 1] * x0[i]) / (x0[i + 1] - x0[i])
  }
  return(p)
}
