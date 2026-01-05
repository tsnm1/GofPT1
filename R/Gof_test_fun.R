library(MASS)
library(glmnet)
library(LassoSIR)
library(VariableScreening)
library(psych)
library(caret)

#' @title A Two-Step Projection-Based Goodness-of-Fit Test for Ultra-High Dimensional Sparse Regressions
#'
#' @description The function implements a novel two-step strategy for testing the goodness-of-fit of para metric sparse regression models in ultra-high dimensional settings,
#' where the predictor dimension far exceeds the sample size. It first constructs multiple test statistics based on projected predictors from distinct projections and
#' second employs powerful $p$-value combination procedures, such as the minimum $p$-value and the Fisher combination of $p$-value,
#' to form our final tests and enhance power.
#' @param X Input matrix with \code{n} rows, each a \code{p}-dimensional observation vector.
#' @param Y Response vector.
#' @param fam Family type for GLM Models fitting. Must be "gaussian", "binomial" or "poisson".
#' @param penalize \code{TRUE} if penalization should be used in glmnet function when fitting the GLM models (see Details below).
#' #' @details
#'   This function tests if the conditional mean of \code{y} given \code{X} could be
#'   originating from a GLM family specified by the user via \code{fam}.
#'
#'   The function works by splitting the data into data1 and data2,
#'   and computes a GLM fit and choose some projections based "LassoSIR" on both parts.
#'   If \code{penalize == TRUE}, these fits use \code{cv.glmnet} from package \code{glmnet},
#'   otherwise they use \code{glmnet} with penalty set to 0.
#'   If \code{fam=="binomial"}, the param categorical in "LassoSIR" function should be TRUE,FALSE otherwise.
#'
#' @return the function outputs two p-values(the minimump-value and the Fisher combination of p-values)
#' @export
#'
#' @examples
#' # Simulation example  H11
#' set.seed(123)
#' n <- 300
#' p <- 50
#' pho <- 0.4
#' mu <- rep(0, p)
#' v <- pho^(0:(p - 1))
#' sigma <- toeplitz(v)
#' x <- mvrnorm(n, mu, sigma)
#' beta0 <- c(rep(1, 5), rep(0, p - 5))
#'
#' a <- 0
#' y <- x %*% beta0 + a * 0.1 * (x %*% beta0)^2 + rnorm(n)
#' output11 <- Gof_CPB_test(x, y, family = "gaussian")
#' a <- 1
#' y <- x %*% beta0 + a * 0.1 * (x %*% beta0)^2 + rnorm(n)
#' output13 <- Gof_CPB_test(x, y, family = "gaussian")
#'
#' # Real data analysis
#' dim(data_crime)
#' y <- data_crime[, 100]
#' x <- data_crime[, -100]
#' output15 <- Gof_CPB_test(x, y, family = "gaussian")
#'
#' # Simulation example # H21
#' set.seed(123)
#' n <- 300
#' p <- 50
#' pho <- 0.4
#' mu <- rep(0, p)
#' v <- pho^(0:(p - 1))
#' sigma <- toeplitz(v)
#' x <- mvrnorm(n, mu, sigma)
#' beta0 <- c(rep(1, 5), rep(0, p - 5))
#'
#' a <- 0
#' z <- x %*% beta0 + a * 0.2 * (x %*% beta0)^2
#' pr <- 1 / (1 + exp(-z))
#' y <- matrix(rbinom(n, 1, pr), ncol = 1)
#' output21 <- Gof_CPB_test(x, y, family = "binomial")
#' a <- 1
#' z <- x %*% beta0 + a * 0.2 * (x %*% beta0)^2
#' pr <- 1 / (1 + exp(-z))
#' y <- matrix(rbinom(n, 1, pr), ncol = 1)
#' output23 <- Gof_CPB_test(x, y, family = "binomial")
#'
#' # Real data analysis
#' dim(data_AML)
#' y <- data_AML$ELN_binary
#' x <- data_AML[, 1:(ncol(data_AML) - 2)]
#' x <- x[, -which(is.na(apply(data_AML, 2, sum)))]
#' output25 <- Gof_CPB_test(x, y, family = "binomial")
#'
Gof_CPB_test <- function(x, y, family = c("gaussian", "binomial", "poisson"),
                         penalize = TRUE, c_h = 1) {
  family <- match.arg(family)

  # Input Checks
  # if (!is.matrix(x) || ncol(x) < 1) stop("x should be a matrix with at least one column.")
  x <- tryCatch(as.matrix(x), error = function(e) stop("x must be a matrix or a data frame."))
  if (ncol(x) < 1) stop("X should have at least one column.")


  n <- nrow(x)
  p <- ncol(x)
  if (length(y) != n) stop("y must have nrow(x) components.")

  # Safety check strictly from original code
  if ((p >= n - 1) && (isFALSE(penalize))) {
    stop("When penalize=FALSE, you must have ncol(x) < nrow(x)-1. Try setting penalize=TRUE.")
  }
  y <- as.numeric(y)

  # Data Splitting
  # --- Data Splitting ---
  split <- split_data(x, y)
  n1 <- nrow(split$x1)
  n2 <- nrow(split$x2)

  # --- Step 1: Estimation ---
  # Using penalize parameter from top-level function
  est1 <- get_residuals_and_beta(split$x1, split$y1, family, penalize = penalize)
  est2 <- get_residuals_and_beta(split$x2, split$y2, family, penalize = penalize)

  # --- Step 2: SDR Projections ---
  sdr1 <- get_sdr_projections(split$x1, list(U = est1$U, original_y = split$y1), family, methods = "local")
  sdr2 <- get_sdr_projections(split$x2, list(U = est2$U, original_y = split$y2), family, methods = "local")

  # --- Step 3: Testing (Dual-parameter architecture) ---
  # Using U_pro and beta_pro as the dual projection targets
  # Cross-fitting: Test split 1 using projections from split 2
  h1 <- c_h * (n1^(-2 / 9))
  h2 <- c_h * (n2^(-2 / 9))
  p_pls1 <- calc_pls_pvals(split$x1, est1$U, sdr2$U_pro, sdr2$y_pro, h1)
  p_pls2 <- calc_pls_pvals(split$x2, est2$U, sdr1$U_pro, sdr1$y_pro, h2)

  # --- Step 4: Combination (Strictly following original logic) ---
  res1 <- combine_pvals(p_pls1)
  res2 <- combine_pvals(p_pls2)

  # Cauchy Combination of the aggregates
  cauchy_comb <- function(p) {
    p <- p[!is.na(p)]
    if (length(p) == 0) {
      return(NA)
    }
    # Truncation logic as requested
    p[p > 0.999999] <- 0.999999
    p[p < 1e-16] <- 1e-16
    1 - pcauchy(mean(tan((0.5 - p) * pi)))
  }

  pval_cauchy_fisher <- cauchy_comb(c(res1["fisher"], res2["fisher"]))
  pval_cauchy_min <- cauchy_comb(c(res1["min_p"], res2["min_p"]))

  return(list(
    pval_cauchy_fisher = pval_cauchy_fisher,
    pval_cauchy_min = pval_cauchy_min
  ))
}


#' @title "Asymptotic Distribution-Free Tests for Ultra-high Dimensional Parametric Regressions via Projected Empirical Processes and p-value Combination"
#'
#' @description This function constructs a Cramer-von Mises type test based on a martingale-transformed,
#' projected residual-marked empirical process. Furthermore, it proposes a novel hybrid test that aggregates empirical process-based tests and local smoothing tests using Cauchy combination.
#' @param X Input matrix with \code{n} rows, each a \code{p}-dimensional observation vector.
#' @param Y Response vector.
#' @param fam Family type for GLM model fitting. Must be "gaussian" or "binomial".
#'
#' @return A list containing various aggregated p-values, where \code{pval_cauchy_hybrid} is the
#' primary result representing the combination of both Martingale and PLS tests.
#'
#' @export
#'
#' @examples
#' # Simulation example # H11
#' set.seed(123)
#' n <- 300
#' p <- 50
#' pho <- 0.4
#' mu <- rep(0, p)
#' v <- pho^(0:(p - 1))
#' sigma <- toeplitz(v)
#' x <- mvrnorm(n, mu, sigma)
#' beta0 <- c(rep(1, 5), rep(0, p - 5))
#'
#' a <- 0
#' y <- x %*% beta0 + a * 0.1 * (x %*% beta0)^2 + rnorm(n)
#' output12 <- Gof_Pcvm_test(x, y, family = "gaussian")
#' a <- 1
#' y <- x %*% beta0 + a * 0.1 * (x %*% beta0)^2 + rnorm(n)
#' output14 <- Gof_Pcvm_test(x, y, family = "gaussian")
#'
#' # Real data analysis
#' dim(data_crime)
#' y <- data_crime[, 100]
#' x <- data_crime[, -100]
#' output16 <- Gof_Pcvm_test(x, y, family = "gaussian")
#'
#' # Simulation example  # H21
#' set.seed(123)
#' n <- 300
#' p <- 50
#' pho <- 0.4
#' mu <- rep(0, p)
#' v <- pho^(0:(p - 1))
#' sigma <- toeplitz(v)
#' x <- mvrnorm(n, mu, sigma)
#' beta0 <- c(rep(1, 5), rep(0, p - 5))
#'
#' a <- 0
#' z <- x %*% beta0 + a * 0.2 * (x %*% beta0)^2
#' pr <- 1 / (1 + exp(-z))
#' y <- matrix(rbinom(n, 1, pr), ncol = 1)
#' output22 <- Gof_Pcvm_test(x, y, family = "binomial")
#' a <- 1
#' z <- x %*% beta0 + a * 0.2 * (x %*% beta0)^2
#' pr <- 1 / (1 + exp(-z))
#' y <- matrix(rbinom(n, 1, pr), ncol = 1)
#' output24 <- Gof_Pcvm_test(x, y, family = "binomial")
#'
#' # Real data analysis
#' dim(data_AML)
#' y <- data_AML$ELN_binary
#' x <- data_AML[, 1:(ncol(data_AML) - 2)]
#' x <- x[, -which(is.na(apply(data_AML, 2, sum)))]
#' output26 <- Gof_Pcvm_test(x, y, family = "binomial")
#'
Gof_Pcvm_test <- function(x, y, family = c("gaussian", "binomial")) {
  family <- match.arg(family)

  # --- Data Splitting ---
  split <- split_data(x, y)

  # --- Step 1: Estimation ---
  # penalize=NULL triggers the Hybrid logic (Caret, GLM/Lasso switch)
  # return(list(U = as.matrix(U1), beta_pro = matrix(beta1_pro, ncol = 1),
  #             beta_hat = beta1_hat, intercept = intercept1, deri_link = deri))
  est1 <- get_residuals_and_beta(split$x1, split$y1, family, penalize = T)
  est2 <- get_residuals_and_beta(split$x2, split$y2, family, penalize = T)

  # --- Step 2: SDR Projections ---
  # return(list(
  #   U_pro = run_sdr(y_list$U, FALSE),
  #   y_pro = run_sdr(y_list$original_y, family == "binomial")
  # ))
  # return(list(U_pro = sir_Upro, y_pro = sir_ypro))
  sdr1 <- get_sdr_projections(split$x1, list(U = est1$U, original_y = split$y1), family, methods = "global")
  sdr2 <- get_sdr_projections(split$x2, list(U = est2$U, original_y = split$y2), family, methods = "global")

  # --- Step 4: Testing ---
  # Martingale Test (Pcvm)
  p_mart1 <- calc_martingale_pvals(split$x1, est1, sdr2$U_pro, est2$beta_pro, family)
  p_mart2 <- calc_martingale_pvals(split$x2, est2, sdr1$U_pro, est1$beta_pro, family)

  # PLS Test (Local Smoothing)
  h1 <- nrow(split$x1)^(-2 / 9)
  h2 <- nrow(split$x2)^(-2 / 9)
  p_pls1 <- calc_pls_pvals(split$x1, est1$U, sdr2$U_pro, est2$beta_pro, h1)
  p_pls2 <- calc_pls_pvals(split$x2, est2$U, sdr1$U_pro, est1$beta_pro, h2)

  # --- Step 5: Combination ---
  res_mart_1 <- combine_pvals(p_mart1)
  res_mart_2 <- combine_pvals(p_mart2)
  res_mart_all <- combine_pvals(c(p_mart1, p_mart2))

  # Final Hybrid Combination
  res_hybrid <- combine_pvals(c(p_mart1, p_mart2, p_pls1, p_pls2))

  result <- list(
    # pval_cauchy1_Pcvm     = res_mart_1[["cauchy"]],
    # pval_cauchy2_Pcvm     = res_mart_2[["cauchy"]],
    pval_cauchy_Pcvm      = res_mart_all[["cauchy"]],
    pval_cauchy_hybrid    = res_hybrid[["cauchy"]]
  )
  return(result)
}
