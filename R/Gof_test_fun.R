library(MASS)
library(glmnet)
library(LassoSIR)
library(VariableScreening)
library(psych)
library(caret)

#' @title A Two-Step Projection-Based Goodness-of-Fit Test for Ultra-High Dimensional Sparse Regressions
#'
#' The function can test good
#' @param X Input matrix with \code{n} rows, each a \code{p}-dimensional observation vector.
#' @param Y Response vector.
#' @param fam Family type for GLM Models fitting. Must be "gaussian", "binomial" or "poisson".
#' @param penalize \code{TRUE} if penalization should be used in glmnet function when fitting the GLM models (see Details below).
#'
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
#' set.seed(22)
#' X <- matrix(rnorm(300 * 30), 300, 30)
#' z <- X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]
#' pr <- 1 / (1 + exp(-z))
#' Y <- rbinom(nrow(X), 1, pr)
#' output <- Gof_CPB_test(X, Y, fam = "binomial", TRUE)
#' output
Gof_CPB_test <- function(X, Y, fam = c("gaussian", "binomial", "poisson"), penalize = TRUE) {
  fam <- match.arg(fam)

  # Input Checks
  if (!is.matrix(X) || ncol(X) < 1) stop("X should be a matrix with at least one column.")
  n <- nrow(X)
  p <- ncol(X)
  if (length(Y) != n) stop("Y must have nrow(X) components.")

  # Safety check strictly from original code
  if ((p >= n - 1) && (isFALSE(penalize))) {
    stop("When penalize=FALSE, you must have ncol(X) < nrow(X)-1. Try setting penalize=TRUE.")
  }
  Y <- as.numeric(Y)

  # Data Splitting
  split <- split_data(X, Y)
  x1 <- split$x1
  y1 <- split$y1
  x2 <- split$x2
  y2 <- split$y2

  # --- Step 1: Estimation ---
  # CPB uses strict GLMNET logic via the helper
  est1 <- get_residuals_and_beta(x1, y1, fam, penalize = penalize)
  est2 <- get_residuals_and_beta(x2, y2, fam, penalize = penalize)

  # --- Step 2: SDR Projections ---
  sdr1 <- get_sdr_projections(x1, list(U = est1$U, original_y = y1), fam)
  sdr2 <- get_sdr_projections(x2, list(U = est2$U, original_y = y2), fam)

  # --- Step 3: Projection Sets (CPB uses U and Y projections) ---
  pro1 <- cbind(sdr1$U_pro, sdr1$y_pro)
  pro2 <- cbind(sdr2$U_pro, sdr2$y_pro)

  # Clean NA columns
  pro1 <- pro1[, !apply(pro1, 2, function(v) any(is.na(v))), drop = FALSE]
  pro2 <- pro2[, !apply(pro2, 2, function(v) any(is.na(v))), drop = FALSE]

  # --- Step 4: PLS Testing ---
  p_pls1 <- calc_pls_pvals(x1, est1$U, pro2)
  p_pls2 <- calc_pls_pvals(x2, est2$U, pro1)

  # --- Step 5: Combination ---
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
#' @description This function implements a novel hybrid goodness-of-fit test for sparse parametric regression models
#' in ultra-high dimensional settings. It aggregates p-values from two distinct approaches—
#' martingale-transformed empirical process-based tests and local smoothing tests—using the Cauchy
#' combination method to ensure robust power against various alternatives.
#'
#' @param X Input matrix with \code{n} rows, each a \code{p}-dimensional observation vector.
#' @param Y Response vector.
#' @param fam Family type for GLM model fitting. Must be "gaussian" or "binomial".
#'
#' @return A list containing various aggregated p-values, where \code{pval_cauchy_hybrid} is the
#' primary result representing the combination of both Martingale and PLS tests.
#'
#' @export
hybrid_test <- function(x, y, family = c("gaussian", "binomial")) {
  family <- match.arg(family)

  # --- Data Splitting ---
  split <- split_data(x, y)
  x1 <- split$x1
  y1 <- split$y1
  x2 <- split$x2
  y2 <- split$y2

  # --- Step 1: Estimation ---
  # penalize=NULL triggers the Hybrid logic (Caret, GLM/Lasso switch)
  est1 <- get_residuals_and_beta(x1, y1, family, penalize = NULL)
  est2 <- get_residuals_and_beta(x2, y2, family, penalize = NULL)

  # --- Step 2: SDR Projections ---
  sdr1 <- get_sdr_projections(x1, list(U = est1$U, original_y = y1), family)
  sdr2 <- get_sdr_projections(x2, list(U = est2$U, original_y = y2), family)

  # --- Step 3: Projection Sets ---
  # Hybrid uses U projections and the Beta direction
  pro1 <- cbind(sdr1$U_pro, est1$beta_pro)
  pro2 <- cbind(sdr2$U_pro, est2$beta_pro)

  # Clean NA columns
  pro1 <- pro1[, !apply(pro1, 2, function(v) any(is.na(v))), drop = FALSE]
  pro2 <- pro2[, !apply(pro2, 2, function(v) any(is.na(v))), drop = FALSE]

  # --- Step 4: Testing ---
  # Martingale Test
  p_mart1 <- calc_martingale_pvals(x1, est1, pro2, family)
  p_mart2 <- calc_martingale_pvals(x2, est2, pro1, family)

  # PLS Test
  p_pls1 <- calc_pls_pvals(x1, est1$U, pro2)
  p_pls2 <- calc_pls_pvals(x2, est2$U, pro1)

  # --- Step 5: Combination ---
  res_mart_1 <- combine_pvals(p_mart1)
  res_mart_2 <- combine_pvals(p_mart2)
  res_mart_all <- combine_pvals(c(p_mart1, p_mart2))

  res_pls_1 <- combine_pvals(p_pls1)
  res_pls_2 <- combine_pvals(p_pls2)
  res_pls_all <- combine_pvals(c(p_pls1, p_pls2))

  # Final Hybrid Combination
  res_hybrid <- combine_pvals(c(p_mart1, p_mart2, p_pls1, p_pls2))

  result <- list(
    pval_cauchy1_Pcvm     = res_mart_1[["cauchy"]],
    pval_cauchy2_Pcvm     = res_mart_2[["cauchy"]],
    pval_cauchy_Pcvm      = res_mart_all[["cauchy"]],
    pval_fisher1_PLS      = res_pls_1[["fisher"]],
    pval_fisher2_PLS      = res_pls_2[["fisher"]],
    pval_fisher_cauchy    = res_pls_all[["cauchy"]],
    pval_cauchy_hybrid    = res_hybrid[["cauchy"]]
  )
  return(result)
}
