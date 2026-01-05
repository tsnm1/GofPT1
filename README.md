
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GofPT1

<!-- badges: start -->

<!-- badges: end -->

This R package implements the methodologies developed in the following two works:

* **A Two-Step Projection-Based Goodness-of-Fit Test for Ultra-High Dimensional Sparse Regressions**: 
This paper proposes a two-step strategy that first constructs multiple test statistics based on projected predictors from distinct projections to mitigate the dimensionality problem. 
The second step employs p-value combination methods (e.g., the minimum p-value and the Fisher combination) to form the final test statistics. 
This projection-based approach significantly mitigates the dimensionality problem, enabling the tests to detect local alternatives converging to the null at the rate as if the predictor were univariate.

* **Asymptotic Distribution-Free Tests for Ultra-high Dimensional Parametric Regressions via Projected Empirical Processes and p-value Combination**:
This paper constructs a Cramer-von Mises type test based on a martingale-transformed, 
projected residual-marked empirical process, which renders the test asymptotically distribution-free. 
Furthermore, it proposes a novel hybrid test that aggregates empirical process-based and local smoothing tests using Cauchy combination, 
which is powerful against both low-frequency and high-frequency alternatives.


## Installation

You can install the development version of GofPT1 like so:

``` r
# install.packages("devtools")
devtools::install_github("tsnm1/GofPT1")
```

## Functions Overview

- **Gof_CPB_test(X, Y, fam, penalize)**: A Two-Step Projection-Based
  Goodness-of-Fit Test for Ultra-High Dimensional Sparse Regressions

- **Gof_Pcvm_test(x, y, family)**: Asymptotic Distribution-Free Tests for
  Ultra-high Dimensional Parametric Regressions via Projected Empirical
  Processes and p-value Combination

## Example

This section provides basic examples demonstrating how to use *GofPT1* to test the goodness-of-fit for high-dimensional linear and logistic regression models.

First, load the necessary libraries:

``` r
library(GofPT1)
library(glmnet) # Required for estimation
#> Loading required package: Matrix
#> Loaded glmnet 4.1-10
library(MASS) # Required for matrix operations
```

1. Linear Regression Models 

In this section, we apply the tests to continuous response variables.

``` r

# Simulation example 
# We simulate a high-dimensional dataset with n=300 and p=50.

set.seed(123)
n <- 300; p <- 50; pho <- 0.4
mu <- rep(0, p);v <- pho^(0:(p-1))
sigma <- toeplitz(v)
x <- mvrnorm(n, mu, sigma)
beta0 <- c(rep(1,5),rep(0,p-5))

# Case 1 H11 (Null Hypothesis H0): The true model is linear. 
a <- 0
y <- x %*% beta0 + a * 0.1*(x %*% beta0)^2 + rnorm(n)   
output11 <- Gof_CPB_test(x, y, fam = "gaussian")
output12 <- Gof_Pcvm_test(x, y, family = "gaussian")

# Case 2 (Alternative Hypothesis H1): We introduce a quadratic term to the response generation.
a <- 1
y <- x %*% beta0 + a * 0.1*(x %*% beta0)^2 + rnorm(n)   
output13 <- Gof_CPB_test(x, y, fam = "gaussian")
output14 <- Gof_Pcvm_test(x, y, family = "gaussian")


# Real data analysis: Crime Data 
dim(data_crime)
#> [1] 1994  100
y <- data_crime[, 100]
x <- data_crime[, -100]
output15 <- Gof_CPB_test(x, y, fam = "gaussian")
output16 <- Gof_Pcvm_test(x, y, fam = "gaussian")
```

2. Logistic Regression Models 

Here, we apply the tests to binary response variables.

``` r
# Simulation example 
# Similar to the linear case, we simulate data under both the null hypothesis (standard logistic model) and the alternative hypothesis (logistic model with quadratic terms).
set.seed(123)
n <- 300; p <- 50; pho <- 0.4
mu <- rep(0, p);v <- pho^(0:(p-1))
sigma <- toeplitz(v)
x <- mvrnorm(n, mu, sigma)
beta0 <- c(rep(1,5),rep(0,p-5))

# H21
a <- 0
z <- x %*% beta0 + a*0.2*(x %*% beta0)^2
pr <- 1/(1 + exp(-z)) 
y <- matrix(rbinom(n, 1, pr),ncol = 1)
output21 <- Gof_CPB_test(x, y, fam = "binomial")
output22 <- Gof_Pcvm_test(x, y, fam = "binomial")

a <- 1
z <- x %*% beta0 + a*0.2*(x %*% beta0)^2
pr <- 1/(1 + exp(-z)) 
y <- matrix(rbinom(n, 1, pr),ncol = 1)
output23 <- Gof_CPB_test(x, y, fam = "binomial")
output24 <- Gof_Pcvm_test(x, y, fam = "binomial")

# Real data analysis: AML Data
dim(data_AML)
#> [1]   444 22845
y <- data_AML$ELN_binary
x <- data_AML[,1:(ncol(data_AML)-2)]
x <- x[,-which(is.na(apply(data_AML,2,sum)))]
output25 <- Gof_CPB_test(x, y, fam = "binomial")
output26 <- Gof_Pcvm_test(x, y, fam = "binomial")
```

## References

[1] [TAN F, LIU J, PENG H, ZHU L. A Two-Step Projection-Based Goodness-of-Fit Test for Ultra-High Dimensional Sparse Regressions[EB/OL]. 2025. arXiv:2412.10721..](https://arxiv.org/abs/2412.10721v3)

[2] [TAN F, TANG S, ZHU L. Asymptotic Distribution-Free Tests for Ultra-high Dimensional Parametric Regressions via Projected Empirical Processes and p-value Combination[EB/OL]. 2026. arXiv:2601.00541..](https://arxiv.org/abs/2601.00541v1)
