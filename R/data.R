#' Communities and Crime Dataset
#'
#' A dataset containing socio-economic and law enforcement data for 1,994 communities
#' across the United States, used to predict violent crime rates.
#'
#' @format A numeric matrix with 1,994 rows and 100 columns. The last column
#' is the response variable (per capita violent crime rate).
#' \describe{
#'   \item{Y}{Response variable: per capita violent crime rate.}
#'   \item{X1-X99}{99 predictors describing demographic and law-enforcement characteristics.}
#' }
#' @source \url{https://archive.ics.uci.edu/dataset/183/communities+and+crime}
#' @usage data(data_crime)
"data_crime"


#' Acute Myeloid Leukemia (AML) RNA-Seq Data
#'
#' Cleaned RNA-Seq expression profiles of 444 patients with Acute Myeloid Leukemia,
#' providing two clinical response variables for testing.
#'
#' @format A list with two components:
#' \describe{
#'   \item{x}{A numeric matrix of 22,834 gene expression levels (ultra-high dimensional).}
#'   \item{y}{A data frame with 444 observations on 2 variables:
#'     \itemize{
#'       \item \code{eln_risk}: 1 for high-risk, 0 for non-high-risk (ELN2017 criteria).
#'       \item \code{de_novo}: 1 for de novo AML, 0 for secondary or other types.
#'     }
#'   }
#' }
#' @references Bottomly D, Long N, Schultz A R, et al. Integrative analysis of drug response and clinical outcome in acute myeloid leukemia[J]. Cancer cell, 2022, 40(8): 850-864. e9.
#' @source \url{https://www.cbioportal.org/study/summary?id=aml_ohsu_2022}
#' @usage data(data_AML)
"data_AML"
