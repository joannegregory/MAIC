#' Summarize baseline characteristics before and after weighting
#'
#' Summarize baseline characteristics before and after weighting
#'
#' @param analysis_data  A data frame containing data with propensity weights
#'   (derived from \code{\link{analysis_dataset}}).
#' @param matching_vars A character vector giving the names of the covariates used
#'   in matching.
#'
#' @return A table with three rows summarising the baseline characteristics for
#'   weighted intervention data, unweighted intervention data and comparator
#'   data.
#'
#' @seealso \code{\link{analysis_dataset}}
#' @export
summarize_baselines <- function(analysis_data, matching_vars){
  # summarize binary variables as n/N (%)

  # summarize continuous variables as mean (SD)

  # summarize number of patients/ ESS

  # Combine above

  # Return a table of baseline characteristics across the levels:  intervention
  # weighted, intervention unweighted and comparator
}

#' Summarize ratio of variance of baseline characteristics before and after weighting
#'
#' Summarize ratio of variance of baseline characteristics before and after
#' weighting. The purpose of this metric is to assess whether re-weighting has
#' affected the shape of the distribution as well as the location of the
#' distribution for each matching variable.
#'
#' @param analysis_data  A data frame containing data with weights (derived from
#'   \code{\link{analysis_dataset}}).
#' @param matching_vars A character vector giving the names of the covariates used
#'   in matching.
#'
#' @return Table of the ratio of variance before and after matching
#'
#' @seealso \code{\link{analysis_dataset}}
#' @export
summarize_baseline_variance <- function(analysis_data, matching_vars){
  # calculate variance of continuous variables before and after matching (with
  # and without weights)

  # calculate variance of binary variables before and after matching (with and
  # without weights)

  # calculate the ratio of the variance

  # Return a table of the ratio of variance before and after matching
}
