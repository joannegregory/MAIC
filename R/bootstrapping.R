#' Bootstrapping for MAIC weighted hazard ratios
#'
#' A function required for the "statistic" argument in the \code{\link{boot}} function.
#' Performs MAIC weighting using {\link{estimate_weights}} and returns a weighted hazard ratio (HR) from a Cox proportional hazards model.
#' @param intervention_data  A data frame containing containing individual patient data from the intervention study.
#' @param comparator_data A data frame containing pseudo individual patient data from the comparator study.
#'  The outcome variables names must match intervention_data.
#' @param matching A character vector giving the names of the covariates to use
#'   in matching. These names must match the column names in intervention_data
#'   and comparator_data.
#' @param i Index used to select a sample within \code{\link{boot}}.
#' @param model A model formula in the form 'Surv(Time, Event==1) ~ ARM'.
#'   Variable names need to match the corresponding columns in intervention_data.
#'
#' @return The HR as a numeric value.
#'
#' @seealso \code{\link{estimate_weights}}, \code{\link{boot}}
#' @export

bootstrap_HR <- function(intervention_data, comparator_data, matching, i, model){

  # Samples the data
  bootstrap_data <- intervention_data[i,]

  # Estimates weights
  perform_wt <- estimate_weights(intervention_data=bootstrap_data, matching_vars=matching, comparator_data=comparator_data)

  # survival data stat
  cox_model <- survival::coxph(model, data = perform_wt$analysis_data, weights = wt)
  HR <- exp(cox_model$coefficients)
}



#' Bootstrapping for MAIC weighted odds ratios
#'
#' A function required for the "statistic" argument in the \code{\link{boot}} function.
#' Performs MAIC weighting using {\link{estimate_weights}} and returns a weighted odds ratio (OR) from a binomial generalised linear model.
#'
#' @param intervention_data  A data frame containing containing individual patient data from the intervention study.
#' @param comparator_data A data frame containing pseudo individual patient data from the comparator study.
#'  The outcome variables names must match intervention_data.
#' @param matching A character vector giving the names of the covariates to use
#'   in matching. These names must match the column names in intervention_data
#'   and comparator_data.
#' @param i Index used to select a sample within \code{\link{boot}}.
#' @param model A model formula in the form 'endpoint ~ treatment_var'.
#'   Variable names need to match the corresponding columns in intervention_data.
#' @return The OR as a numeric value.
#'
#' @seealso \code{\link{estimate_weights}}, \code{\link{boot}}
#' @export
bootstrap_OR <- function(intervention_data, comparator_data, matching, i, model){

  # Samples the data
  bootstrap_data <- intervention_data[i,]

  # Estimates weights
  perform_wt <- estimate_weights(intervention_data=bootstrap_data, matching_vars=matching,  comparator_data=comparator_data)

  # Perform logistic regression and extract the OR estimate
  logistic.regr <- suppressWarnings(glm(formula = model, family=binomial(link="logit"), data = perform_wt$analysis_data, weight = wt))
  OR <- exp(as.numeric(coef(logistic.regr)[2]))
}

