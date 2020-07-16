#' Bootstrapping for MAIC weighted hazard ratios
#'
#' A function required for the "statistic" argument in the \code{\link{boot}} function.
#' Performs MAIC weighting using {\link{estimate_weights}} and returns a weighted hazard ratio (HR) from a Cox proportional hazards model.
#' @param intervention_data  A data frame containing individual patient data from the intervention study.
#' @param matching A character vector giving the names of the covariates to use
#'   in matching. These names must match the column names in intervention_data.
#' @param i Index used to select a sample within \code{\link{boot}}.
#' @param model A model formula in the form 'Surv(Time, Event==1) ~ ARM'.
#'   Variable names need to match the corresponding columns in intervention_data.
<<<<<<< HEAD
#' @param comparator_data A data frame containing pseudo individual patient data from the comparator study needed to derive the relative treatment effect.
=======
#'   #' @param comparator_data A data frame containing pseudo individual patient data from the comparator study needed to derive the relative treatment effect.
>>>>>>> b4d1091cb48d6167115dd4b0ea6f1ded7a736654
#'  The outcome variables names must match intervention_data.

#'
#' @return The HR as a numeric value.
#'
#' @seealso \code{\link{estimate_weights}}, \code{\link{boot}}
#'
#' @example inst/examples/MAICexample.R
#'
#' @export

bootstrap_HR <- function(intervention_data, matching, i, model, comparator_data){

  # Samples the data
  bootstrap_data <- intervention_data[i,]

  # Estimates weights
  perform_wt <- estimate_weights(intervention_data=bootstrap_data, matching_vars=matching)

  # Give comparator data weights of 1
  comparator_data_wts <- comparator_data %>% dplyr::mutate(wt=1, wt_rs=1, ARM="Comparator")

  # Add the comparator data
  combined_data <- dplyr::bind_rows(perform_wt$analysis_data, comparator_data_wts)
  combined_data$ARM <- relevel(as.factor(combined_data$ARM), ref="Comparator")

  # survival data stat
  cox_model <- survival::coxph(model, data = combined_data, weights = wt)
  HR <- exp(cox_model$coefficients)
}




#' Bootstrapping for MAIC weighted odds ratios
#'
#' A function required for the "statistic" argument in the \code{\link{boot}} function.
#' Performs MAIC weighting using {\link{estimate_weights}} and returns a weighted odds ratio (OR) from a binomial generalised linear model.
#'
#' @param intervention_data  A data frame containing individual patient data from the intervention study.

#' @param matching A character vector giving the names of the covariates to use
#'   in matching. These names must match the column names in intervention_data.
#' @param i Index used to select a sample within \code{\link{boot}}.
#' @param model A model formula in the form 'endpoint ~ treatment_var'.
#'   Variable names need to match the corresponding columns in intervention_data.
#' @param comparator_data A data frame containing pseudo individual patient data from the comparator study needed to derive the relative treatment effect.
#'  The outcome variables names must match intervention_data.
<<<<<<< HEAD
#' @return The OR as a numeric value.
=======
#'  #' @return The OR as a numeric value.
>>>>>>> b4d1091cb48d6167115dd4b0ea6f1ded7a736654
#'
#' @seealso \code{\link{estimate_weights}}, \code{\link{boot}}
#'
#' @example inst/examples/MAICexample.R
#'
#' @export
bootstrap_OR <- function(intervention_data, matching, i, model, comparator_data){

  # Samples the data
  bootstrap_data <- intervention_data[i,]

  # Estimates weights
  perform_wt <- estimate_weights(intervention_data=bootstrap_data, matching_vars=matching)

  # Give comparator data weights of 1
  comparator_data_wts <- comparator_data %>% dplyr::mutate(wt=1, wt_rs=1, ARM="Comparator")

  # Add the comparator data
  combined_data <- dplyr::bind_rows(perform_wt$analysis_data, comparator_data_wts)
  combined_data$ARM <- relevel(as.factor(combined_data$ARM), ref="Comparator")

  # Perform logistic regression and extract the OR estimate
  logistic.regr <- suppressWarnings(glm(formula = model, family=binomial(link="logit"), data = combined_data, weight = wt))
  OR <- exp(as.numeric(coef(logistic.regr)[2]))
}
