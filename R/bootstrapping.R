
#' Bootstrapping for MAIC propensity weighted hazard ratios HELP PAGE NEEDS UPDATING
#'
#' Bootstrapping of the hazard ratio to capture the uncertainty in the
#' estimation of MAIC propensity weights using the following general process:
#' \enumerate{
#'    \item Sample patients with replacement from the intervention study
#'    \item Estimate propensity weights for the current sample of patients
#'    \item Calculate a weighted hazard ratio
#'    \item Repeat the steps above N times. Approx. 1000 simulations are usually
#'    sufficient
#' }
#' bootstrapped samples
#'
#' @param analysis_data  A data frame containing data with weights (derived from
#'   \code{\link{analysis_dataset}}).
#' @param n_sim  Number of simulations. Default = 1000.
#'
#' @return A data frame of n_sim bootstrapped HR values.
#'
#' @seealso \code{\link{analysis_dataset}}, \code{\link{estimate_weights}}, \code{\link{HR_summary}}
#' @export

boostrap_HR <- function(intervention_data, comparator_data, matching, i, model){

  # Samples the data
  bootstrap_data <- intervention_data[i,]

  # Estimates weights
  perform_wt <- estimate_weights(intervention_data=bootstrap_data, matching_vars=matching, comparator_data=comparator_data)

  # survival data stat
  cox_model <- coxph(model, data = perform_wt$analysis_data, weights = wt)
  HR <- exp(cox_model$coefficients)
}




#' Bootstrapping for MAIC propensity weighted odds ratios HELP PAGE NEEDS UPDATING
#'
#' Bootstrapping of the odds ratio to capture the uncertainty in the
#' estimation of MAIC propensity weights using the following general process:
#' \enumerate{
#'    \item Sample patients with replacement from the intervention study
#'    \item Estimate propensity weights for the current sample of patients
#'    \item Calculate a weighted odds ratio
#'    \item Repeat the steps above N times. Approx. 1000 simulations are usually
#'    sufficient
#' }
#' Refer to \code{\link{OR_summary}} for alternative approaches to utilising
#' bootstrapped samples
#'
#' @param intervention_data
#' @param comparator_data
#' @param matching
#' @param i
#' @param model A model formula in the form 'endpoint ~ treatment_var'.
#'   Variable names need to match the corresponding columns in intervention_data
#'
#' @return A data frame of n_sim bootstrapped OR values.
#'
#' @seealso \code{\link{analysis_dataset}}, \code{\link{estimate_weights}}, \code{\link{OR_summary}}
#' @export
boostrap_OR <- function(intervention_data, comparator_data, matching, i, model){

  # Samples the data
  bootstrap_data <- intervention_data[i,]

  # Estimates weights
  perform_wt <- estimate_weights(intervention_data=bootstrap_data, matching_vars=matching,  comparator_data=comparator_data)

  # Perform logistic regression and extract the OR estimate
  logistic.regr <- suppressWarnings(glm(formula = model, family=binomial(link="logit"), data = perform_wt$analysis_data, weight = wt))
  OR <- exp(as.numeric(coef(logistic.regr)[2]))
}
