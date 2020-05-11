
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

bootstrap_HR <- function(intervention_data, comparator_data, matching, i){

  # Samples the data
  bootstrap_data <- intervention_data[i,]

  # Estimates weights
  perform_wt <- estimate_weights(intervention_data=bootstrap_data, matching_vars=matching, comparator_data=comparator_data)

  # survival data stat
  cox_model <- coxph(Surv(Time, Event==1) ~ ARM, data = perform_wt$analysis_data, weights = wt)
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
#' @param endpoint The variable name for the binary endpoint.
#'   Variable names need to match the corresponding columns in intervention_data
#'
#' @return A data frame of n_sim bootstrapped OR values.
#'
#' @seealso \code{\link{analysis_dataset}}, \code{\link{estimate_weights}}, \code{\link{OR_summary}}
#' @export
bootstrap_OR <- function(intervention_data, comparator_data, matching, i, endpoint){

  # Samples the data
  bootstrap_data <- intervention_data[i,]

  # Estimates weights
  perform_wt <- estimate_weights(intervention_data=bootstrap_data, matching_vars=matching,  comparator_data=comparator_data)

  # Perform logistic regression and extract the OR estimate
  logistic.regr <- suppressWarnings(glm(formula = endpoint~ARM, family=binomial(link="logit"), data = perform_wt$analysis_data, weight = wt))
  OR <- exp(as.numeric(coef(logistic.regr)[2]))
}



#' Bootstrapping for MAIC propensity weighted relative risk HELP PAGE NEEDS UPDATING
#'
#' Bootstrapping of the relative risk to capture the uncertainty in the
#' estimation of MAIC propensity weights using the following general process:
#' \enumerate{
#'    \item Sample patients with replacement from the intervention study
#'    \item Estimate propensity weights for the current sample of patients
#'    \item Calculate a weighted relative risk
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
#' @param endpoint The variable name for the binary endpoint.
#'   Variable names need to match the corresponding columns in intervention_data
#'
#' @return A data frame of n_sim bootstrapped OR values.
#'
#' @seealso \code{\link{analysis_dataset}}, \code{\link{estimate_weights}}, \code{\link{OR_summary}}
#' @export
bootstrap_RR <- function(intervention_data, comparator_data, matching, i, endpoint){

  # Samples the data
  bootstrap_data <- intervention_data[i,]

  # Estimates weights
  perform_wt <- estimate_weights(intervention_data=bootstrap_data, matching_vars=matching,  comparator_data=comparator_data)

  # Perform logistic regression and extract the RR estimate
  poisson_regr <- suppressWarnings(glm(formula = endpoint~ARM, family=poisson(link="log"), data = perform_wt$analysis_data, weight = wt))
  RR <- exp(as.numeric(coef(poisson_regr)[2]))
}

