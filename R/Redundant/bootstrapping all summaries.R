
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

boostrap_HR <- function(intervention_data, i,  matching_vars, comparator_data, binary_var){

  # Samples the data
  bootstrap_data <- intervention_data[i,]

  # Estimates weights
  perform_wt <- estimate_weights(intervention_data=bootstrap_data, matching_vars,  comparator_data)

  # survival data stats
  fit.A <- surv_fit(Surv(Time, Event) ~ 1, data = perform_wt$intervention_wt_data, weights=perform_wt$intervention_wt_data$wt)
  median_A<-surv_median(fit.A)$median
  #
  fit.B <- surv_fit(Surv(Time, Event) ~ 1, data = perform_wt$comparator_wt_data)
  median_B<-surv_median(fit.B)$median

  cox_model <- coxph(Surv(Time, Event==1) ~ ARM, data = perform_wt$analysis_data, weights = wt)
  HR <- exp(cox_model$coefficients)

  c("HR"=HR, "median_A"=median_A, "median_B"=median_B)
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
#' @param analysis_data  A data frame containing data with weights (derived from
#'   \code{\link{analysis_dataset}}).
#' @param n_sim  Number of simulations. Default = 1000.
#'
#' @return A data frame of n_sim bootstrapped OR values.
#'
#' @seealso \code{\link{analysis_dataset}}, \code{\link{estimate_weights}}, \code{\link{OR_summary}}
#' @export
boostrap_OR <- function(intervention_data, i,  matching_vars, comparator_data, binary_var){


  # Samples the data
  bootstrap_data <- intervention_data[i,]

  # Estimates weights
  perform_wt <- estimate_weights(intervention_data=bootstrap_data, matching_vars=matching_vars,  comparator_data=comparator_data)


  # binary data stats - note we can probably get rid of some of these, just to show that the OR works using a logitsic regression or manually shall we also add RR?
  prop_A<-weighted.mean(perform_wt$intervention_wt_data[,binary_var],perform_wt$intervention_wt_data$wt)
  prop_B<-mean(perform_wt$comparator_wt_data[,binary_var])
  logistic.regr_RR <- suppressWarnings(glm(Binary_event~ARM, family=poisson(link="log"), data = perform_wt$analysis_data, weight = wt))
  RR <- exp(as.numeric(coef(logistic.regr_RR)[2]))
  RR_test <- prop_A/prop_B

  odds_A <- weighted.mean(perform_wt$intervention_wt_data[,binary_var],perform_wt$intervention_wt_data$wt)/weighted.mean(1-perform_wt$intervention_wt_data[,binary_var],perform_wt$intervention_wt_data$wt)
  odds_B<-mean(perform_wt$comparator_wt_data[,binary_var])/mean(1-perform_wt$comparator_wt_data[,binary_var])
  logistic.regr <- suppressWarnings(glm(Binary_event~ARM, family=binomial(link="logit"), data = perform_wt$analysis_data, weight = wt))
  OR <- exp(as.numeric(coef(logistic.regr)[2]))
  OR_test <- odds_A/odds_B


  c("OR" = OR,  "RR" = RR, "odds_A"=odds_A, "odds_B"=odds_B, "prop_A" = prop_A, "prop_B" = prop_B, "OR_test"=OR_test, "RR_test" = RR_test)
}
