
#' Bootstrapping for MAIC propensity weighted hazard ratios
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
#' Refer to \code{\link{HR_summary}} for alternative approaches to utilising
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
bootstrap_HR <- function(analysis_dataset, n_sim = 1000){

  # Samples patients with replacement from intervention arm of analysis_dataset

  # perform matching using estimate_weights

  # calculate the HR

  # repeat above steps above n_sim times (in a loop)

  # Return data frame of bootstrapped HRs
}
