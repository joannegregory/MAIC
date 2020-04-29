## These functions perform survival analyses on the MAIC dataset including KM summary and naive hazard ratios


#' Summarize Kaplan-Meier estimate of survival
#'
#' Return a summary of the Kaplan-Meier estimate of survival based on
#' \code{\link[survival]{survfit}} including number of patients, number of
#' events and median time to event
#'
#' @param survfit_object  A survival model object derived from
#'   \code{\link{survfit}}. This can be either a weighted or unweighted estimate.
#'
#' @return A summary table of the Kaplan-Meier estimate including number of patients, number of
#' events and median time to event
#'
#' @export
KM_summary <- function(survfit_object){

  # KM summary of median, (number of patients, events and median survival)

  # Returns a data frame
}

#' Summary table of hazard ratio estimates
#'
#' Returns table of hazard ratios with 95\% CI: no weighting (naive Cox
#' proportional hazards model), weighted Cox model using fixed weights and a
#' weighted Cox model using bootstrapped weights
#'
#' @param analysis_data A data frame containing data with weights (derived from
#'   \code{\link{analysis_dataset}}).
#' @param bootstrap Logical. If TRUE (Default) perform bootstrapping of the
#'   propensity weights.
#' @param bootstrap_method. A character string specifying the bootstrap method
#'   to be used. Default = 'percentile'. This is ignored if bootstrap = FALSE.
#' @param ... Additional arguments to be passed to the bootstrapping functions.
#'   For example the number of simulations.
#'
#' @return A summary table of hazard ratio estimates with 95\% CI
#'
#' @seealso \code{\link{analysis_dataset}}, \code{\link{bootstrap_HR}}
#' @export
HR_summary <- function(analysis_data, bootstrap = TRUE, bootstrap_method = 'percentile', ...){

  # Calculate naive HR from cox model

  # Calculate weighted HR from cox model

  if(bootstrap) {
    #perform bootstrapping to get a sample of HR values
    bootstrapped_HR <- bootstrap_HR(MAIC_analysis_dataset, n_sim = 1000)

    if(bootstrap_method == 'percentile'){
      # Calculate median and 95% CI based on percentiles of bootstrapped HR values
    }

    if(bootstrap_method == 'TBC'){ #Alternative method to be confirmed
      # Calculate alternative bootstrapped HR
    }

  }
  # Combine the above

  # Return table of hazard ratios
}


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
