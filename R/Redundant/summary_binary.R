# Function to summarize results in terms of ORs

#' Summary table of odds ratio estimates
#'
#' Returns table of odds ratios with 95\% CI: no weighting (naive logistic
#' regression), weighted logistic regression model using fixed weights and a
#' weighted logistic regression using bootstrapped weights
#'
#' @param analysis_data  A data frame containing data with weights (derived from
#'   \code{\link{analysis_dataset}}).
#' @param bootstrap Logical. If TRUE (Default) perform bootstrapping of the
#'   propensity weights.
#' @param bootstrap_method A character string specifying the bootstrap method
#'   to be used. Default = 'percentile'. This is ignored if bootstrap = FALSE.
#' @param ... Additional arguments to be passed to the bootstrapping functions.
#'   For example the number of simulations.
#'
#' @return A data frame of ORs with 95\% CI.
#'
#' @seealso \code{\link{analysis_dataset}}, \code{\link{bootstrap_OR}}
#' @export
OR_summary <- function(analysis_data, bootstrap = TRUE, bootstrap_method = 'percentile', ...){

  # Calculate naive OR and 95% CI from a logistic regression model

  # Calculate weighted OR and 95% CI from a weighted logistic regression model

  if(bootstrap) {
    #perform bootstrapping to get a sample of OR values
    bootstrapped_OR <- bootstrap_OR(analysis_dataset, n_sim = 1000)

    if(bootstrap_method == 'percentile'){
      # Calculate median and 95% CI based on percentiles of bootstrapped OR values
    }

    if(bootstrap_method == 'TBC'){ #Alternative method to be confirmed
      # Calculate alternative bootstrapped OR
    }

  }

  # Combine the above results into a data frame

  # Return a data frame of ORs
}


#' Bootstrapping for MAIC propensity weighted odds ratios
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
#' @return A data frame n_sim bootstrapped OR values
#'
#' @seealso \code{\link{analysis_dataset}}, \code{\link{estimate_weights}}, \code{\link{OR_summary}}
#' @export
bootstrap_OR <- function(analysis_data, n_sim = 1000){

  # Sample patients with replacement from the intervention arm of analysis_data

  # Perform the weighting using estimate_weights

  # Calculate the OR using a logistic regression model

  # repeat above steps above n_sim times (in a loop)

  # Return data frame of n_sim bootstrap ORs
}
