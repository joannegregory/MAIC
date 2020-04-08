# Functions for the estimation of propensity weights


# Internal functions - Not exported ---------------------------------------
# Objective function
objfn <- function(a1, X){
  sum(exp(X %*% a1))
}

# Gradient function
gradfn <- function(a1, X){
  colSums(sweep(X, 1, exp(X %*% a1), "*"))
}

# External functions ------------------------------------------------------

#' Estimate MAIC propensity weights
#'
#' Estimate propensity weights for matched-adjusted indirect comparison (MAIC)
#'
#' @param intervention_data A data frame containing individual patient data from the intervention study.
#' @param baseline_comparator  A data frame of average covariate values from the comparator study.
#' @param matching_vars A character vector giving the names of the covariates to use
#'   in matching. These names must match the column names in intervention_data
#'   and baseline_comparator.
#' @param ... Additional arguments to be passed to optimisation functions such
#'   as the method for maximum likelihood optimisation. The default is method =
#'   "BFGS". Refer to \code{\link[stats]{optim}} for options.
#'
#' @return A data frame named intervention_wts including the intervention data
#'   with additional columns named WT (weights) and WT_RS (rescaled weights).
#'
#' @references NICE DSU TECHNICAL SUPPORT DOCUMENT 18: METHODS FOR
#'   POPULATION-ADJUSTED INDIRECT COMPARISONS IN SUBMSISSIONS TO NICE, REPORT BY
#'   THE DECISION SUPPORT UNIT, December 2016
#'
#' @export
estimate_weights <- function(intervention_data, baseline_comparator , matching_vars, ...){
  #Centre covariates

  #Estimate weights

  #Rescale weights

  #Return a data frame with weights attached intervention_wts
}

# Functions for summarizing the weights ---------------------------------

#' Estimate effective sample size
#'
#' Estimate the effective sample size (ESS).
#'
#' @param intervention_wts A data frame containing individual patient data from
#'   the intervention study, including a column containing the weights (derived
#'   using estimate_weights).
#' @param wt_col The name of the weights column in the data frame containing the
#'   intervention individual patient data and the MAIC propensity weights. The
#'   default is WT.
#'
#' @return The effective sample size (ESS) as a numeric value.
#'
#' @references NICE DSU TECHNICAL SUPPORT DOCUMENT 18: METHODS FOR
#'   POPULATION-ADJUSTED INDIRECT COMPARISONS IN SUBMSISSIONS TO NICE, REPORT BY
#'   THE DECISION SUPPORT UNIT, December 2016
#'
#' @seealso \code{\link{estimate_weights}}
#' @export
estimate_ess <- function(intervention_wts, wt_col=WT){

  # Estimate ESS [ess <- (sum(intervention_wts$wt_col)^2/sum(intervention_wts$wt_col^2))]

  # Return ESS
}

#' Summarize the weight values
#'
#' Produce a summary of the weights (minimum, maximum, median, mean, sd). Mean
#' and standard deviation are provided for completeness. In practice the
#' distribution of weights may be skewed in which case mean and sd should be
#' interpreted with caution
#'
#' @param intervention_wts A data frame containing individual patient data from
#'   the intervention study, including a column containing the weights (derived
#'   using estimate_weights).
#' @param wt_col The name of the weights column in the data frame containing the
#'   intervention individual patient data and the MAIC propensity weights. The
#'   default is WT.
#'
#' @return A data frame that inlcudes a summary (minimum, maximum, median, mean) of the weights.
#'
#' @seealso \code{\link{estimate_weights}}
#' @export
summarize_wts <- function(intervention_wts, wt_col=WT){

  # Summarize the minimum, maximum, median, mean and SD of the weights

  # Return data frame summarizing the minimum, maximum, median, mean and SD of the weights
}

#' Produce histograms of weights and rescaled weights
#'
#' Produce a plot containing two histograms (one of the weights and one of the recaled weights).
#'
#' @param intervention_wts A data frame containing individual patient data from
#'   the intervention study, including a column containing the weights (derived
#'   using estimate_weights).
#' @param wt_col The name of the weights column in the data frame containing the
#'   intervention individual patient data and the MAIC propensity weights. The
#'   default is WT.
#' @param rs_wt_col The name of the rescaled weights column in the data frame
#'   containing the intervention individual patient data and the MAIC propensity
#'   weights. The default is WT_RS.
#'
#' @return A data frame that inlcudes a summary (minimum, maximum, median, mean, sd) of the weights.
#'
#' @seealso \code{\link{estimate_weights}}
#' @export
hist_wts <- function(intervention_wts, wt_col=WT, rs_wt_col=WT_RS){

  # Plot histogram of weights

  # Plot histogram of rescaled weights

  # Combine histograms onto one plot (one on top of the other)

  # Return ggplot object containing a histogram of the weights and a histogram
  # of the recaled weights
}

#' Produce a data frame of the weights assigned to alternative patient profiles
#'
#' Select the patient characteristics used in the matching and the MAIC weights
#' and output a data frame of unique propensity weight values with the
#' associated summary baseline characteristics. This data frame helps to
#' understand how different patient profiles are contributing to the analyses by
#' illustrating the patient characteristics associated with different weight
#' values. For example, min, max and median weights. This function is most
#' useful when only matching on binary variables as there are fewer unique
#' values.
#'
#' @param intervention_wts A data frame containing individual patient data from
#'   the intervention study, including a column containing the weights (derived
#'   using estimate_weights).
#' @param wt_col The name of the weights column in the data frame containing the
#'   intervention individual patient data and the MAIC propensity weights. The
#'   default is WT.
#' @param matching_vars A character vector giving the names of the covariates to
#'   use in matching. These names must match the column names in
#'   intervention_data and baseline_comparator.
#'
#' @return A data frame that includes a summary of patient characteristics
#'   associated with alternative weight values.
#'
#' @seealso \code{\link{estimate_weights}}
#' @export
profile_wts <- function(intervention_wts, wt_col=WT, matching_vars){

  # Filter the intervention_wts data frame to only contain each unique set of
  # patient characteritics with the weight assigned

  # Sort the data frame by weight (smallest to largest)

  # Return data frame containing a column with a unique set of weights and the
  # set of patient characteristics that each weight was assigned to
}

#' Weight diagnostics
#'
#' Produce a set of useful diagnostic metrics to summarize propensity weights
#' \itemize{
#'   \item ESS
#'   \item Summary statistics of the weights: minimum, maximum, median, mean, sd
#'   \item Patient profile associated with given weight values
#' }
#'
#' This function calls each of estimate_ess, summarize_wts and profile_wts in turn.
#'
#' @param intervention_wts A data frame containing individual patient data from
#'   the intervention study, including a column containing the weights (derived
#'   using estimate_weights).
#' @param wt_col The name of the weights column in the data frame containing the
#'   intervention individual patient data and the MAIC propensity weights. The
#'   default is WT.
#' @param rs_wt_col The name of the rescaled weights column in the data frame
#'   containing the intervention individual patient data and the MAIC propensity
#'   weights. The default is WT_RS.
#' @param matching_vars A character vector giving the names of the covariates to
#'   use in matching.
#'
#' @return A list of data frames to use for weight diagnostics.
#'
#' @seealso \code{\link{estimate_weights}}, \code{\link{estimate_ess}}, \code{\link{summarize_wts}}, \code{\link{profile_wts}}
#' @export
wt_diagnostics <- function(intervention_wts, wt_col=WT, rs_wt_col=WT_RS, matching_vars){

  # Calls functions: estimate_ess, summarize_wts and profile_wts

  # Run each function

  # Return a list of data frames containing the outputs produced by estimate_ess, summarize_wts and profile_wt
}





