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
estimate_weights <- function(intervention_data, matching_vars, ...){

  opt1 <- optim(par = rep(0,dim(intervention_data[, matching_vars])[2]),
                      fn = objfn,
                      gr = gradfn,
                      X = as.matrix(intervention_data[, matching_vars]),
                      method = "BFGS")

  a1 <- opt1$par


  # Calculation of weights.
  WT <- as.vector(exp(as.matrix(intervention_data[, matching_vars]) %*% a1))

  # rescaled weights
  WT_RS <- (WT / sum(WT)) * dim(intervention_data)[1]

  # combine intervention_data with weights
  intervention_wts <- cbind(intervention_data, WT, WT_RS)

  return(intervention_wts)
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

  ess <- (sum(intervention_wts$wt_col)^2/sum(intervention_wts$wt_col^2))

  return(ess)
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
    summary <- intervention_wts %>%
                summarise(
                          min = min(intervention_wts$wt),
                          max = max(intervention_wts$wt),
                          median = median(intervention_wts$wt),
                          mean = mean(intervention_wts$wt),
                          sd = sd(intervention_wts$wt)
                          )
    return(summary)
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
#' @return A histogram plot of the weights and rescaled weights.
#'
#' @seealso \code{\link{estimate_weights}}
#' @export
hist_wts <- function(intervention_wts, wt_col="WT", rs_wt_col="WT_RS"){

  wt_data <- data[,c(wt_col, rs_wt_col)] %>% # select the weights and rescaled weights only
    rename("Weights" = wt_col, "Rescaled weights" = rs_wt_col) %>% # rename the weight columns for histogram output
    gather() # change the data so there is one column of weights and another to define whether the weight is weight/rescaled weight

  hist_plot <- qplot(data = wt_data,
                     value,
                     geom="histogram",
                     xlab = "Histograms of weights and rescaled weights",
                     binwidth=0.05) +
    facet_wrap(~key,  ncol=1) + # creates two hisotgrams (one for weights and one for rescaled weights)
    theme_bw()+
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 16))+
    ylab("Frequency")

  return(hist_plot)
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
#' @param rs_wt_col The name of the rescaled weights column in the data frame
#'   containing the intervention individual patient data and the MAIC propensity
#'   weights. The default is WT_RS.
#' @param matching_vars A character vector giving the names of the covariates to
#'   use in matching. These names must match the column names in
#'   intervention_data and baseline_comparator.
#'
#' @return A data frame that includes a summary of patient characteristics
#'   associated with alternative weight values.
#'
#' @seealso \code{\link{estimate_weights}}
#' @export
profile_wts <- function(intervention_wts, wt_col="WT", rs_wt_col="WT_RS", matching_vars){

  profile_data <- intervention_wts[, c(wt_col, rs_wt_col, matching_vars)] %>%
                  distinct()
  return(profile_data)
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





