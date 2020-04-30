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

#' Estimate MAIC propensity weights HELP PAGE NEEDS UPDATING
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
#' @return A list containing 2 objects. First, a data frame named intervention_wts_data including the intervention data
#'   with additional columns named WT (weights) and WT_RS (rescaled weights). Second, a vector of the matching variables names used.
#'
#' @references NICE DSU TECHNICAL SUPPORT DOCUMENT 18: METHODS FOR
#'   POPULATION-ADJUSTED INDIRECT COMPARISONS IN SUBMSISSIONS TO NICE, REPORT BY
#'   THE DECISION SUPPORT UNIT, December 2016
#'
#' @export
estimate_weights <- function(intervention_data, comparator_data, matching_vars){

  #Basic checks of inputs before proceeding
  #Check intervention data is a data frame
  assertthat::assert_that(
    is.data.frame(intervention_data),
    msg = "intervention_data is expected to be a data frame"
  )
  #Check comparator data is a data frame
  assertthat::assert_that(
    is.data.frame(comparator_data),
    msg = "comparator_data is expected to be a data frame"
  )
  #Check that matching_vars is a character vector
  assertthat::assert_that(
    is.character(matching_vars),
    msg = "matching_vars is expected to be a character vector"
  )
  #Check that all named matching variables are in the intervention dataset
  assertthat::assert_that(
    all(matching_vars %in% colnames(intervention_data)),
    msg = "matching_vars contains variable names that are not in the intervention dataset"
  )

  assertthat::assert_that(
    all(colnames(comparator_data) %in% colnames(intervention_data)),
    msg = "Column names in the comparator dataset do not all match columns in the intervention dataset. Check your inputs"
  )

  # Optimise Q(b) using Newton-Raphson techniques
  opt1 <- optim(par = rep(0,dim(intervention_data[,matching_vars])[2]),
                fn = objfn,
                gr = gradfn,
                X = as.matrix(intervention_data[,matching_vars]),
                method = "BFGS")

  a1 <- opt1$par

  # Calculate weights for intervention data and combine with dataset
  data_with_wts <- dplyr::mutate(intervention_data,
                                 wt = as.vector(exp(as.matrix(intervention_data[,matching_vars]) %*% a1)), # weights
                                 wt_rs = (wt / sum(wt)) * nrow(intervention_data), # rescaled weights
                                 ARM = "Intervention"
  )

  # assign weight=1 to comparator data
  comparator_data_wts <- comparator_data %>%
    dplyr::mutate(wt=1, wt_rs=1, ARM="Comparator")

  # Join comparator data with the intervention data
  all_data <- rbind.fill(data_with_wts, comparator_data_wts)
  all_data$ARM <- relevel(as.factor(all_data$ARM), ref="Comparator")

  # Outputs are:
  #       - the analysis data (intervention PLD, weights and comparator pseudo PLD)
  #       - A charcacter vector with the name of the centered matching variables
  #       - A charcacter vector with the name of the matching variables
  output <- list(
    matching_vars = matching_vars,
    analysis_data = all_data
  )

  return(output)

}

# Functions for summarizing the weights ---------------------------------

#' Estimate effective sample size HELP PAGE NEEDS UPDATING
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
estimate_ess <- function(data, wt="wt"){
  ess <- sum(data[,wt])^2/sum(data[,wt]^2)
  return(ess)
}


#' Summarize the weight values HELP PAGE NEEDS UPDATING
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
summarize_wts <- function(data, wt="wt"){
  summary <- data.frame(
    mean = mean(data[,wt]),
    sd = sd(data[,wt]),
    median = median(data[,wt]),
    min = min(data[,wt]),
    max = max(data[,wt])
  )
  return(summary)
}


#' Produce histograms of weights and rescaled weights HELP PAGE NEEDS UPDATING
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
hist_wts <- function(data, wt_col="wt", rs_wt_col="wt_rs", bin = 30) {

  wt_data <- data %>%
    dplyr::select(c(wt_col, rs_wt_col)) %>% # select only the weights and rescaled weights
    rename("Weights" = wt_col, "Rescaled weights" = rs_wt_col) %>% # rename so for plots
    gather() # weights and rescaled weights in one column for plotting


  hist_plot <- ggplot2::ggplot(wt_data) + ggplot2::geom_histogram(aes(value), bins = bin) +
    ggplot2::facet_wrap(~key,  ncol=1) + # gives the two plots (one on top of the other)
    ggplot2::theme_bw()+
    ggplot2::theme(axis.title = element_text(size = 16),
                   axis.text = element_text(size = 16)) +
    ggplot2::ylab("Frequency")

  return(hist_plot)
}


#' Produce a data frame of the weights assigned to alternative patient profiles HELP PAGE NEEDS UPDATING
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
profile_wts <- function(data, wt="wt", wt_rs="wt_rs", vars){
  profile_data <-  data %>%
    select(vars, wt, wt_rs)

  profile_wts <- profile_data %>%
    distinct() %>%
    arrange(wt)

  return(profile_wts)
}

#' Weight diagnostics HELP PAGE NEEDS UPDATING
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
wt_diagnostics <- function(data, matching_vars, wt="wt"){

  # ESS
  ESS <- estimate_ess(data)

  # Summary
  summ_wts <- summarize_wts(data)

  # Weight profiles
  profile <- profile_wts(data, vars=matching_vars)

  output <- list("ESS" = ESS,
                 "Summary_of_weights" = summ_wts,
                 "Weight_profiles" = profile
  )
  return(output)
}






