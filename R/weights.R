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
#' Estimate propensity weights for matching-adjusted indirect comparison (MAIC).
#'
#' @param intervention_data A data frame containing individual patient data from the intervention study.
#' @param comparator_data  A data frame containing pseudo individual patient data from the comparator study.
#'  The outcome variables names must match intervention_data.
#' @param matching_vars A character vector giving the names of the covariates to use
#'   in matching. These names must match the column names in intervention_data.
#' @param ... Additional arguments to be passed to optimisation functions such
#'   as the method for maximum likelihood optimisation. The default is method =
#'   "BFGS". Refer to \code{\link[stats]{optim}} for options.
#'
#' @return A list containing 2 objects. First, a data frame named analysis_data containing intervention_data and comparator_data
#'   with additional columns named wt (weights) and wt_rs (rescaled weights) where the weights are 1 for the comparator data.
#'   Second, a vector called matching_vars of the matching variables names used.
#' @references NICE DSU TECHNICAL SUPPORT DOCUMENT 18: METHODS FOR
#'   POPULATION-ADJUSTED INDIRECT COMPARISONS IN SUBMSISSIONS TO NICE, REPORT BY
#'   THE DECISION SUPPORT UNIT, December 2016
#' @seealso \code{\link{optim}}
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
  #       - A character  vector with the name of the matching variables
  output <- list(
    matching_vars = matching_vars,
    analysis_data = all_data
  )

  return(output)

}

# Functions for summarizing the weights ---------------------------------

#' Estimate effective sample size
#'
#' Estimate the effective sample size (ESS).
#'
#' @param data A data frame containing individual patient data from
#'   the intervention study, including a column containing the weights (derived
#'   using \code{\link{estimate_weights}}).
#' @param wt_col The name of the weights column in the data frame containing the
#'   intervention individual patient data and the MAIC propensity weights. The
#'   default is wt.
#'
#' @return The effective sample size (ESS) as a numeric value.
#'
#' @references NICE DSU TECHNICAL SUPPORT DOCUMENT 18: METHODS FOR
#'   POPULATION-ADJUSTED INDIRECT COMPARISONS IN SUBMSISSIONS TO NICE, REPORT BY
#'   THE DECISION SUPPORT UNIT, December 2016
#'
#' @seealso \code{\link{estimate_weights}}
#' @export
estimate_ess <- function(data, wt_col="wt"){
  ess <- sum(data[,wt_col])^2/sum(data[,wt_col]^2)
  return(ess)
}


#' Summarize the weight values
#'
#' Produce a summary of the weights (minimum, maximum, median, mean, SD). Mean
#' and standard deviation are provided for completeness. In practice the
#' distribution of weights may be skewed in which case mean and SD should be
#' interpreted with caution.
#'
#' @param data A data frame containing individual patient data from
#'   the intervention study, including a column containing the weights (derived
#'   using \code{\link{estimate_weights}}).
#' @param wt_col The name of the weights column in the data frame containing the
#'   intervention individual patient data and the MAIC propensity weights. The
#'   default is wt.
#' @param rs_wt_col The name of the rescaled weights column in the data frame
#'   containing the intervention individual patient data and the MAIC propensity
#'   weights. The default is wt_rs.
#'
#' @return A data frame that includes a summary (minimum, maximum, median, mean) of the weights and rescaled weights.
#'
#' @seealso \code{\link{estimate_weights}}
#' @export
summarize_wts <- function(data, wt_col="wt", rs_wt_col="wt_rs"){
  summary <- data.frame(
    type = c("Weights", "Rescaled weights"),
    mean = c(mean(data[,wt_col]), mean(data[,rs_wt_col])),
    sd = c(sd(data[,wt_col]), sd(data[,rs_wt_col])),
    median = c(median(data[,wt_col]), median(data[,rs_wt_col])),
    min = c(min(data[,wt_col]), min(data[,rs_wt_col])),
    max = c(max(data[,wt_col]), max(data[,rs_wt_col]))
  )
  return(summary)
}


#' Produce histograms of weights and rescaled weights
#'
#' Produce a plot containing two histograms (one of the weights and one of the rescaled weights).
#'
#' @param data A data frame containing individual patient data from
#'   the intervention study, including a column containing the weights (derived
#'   using \code{\link{estimate_weights}}).
#' @param wt_col The name of the weights column in the data frame containing the
#'   intervention individual patient data and the MAIC propensity weights. The
#'   default is wt.
#' @param rs_wt_col The name of the rescaled weights column in the data frame
#'   containing the intervention individual patient data and the MAIC propensity
#'   weights. The default is wt_rs.
#' @param bin Number of bins to plot histogram. The default is 30.
#'
#' @return A histogram plot of the weights and rescaled weights.
#'
#' @seealso \code{\link{estimate_weights}}
#' @export
hist_wts <- function(data, wt_col="wt", rs_wt_col="wt_rs", bin = 30) {

  wt_data <- data %>%
    dplyr::select(c(wt_col, rs_wt_col)) %>% # select only the weights and rescaled weights
    dplyr::rename("Weights" = wt_col, "Rescaled weights" = rs_wt_col) %>% # rename so for plots
    tidyr::gather() # weights and rescaled weights in one column for plotting


  hist_plot <- ggplot2::ggplot(wt_data) + ggplot2::geom_histogram(ggplot2::aes(value), bins = bin) +
    ggplot2::facet_wrap(~key,  ncol=1) + # gives the two plots (one on top of the other)
    ggplot2::theme_bw()+
    ggplot2::theme(axis.title = ggplot2::element_text(size = 16),
                   axis.text = ggplot2::element_text(size = 16)) +
    ggplot2::ylab("Frequency") +
    ggplot2::xlab("Weight")

  return(hist_plot)
}


#' Produce a data frame of the weights assigned to patient profiles
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
#' @param data A data frame containing individual patient data from
#'   the intervention study, including a column containing the weights (derived
#'   using \code{\link{estimate_weights}}).
#' @param wt_col The name of the weights column in the data frame containing the
#'   intervention individual patient data and the MAIC propensity weights. The
#'   default is wt.
#' @param rs_wt_col The name of the rescaled weights column in the data frame
#'   containing the intervention individual patient data and the MAIC propensity
#'   weights. The default is wt_rs.
#' @param vars A character vector giving the variable names of the baseline
#'   characteristics. These names must match the column names in data.
#'
#' @return A data frame that includes a summary of patient characteristics
#'   associated with each weight value.
#'
#' @seealso \code{\link{estimate_weights}}
#' @export
profile_wts <- function(data, wt_col="wt", wt_rs="wt_rs", vars){
  profile_data <-  data %>%
    dplyr::select(vars, wt_col, wt_rs)

  profile_wts <- profile_data %>%
    dplyr::distinct()

  return(profile_wts)
}

#' Weight diagnostics
#'
#' Produce a set of useful diagnostic metrics to summarize propensity weights
#' \itemize{
#'   \item ESS (\code{\link{estimate_ess}})
#'   \item Summary statistics of the weights: minimum, maximum, median, mean, SD (\code{\link{summarize_wts}})
#'   \item Patient profile associated with weight values (\code{\link{profile_wts}})
#' }
#'
#' @param data A data frame containing individual patient data from
#'   the intervention study, including a column containing the weights (derived
#'   using estimate_weights).
#' @param wt_col The name of the weights column in the data frame containing the
#'   intervention individual patient data and the MAIC propensity weights. The
#'   default is wt.
#' @param rs_wt_col The name of the rescaled weights column in the data frame
#'   containing the intervention individual patient data and the MAIC propensity
#'   weights. The default is wt_rs.
#' @param vars A character vector giving the variable names of the baseline
#'   characteristics. These names must match the column names in data.
#'
#' @return List of the following:
#' \itemize{
#'   \item The effective sample size (ESS) as a numeric value.
#'   \item A data frame that includes a summary (minimum, maximum, median, mean) of the weights and rescaled weights.
#'   \item A data frame that includes a summary of patient characteristics
#'   associated with each weight value.
#' }
#'
#' @seealso \code{\link{estimate_weights}}, \code{\link{estimate_ess}}, \code{\link{summarize_wts}}, \code{\link{profile_wts}}
#' @export
wt_diagnostics <- function(data, wt_col="wt", wt_rs="wt_rs", vars){

  # ESS
  ESS <- estimate_ess(data, wt_col)

  # Summary
  summ_wts <- summarize_wts(data, wt_col)

  # Weight profiles
  profile <- profile_wts(data, wt_col, wt_rs, vars)

  output <- list("ESS" = ESS,
                 "Summary_of_weights" = summ_wts,
                 "Weight_profiles" = profile
  )
  return(output)
}






