#' Derive analysis dataset for weighted analyses
#'
#' Creates analysis dataset including weighted intervention data and comparator
#' pseudo data with propensity weights
#'
#' @param intervention_wts A data frame containing individual patient data
#'   from the intervention study, including a column containing the weights
#'   (derived using estimate_weights).
#' @param comparator_data A data frame containing the pseudo-individual patient
#'   data for the comparator intervention.
#' @param wt_col The name of the weights column in the data frame containing the
#'   intervention individual patient data and the MAIC propensity weights. The
#'   default is WT.
#' @param rs_wt_col The name of the rescaled weights column in the data frame
#'   containing the intervention individual patient data and the MAIC propensity
#'   weights. The default is WT_RS.
#'
#' @return A data frame including the intervention data with propensity weights
#'   and comparator psuedo data with weights=1.
#'
#' @seealso \code{\link{estimate_weights}}
#' @export
analysis_dataset <- function(intervention_wts, comparator_data, wt_col=WT, rs_wt_col=WT_RS){

  # Combine intervention_wts and comparator_data and set weights in comparator_data=1

  # Return a data frame named analysis_data containing intervention_wts and comparator_data

}

