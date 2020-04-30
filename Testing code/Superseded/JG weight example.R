estimate_weights <- function(intervention_data, matching_vars, compdata, ...){
  
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
  all_data <- rbind(intervention_wts, compdata)
  output <- list(alldata_wts_data=all_data, matching_vars=matching_vars)
  
  return(output)
}

weight_example <- estimate_weights(intervention_data=intervention, matching_vars=c("Smoke.centered", "ECOG0.centered"))

estimate.ess <- function(data, wt=WT){
  ess <- sum(data$WT)^2/sum(data$WT^2)
  return(ess)
}

data_for_diagnostics <- weight_example$alldata_wts_data %>% filter(trt=="A")
ESS <- estimate.ess(data=data_for_diagnostics)
ESS

all_diagnostics <- function(
  x <- weight_example$alldata_wts_data %>%
    filter(trt==intervention)
  outputlist(estimate.ess(x)
  weight_summ(x)...))
  return(output)
)
all_diagnostics(data=weight_example, intervention="A")

