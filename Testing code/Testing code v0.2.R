# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 2797 exploratory - testing out the code                         #
# Author: SS/ JG (23.04.2020)                                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Load libraries --------------------------------------------------
library(plyr)
library(dplyr)
library(haven)
library(ggplot2)
library(tidyr)
library(boot)
library(survminer)
library(survminer)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Directories -----------------------------------------------------
base_dir <- 'G:/Clients/Roche/2797 Development of R Code for MAIC and mixture cure models/Project'
data_path <- file.path(base_dir,'2 Exploratory/Simulated datasets')



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Load data -------------------------------------------------------

#### Intervention data

# Read in ADaM data and rename variables of interest
adsl <- read.csv(file.path(data_path, "adsl.csv")) %>% # subject level data
  mutate(SEX=ifelse(SEX=="Male",1,0))

adrs <- read.csv(file.path(data_path, "adrs.csv")) %>% # response data
  filter(PARAM=="Response") %>%
  select(USUBJID, ARM, Binary_event=AVAL) ## change "Binary_event" to "response"?

adtte <- read.csv(file.path(data_path, "adtte.csv")) %>% # time to event data
  filter(PARAMCD=="OS") %>%
  mutate(Event=1-CNSR) %>%
  select(USUBJID, ARM, Time=AVAL, Event)

# Combine all intervention data
intervention_input <- join_all(list(adsl, adrs, adtte), type = "full", by=c("USUBJID", "ARM"))

#### Comparator pseudo data

# read in pseudo survival data
comparator_surv <- read.csv(file.path(data_path,"psuedo_IPD.csv"))

# simulate response data based on the known proportion of responders
comparator_n <- 300 # total number of patients in the comparator data
comparator_prop_events <- 0.4 # proportion of responders
comparator_binary <- data.frame("Binary_event"=
                                  c(rep(1,comparator_n*comparator_prop_events),
                                    rep(0, comparator_n*(1-comparator_prop_events))))

# join survival and response comparator data note not a 1:1 relationship - the
# rows do not represent the same observation
comparator_input <- cbind(comparator_surv, comparator_binary)

# Baseline agregate data for the comparator popultion
target_pop <- read.csv(file.path(data_path,"Aggregate data.csv"))
target_pop



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Estimate weights, rescaled weight and join data -----------------

#### center baseline characteristics
# (subtract the aggregate comparator data from the corresponding column of intervention PLD)
intervention_data <- intervention_input %>%
                          mutate(Age_centered = AGE - target_pop$age.mean,
                                 # it is necessary to balance on both mean and standard deviation for continous variables:
                                 Age_squared_centered = (AGE^2) - (target_pop$age.mean^2 + target_pop$age.sd^2),
                                 Sex_centered = SEX - target_pop$prop.male,
                                 Smoke_centered = SMOKE - target_pop$prop.smoke,
                                 ECOG0_centered = ECOG0 - target_pop$prop.ecog0)
head(intervention_data)
names(intervention_data)

# cent_match_cov <- c("Age_centered",
#                     "Age_squared.centered",
#                     "Sex_centered",
#                     "Smoke_centered",
#                     "ECOG0_centered")

cent_match_cov <- names(intervention_data)[13:15]
cent_match_cov

match_cov <- names(intervention_data)[5:7]
match_cov


#### optimisation procedure and calculation of weights

# This is a fuction for Q(b)
objfn <- function(a1, X){
  sum(exp(X %*% a1))
}


# Gradient function => Derivative of Q(b).
gradfn <- function(a1, X){
  colSums(sweep(X, 1, exp(X %*% a1), "*"))
}


# Function to estimate weights, create analysis dataset and output a
# character vector of the variable names
estimate_weights <- function(intervention_data, cent_vars, vars, comparator_data){

  # Optimise Q(b) using Newton-Raphson techniques
  print(opt1 <- optim(par = rep(0,dim(intervention_data[,cent_vars])[2]),
                      fn = objfn,
                      gr = gradfn,
                      X = as.matrix(intervention_data[,cent_vars]),
                      method = "BFGS"))

  a1 <- opt1$par


  # Calculation of weights.
  wt <- as.vector(exp(as.matrix(intervention_data[,cent_vars]) %*% a1))

  # rescaled weights
  wt_rs <- (wt / sum(wt)) * dim(intervention_data)[1]

  # combine data with weights
  data_with_wts <- cbind(intervention_data, wt, wt_rs)

  # assign weight=1 to comparator data
  comparator_data_wts <- comparator_data %>% mutate(wt=1, wt_rs=1)

  # Join comparator data with the intervention data
  all_data <- rbind.fill(data_with_wts, comparator_data_wts)

  # Outputs are:
  #       - the analysis data (intervention PLD, weights and comparator pseudo PLD)
  #       - A charcacter vector with the name of the centered matching variables
  #       - A charcacter vector with the name of the matching variables
  output <- list(analysis_data = all_data,
                 centered_matching_vars = cent_vars,
                 matching_vars = vars
                 )

  return(output)

}


est_weights <- estimate_weights(intervention_data=intervention_data,
                                cent_vars = cent_match_cov,
                                vars = match_cov,
                                comparator_data=comparator_input)

head(est_weights$analysis_data)
est_weights$matching_vars
est_weights$centered_matching_vars

weighted.mean(intervention_wt_data$Binary_event, intervention_wt_data$wt_rs)
weighted.mean(intervention_wt_data$Binary_event, intervention_wt_data$wt)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Weight diagnostics ----------------------------------------------

# Weight diagnostics with the intervention data only:
intervention_wt_data <- est_weights$analysis_data %>%
                        filter(ARM=="A")

#### ESS

# Function to caclulate ESS:
estimate_ess <- function(data, wt=wt){
  ess <- sum(data$wt)^2/sum(data$wt^2)
  return(ess)
}

ESS <- estimate_ess(data=intervention_wt_data)
ESS


#### Weights summary

# Function to caclulate min, max, median, mean and sd of wts:
summarize_wts <- function(data, wt="wt"){
  summary <- data %>%
    summarise(
      min = min(data[,wt]),
      max = max(data[,wt]),
      median = median(data[,wt]),
      mean = mean(data[,wt]),
      sd = sd(data[,wt])
    )
  return(summary)
}

weight_summ <- summarize_wts(data=intervention_wt_data)
weight_summ

weight_rs_summ <- summarize_wts(data=intervention_wt_data, wt=wt_rs)
weight_rs_summ


##### Weight profiles

profile_wts <- function(data, wt=wt, wt_rs=wt_rs, vars){
  profile_data <-  data %>%
                  select(vars, wt, wt_rs)

  profile_wts <- profile_data %>%
                 distinct() %>%
                 arrange(wt)

  return(profile_wts)
}

wts_profile <- profile_wts(data=intervention_wt_data, vars = est_weights$matching_vars)


plot(wts_profile$AGE, wts_profile$wt)
boxplot(wts_profile$SEX , wts_profile$wt)

#### Histograms

# Function to plot a histogram of weights and rescaled weights:
hist_wts <- function(data, wt_col="wt", rs_wt_col="wt_rs", bin_width=NULL) {

wt_data <- data[,c(wt_col, rs_wt_col)] %>% # select only the weights and rescaled weights
            rename("Weights" = wt_col, "Rescaled weights" = rs_wt_col) %>% # rename so for plots
            gather() # weights and rescaled weights in one column for plotting

hist_plot <- qplot(data = wt_data,
                    value,
                    geom = "histogram",
                    xlab = "Histograms of weights and rescaled weights",
                    binwidth = bin_width
                   ) +
              facet_wrap(~key,  ncol=1) + # gives the two plots (one on top of the other)
              theme_bw()+
              theme(axis.title = element_text(size = 16),
                    axis.text = element_text(size = 16))+
              ylab("Frequency")

return(hist_plot)
}

histogram <- hist_wts(data=intervention_wt_data)
histogram

histogram_bw1 <- hist_wts(data=intervention_wt_data, bin_width=1)
histogram_bw1



#### All weight diagnostics

##### FUNCTION NOT WORKING
# Function to combine weight diagnostic functions above:
all_wt_diagnostics <- function(data, # analysis data from estimate_weights
                               #arm,
                               matching_vars,
                               wt=wt,
                               ...){

  # intervention_data <- data %>% filter(ARM==arm)

    # ESS
  ESS <- estimate_ess(data)

  # Summary
  summ_wts <- summarize_wts(data)

  # Weight profiles
  profile <- profile_wts(data, vars=matching_vars)

  # Histogram of weights
  hist_plot <- hist_wts(data, ...)

  output <- list("ESS" = ESS,
                 "Summary_of_weights" = summ_wts,
                 "Histogram_of_weights" = hist_plot,
                  "Weight_profiles" = profile
                 )
  return(output)
}


diagnostics <- all_wt_diagnostics(data=intervention_wt_data, matching_vars = est_weights$matching_vars, bin_width=0.1)

diagnostics$ESS
diagnostics$Summary_of_weights
diagnostics$Weight_profiles
diagnostics$Histogram_of_weights

# all_wt_diagnostics(data = est.weights$analysis_data,
#                    arm = "A",
#                    vars = est.weights$matching_vars)





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Bootstrapping ------------------------------------------------

boostrap_OR <- function(intervention_data, i,  cent_vars, comparator_data, binary_var){

  # Samples the data
  #i=c(1:nrow(intervention_data))
  bootstrap_data <- intervention_data[i,]

  # Estimates weights
  print(opt1 <- optim(par = rep(0,dim(bootstrap_data[, cent_vars])[2]),
                      fn = objfn,
                      gr = gradfn,
                      X = as.matrix(bootstrap_data[, cent_vars]),
                      method = "BFGS"))

  a1 <- opt1$par

  weights <- as.vector(exp(as.matrix(bootstrap_data[, cent_vars]) %*% a1))

  bootstrap_dat_wt <- cbind(bootstrap_data, weights)

  #Combines intervention dataset with comparator_data
  analysis_data <- bootstrap_dat_wt %>%
                   rbind.fill(comparator_data %>% mutate(ARM="Comparator",weights=1))
  #sort out reference

  # binary data
  prop_A<-weighted.mean(bootstrap_dat_wt[,binary_var],bootstrap_dat_wt$weights)
  prop_B<-mean(comparator_data[,binary_var])
  odds_A <- weighted.mean(bootstrap_dat_wt[,binary_var],bootstrap_dat_wt$weights)/weighted.mean(1-bootstrap_dat_wt[,binary_var],bootstrap_dat_wt$weights)
  odds_B<-mean(comparator_data[,binary_var])/mean(1-comparator_data[,binary_var])
  logistic.regr <- suppressWarnings(glm(Binary_event~ARM, family=binomial(link="logit"), data = analysis_data, weights = weights))
  OR <- exp(as.numeric(coef(logistic.regr)[2]))
  OR1 <- odds_A/odds_B

  log_OR <- as.numeric(coef(logistic.regr)[2])
  # sort out refereces

  c("OR" = OR, "Log_OR" = log_OR, "prop_A" = prop_A, "prop_B" = prop_B,"odds_A"=odds_A, "odds_B"=odds_B, "OR1"=OR1)
}


OR_bootstraps <- boot(intervention_data, boostrap_OR, R=10, cent_vars = cent_match_cov, comparator_data=comparator_input, binary_var="Binary_event")

boot_data <- as.data.frame(OR_bootstraps$t)
colnames(boot_data) <- colnames(t(as.data.frame(OR_bootstraps$t0)))

head(boot_data)
names(boot_data)
OR_bootstraps0
summary(OR_bootstraps)

hist(boot_data$OR, main = "",xlab = "Boostrapped OR")
abline(v= quantile(boot_data$OR,probs = c(0.025,0.5,0.975)), lty=2)

OR.LCI <- data.frame(LCI = quantile(boot_data$OR, probs = c(0.025)), row.names=NULL)
OR.UCI <- data.frame(LCI = quantile(boot_data$OR, probs = c(0.975)), row.names=NULL)
paste0(OR.LCI, OR.UCI)
boot.ci(boot.out = bootstraps,  type=c("norm")) #takes thew first value

boot.ci(boot.out = bootstraps, index=1, type=c("norm")) # takes specific values
boot.ci(boot.out = bootstraps, index=2, type=c("norm")) # takes specific values
boot.ci(boot.out = bootstraps, index=3, type=c("norm")) # takes specific values
boot.ci(boot.out = bootstraps, index=4, type=c("norm"))
# error - need to investigate
boot.ci(boot.out = bootstraps,index=3, type=c("bca"))

# Bootstrap function

boostrap.func <- function(intervention_data, i, vars, comparator_data, binary_var){
  # Samples the data
  bootstrap.data <- intervention_data[i,]

  # Performs weighing
  print(opt1 <- optim(par = rep(0,dim(bootstrap.data[, vars])[2]),
                      fn = objfn,
                      gr = gradfn,
                      X = as.matrix(bootstrap.data[, vars]),
                      method = "BFGS"))

  a1 <- opt1$par

  weights <- as.vector(exp(as.matrix(bootstrap.data[, vars]) %*% a1))

  bootstrap_dat_wt <- cbind(bootstrap.data, weights)

  #Combines intervention dataset with comparator_data
  analysis_data <- bootstrap_dat_wt %>%
    #select(ARM, Binary_event, weights) %>%
    rbind.fill(comparator_data)

  # binary data
  prop_A<-weighted.mean(bootstrap_dat_wt[,binary_var],bootstrap_dat_wt$weights)
  prop_B<-mean(comparator_data[,binary_var])
  odds_A <- weighted.mean(bootstrap_dat_wt[,binary_var],bootstrap_dat_wt$weights)/weighted.mean(1-bootstrap_dat_wt[,binary_var],bootstrap_dat_wt$weights)
  odds_B<-mean(comparator_data[,binary_var])/mean(1-comparator_data[,binary_var])
  logistic.regr <- suppressWarnings(glm(Binary_event~ARM, family=binomial(link="logit"), data = analysis_data, weights = weights))
  OR <- exp(as.numeric(coef(logistic.regr)["ARMComparator"]))
  OR_test <- odds_A/odds_B

  # survival data
  fit.A <- surv_fit(Surv(Time, Event) ~ 1, data = bootstrap_dat_wt, weights=weights)
  median_A<-surv_median(fit.A)$median
  #
  fit.B <- surv_fit(Surv(Time, Event) ~ 1, data = comparator_data)
  median_B<-surv_median(fit.B)$median

  cox_model <- coxph(Surv(Time, Event==1) ~ ARM, data = analysis_data, weights = weights)
  HR <- exp(cox_model$coefficients)

  c("OR"=OR, "prop_A"=prop_A, "prop_B"=prop_B,"HR"=HR, "median_A"=median_A, "median_B"=median_B)
}









# summaries ---------------------------------------------------------------

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Baseline summaries ------------------------------------------------

# To be included in the vignette, not as a function in the package

# Set weighted and unweighted intervention data
baseline_analysis_data <- intervention_wt_data %>% mutate(Treatment=paste0(ARM, "_matched")) %>%
  rbind(intervention_wt_data %>% mutate(Treatment=paste0(ARM, "_unadjusted"), wt = 1, wt_rs = 1)) %>%
  select(Treatment, est_weights$matching_vars, wt, wt_rs)

# Renames target population cols to match est_weights$matching_vars
est_weights$matching_vars
names(target_pop)
target_pop_standard <- target_pop %>%
  #EDIT
  rename(N=N,
         Treatment=ARM,
         AGE=age.mean,
         SEX=prop.male,
         SMOKE=prop.smoke,
         ECOG0=prop.ecog0
  ) %>%
  select(N, Treatment, est_weights$matching_vars)

# Summerises the baseline characteristics
Baseline_summary <- baseline_analysis_data %>%
  dplyr::group_by(Treatment) %>%
  summarise_each(funs(weighted.mean(., wt_rs)),-wt_rs) %>%
  rbind(target_pop_standard  %>% select(Treatment, est_weights$matching_vars) )

Baseline_summary_n <- baseline_analysis_data %>%
  dplyr::group_by(Treatment) %>%
  dplyr::summarise(
    'N' = n()) %>%
  rbind(target_pop_standard  %>% select(N, Treatment)) %>%
  rename(`N/ESS`=N)


Baseline_summary_all <- full_join(Baseline_summary_n, Baseline_summary, by="Treatment")

Baseline_summary_all2 <- cbind(Baseline_summary_all %>% select(-c(est_weights$matching_vars)),
                               lapply(Baseline_summary_all %>% select(est_weights$matching_vars), sprintf, fmt = "%.2f") %>% as.data.frame())


# replace N with ESS
Baseline_summary_all2$`N/ESS`[Baseline_summary_all2$Treatment == "B"] <- ESS



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Survival summaries ------------------------------------------------


