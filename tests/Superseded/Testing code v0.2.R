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
library(survival)



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
# Center baseline characteristics, estimate weights, rescaled weights and join data -----------------

#### center baseline characteristics
# (subtract the aggregate comparator data from the corresponding column of intervention PLD)
names(intervention_data)
intervention_data <- intervention_input %>%
                          mutate(Age_centered = AGE - target_pop$age.mean,
                                 # it is necessary to balance on both mean and standard deviation for continous variables:
                                 Age_squared_centered = (AGE^2) - (target_pop$age.mean^2 + target_pop$age.sd^2),
                                 Sex_centered = SEX - target_pop$prop.male,
                                 Smoke_centered = SMOKE - target_pop$prop.smoke,
                                 ECOG0_centered = ECOG0 - target_pop$prop.ecog0)
head(intervention_data)

# Set matching covariates
cent_match_cov <- c("Age_centered",
                    "Age_squared_centered",
                    "Sex_centered",
                    "Smoke_centered",
                    "ECOG0_centered")


match_cov <- c("AGE",
               "SEX",
               "SMOKE",
               "ECOG0")


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
estimate_weights <- function(intervention_data, cent_vars, comparator_data){

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
  data_with_wts <- cbind(intervention_data, wt, wt_rs) %>%
    mutate(ARM="Intervention")

  # assign weight=1 to comparator data
  comparator_data_wts <- comparator_data %>% mutate(wt=1, wt_rs=1, ARM="Comparator")

  # Join comparator data with the intervention data
  all_data <- rbind.fill(data_with_wts, comparator_data_wts)
  all_data$ARM <- relevel(as.factor(all_data$ARM), ref="Comparator")



  # Outputs are:
  #       - the analysis data (intervention PLD, weights and comparator pseudo PLD)
  #       - A charcacter vector with the name of the centered matching variables
  #       - A charcacter vector with the name of the matching variables
  output <- list(analysis_data = all_data,
                 centered_matching_vars = cent_vars,
                 intervention_wt_data=data_with_wts,
                 comparator_wt_data = comparator_data_wts
                 )

  return(output)

}


est_weights <- estimate_weights(intervention_data=intervention_data,
                                cent_vars = cent_match_cov,
                                comparator_data=comparator_input)

head(est_weights$analysis_data)
head(est_weights$intervention_wt_data)
head(est_weights$comparator_wt_data)
est_weights$centered_matching_vars


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Weight diagnostics ----------------------------------------------


#### ESS

# Function to caclulate ESS:
estimate_ess <- function(data, wt=wt){
  ess <- sum(data$wt)^2/sum(data$wt^2)
  return(ess)
}

ESS <- estimate_ess(data=est_weights$intervention_wt_data)
ESS


#### Weights summary

# Function to caclulate min, max, median, mean and sd of wts:
summarize_wts <- function(data, wt="wt"){
  summary <- data.frame(
      min = min(data[,wt]),
      max = max(data[,wt]),
      median = median(data[,wt]),
      mean = mean(data[,wt]),
      sd = sd(data[,wt])
    )
  return(summary)
}

weight_summ <- summarize_wts(data=est_weights$intervention_wt_data)
weight_summ

weight_rs_summ <- summarize_wts(data=est_weights$intervention_wt_data, wt="wt_rs")
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

wts_profile <- profile_wts(data=est_weights$intervention_wt_data, vars = match_cov)

# worth adding something like this?
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

histogram <- hist_wts(data=est_weights$intervention_wt_data)
histogram

histogram_bw <- hist_wts(data=est_weights$intervention_wt_data, bin_width=0.1)
histogram_bw



#### All weight diagnostics

# Function to combine weight diagnostic functions above:
all_wt_diagnostics <- function(data, # analysis data from estimate_weights
                               #arm,
                               matching_vars,
                               wt="wt",
                               ...){

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


diagnostics <- all_wt_diagnostics(data=est_weights$intervention_wt_data, matching_vars = match_cov, bin_width=0.1)

diagnostics$ESS
diagnostics$Summary_of_weights
diagnostics$Weight_profiles
diagnostics$Histogram_of_weights





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Bootstrapping ------------------------------------------------
# OR bootstraps-------------------------------------------------------------------------

boostrap_OR <- function(intervention_data, i,  cent_vars, comparator_data, binary_var){

  # Samples the data
  bootstrap_data <- intervention_data[i,]

  # Estimates weights
  perform_wt <- estimate_weights(intervention_data=bootstrap_data, cent_vars,  comparator_data)


  # binary data stats - note we can probably get rid of some of these, just to show that the OR works using a logitsic regression or manually shall we also add RR?
  prop_A<-weighted.mean(perform_wt$intervention_wt_data[,binary_var],perform_wt$intervention_wt_data$wt)
  prop_B<-mean(perform_wt$comparator_wt_data[,binary_var])
  logistic.regr_RR <- suppressWarnings(glm(Binary_event~ARM, family=poisson(link="log"), data = perform_wt$analysis_data, weight = wt))
  RR <- exp(as.numeric(coef(logistic.regr_RR)[2]))
  RR_test <- prop_A/prop_B

  odds_A <- weighted.mean(perform_wt$intervention_wt_data[,binary_var],perform_wt$intervention_wt_data$wt)/weighted.mean(1-perform_wt$intervention_wt_data[,binary_var],perform_wt$intervention_wt_data$wt)
  odds_B<-mean(perform_wt$comparator_wt_data[,binary_var])/mean(1-perform_wt$comparator_wt_data[,binary_var])
  logistic.regr <- suppressWarnings(glm(Binary_event~ARM, family=binomial(link="logit"), data = perform_wt$analysis_data, weight = wt))
  OR <- exp(as.numeric(coef(logistic.regr)[2]))
  OR_test <- odds_A/odds_B


  c("OR" = OR,  "RR" = RR, "odds_A"=odds_A, "odds_B"=odds_B, "prop_A" = prop_A, "prop_B" = prop_B, "OR_test"=OR_test, "RR_test" = RR_test)
}



# for Richard - to show how the function works
#test <- boostrap_OR(intervention_data=intervention_data, i=c(1:nrow(intervention_data)), cent_vars = cent_match_cov, comparator_data=comparator_input, binary_var="Binary_event")

OR_bootstraps <- boot(intervention_data, boostrap_OR, R=1000, cent_vars = cent_match_cov, comparator_data=comparator_input, binary_var="Binary_event")

# Bootstrap estimates
boot_data <- as.data.frame(OR_bootstraps$t)
colnames(boot_data) <- colnames(t(as.data.frame(OR_bootstraps$t0)))
head(boot_data)

# summerise bootstrap estimates
hist(boot_data$OR, main = "",xlab = "Boostrapped OR")
abline(v= quantile(boot_data$OR,probs = c(0.025,0.5,0.975)), lty=2)

OR.median <- quantile(boot_data$OR, probs = c(0.5))

OR.LCI <- quantile(boot_data$OR, probs = c(0.025))
OR.UCI <- quantile(boot_data$OR, probs = c(0.975))
paste0(OR.median, " (", OR.LCI, ",", OR.UCI, ")")

# Normal CI
boot.ci.OR <- boot.ci(boot.out = OR_bootstraps, index=1, type=c("norm")) # takes specific values
boot.ci.RR <- boot.ci(boot.out = OR_bootstraps, index=2, type=c("norm")) # takes specific values

# BCA CI
boot.ci.OR.BCA <- boot.ci(boot.out = OR_bootstraps, index=1, type=c("bca"))
boot.ci.RR.BCA <- boot.ci(boot.out = OR_bootstraps, index=2, type=c("bca"))

# Bootstrap CI function
boot.ci.OR$t0
boot.ci.OR$normal[2:3]

boot.ci.OR.BCA$t0
boot.ci.OR.BCA$bca[4:5]


# HR bootstraps -----------------------------------------------------------


boostrap_HR <- function(intervention_data, i,  cent_vars, comparator_data, binary_var){

    # Samples the data
    bootstrap_data <- intervention_data[i,]

    # Estimates weights
    perform_wt <- estimate_weights(intervention_data=bootstrap_data, cent_vars,  comparator_data)

    # survival data stats
    fit.A <- surv_fit(Surv(Time, Event) ~ 1, data = perform_wt$intervention_wt_data, weights=perform_wt$intervention_wt_data$wt)
    median_A<-surv_median(fit.A)$median
    #
    fit.B <- surv_fit(Surv(Time, Event) ~ 1, data = perform_wt$comparator_wt_data)
    median_B<-surv_median(fit.B)$median

    cox_model <- coxph(Surv(Time, Event==1) ~ ARM, data = perform_wt$analysis_data, weights = wt)
    HR <- exp(cox_model$coefficients)

    c("HR"=HR, "median_A"=median_A, "median_B"=median_B)
    }

  # for Richard - to show how the function works
  test_HR <- boostrap_HR(intervention_data=intervention_data, i=c(1:nrow(intervention_data)), cent_vars = cent_match_cov, comparator_data=comparator_input)

HR_bootstraps <- boot(intervention_data, boostrap_HR, R=1000, cent_vars = cent_match_cov, comparator_data=comparator_input, binary_var="Binary_event")

  # Bootstrap estimates
  boot_data <- as.data.frame(HR_bootstraps$t)
  colnames(boot_data) <- colnames(t(as.data.frame(HR_bootstraps$t0)))
  head(boot_data)

  # summerise bootstrap estimates
  hist(boot_data$HR, main = "",xlab = "Boostrapped HR")
  abline(v= quantile(boot_data$HR,probs = c(0.025,0.5,0.975)), lty=2)

  HR.median <- quantile(boot_data$HR, probs = c(0.5))

  HR.LCI <- quantile(boot_data$HR, probs = c(0.025))
  HR.UCI <- quantile(boot_data$HR, probs = c(0.975))
  paste0(HR.median, " (", HR.LCI, ",", HR.UCI, ")")

  # Normal CI
  boot.ci.HR <- boot.ci(boot.out = HR_bootstraps, index=1, type=c("norm")) # takes specific values
  boot.ci.RR <- boot.ci(boot.out = HR_bootstraps, index=2, type=c("norm")) # takes specific values
  # BCA CI

  boot.ci.HR.BCA <- boot.ci(boot.out = HR_bootstraps, index=1, type=c("bca"))
  boot.ci.RR.BCA <- boot.ci(boot.out = HR_bootstraps, index=2, type=c("bca"))

  # Bootstrap CI function
  boot.ci.HR$t0
  boot.ci.HR$normal[2:3]

  boot.ci.HR.BCA$t0
  boot.ci.HR.BCA$bca[4:5]




# summaries ---------------------------------------------------------------

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Baseline summaries ------------------------------------------------

# To be included in the vignette, not as a function in the package

# Set weighted and unweighted intervention data
baseline_analysis_data <- est_weights$intervention_wt_data %>% mutate(Treatment=paste0(ARM, "_matched")) %>%
  rbind(est_weights$intervention_wt_data %>% mutate(Treatment=paste0(ARM, "_unadjusted"), wt = 1, wt_rs = 1)) %>%
  select(Treatment, match_cov, wt, wt_rs)

# Renames target population cols to match match_cov
match_cov
names(target_pop)
target_pop_standard <- target_pop %>%
  #EDIT
  rename(N=N,
         Treatment="ARM",
         AGE=age.mean,
         SEX=prop.male,
         SMOKE=prop.smoke,
         ECOG0=prop.ecog0
  ) %>%
  select(N, Treatment, match_cov)

# Summerises the baseline characteristics
Baseline_summary <- baseline_analysis_data %>%
  dplyr::group_by(Treatment) %>%
  summarise_each(list(~ weighted.mean(., wt)),-c(wt,wt_rs)) %>%
  rbind(target_pop_standard  %>% select(Treatment, match_cov) )

Baseline_summary_n <- baseline_analysis_data %>%
  dplyr::group_by(Treatment) %>%
  dplyr::summarise(
    'N' = n()) %>%
  rbind(target_pop_standard  %>% select(N, Treatment)) %>%
  rename(`N/ESS`=N)


Baseline_summary_all <- full_join(Baseline_summary_n, Baseline_summary, by="Treatment")

Baseline_summary <- cbind(Baseline_summary_all %>% select(-c(match_cov)),
                               lapply(Baseline_summary_all %>% select(match_cov), sprintf, fmt = "%.2f") %>% as.data.frame())


# replace N with ESS
Baseline_summary$`N/ESS`[Baseline_summary$Treatment == "Intervention_unadjusted"] <- ESS
Baseline_summary


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Survival summaries ------------------------------------------------


## HRs
unweighted.cox <- coxph(Surv(Time, Event==1) ~ ARM, data = est_weights$analysis_data)
weighted.cox <- coxph(Surv(Time, Event==1) ~ ARM, data = est_weights$analysis_data, weights = wt)


cox.summ <- rbind(summary(unweighted.cox)$conf.int, summary(weighted.cox)$conf.int) %>%
  as.data.frame() %>%
  select(-`exp(-coef)`) %>% #drop unnecessary column
  rename(HR = `exp(coef)`, HR.low.CI = `lower .95`, HR.upp.CI = `upper .95`) %>%
  mutate(Method = c("Unadjusted", "Cox weighted")) %>%
  rbind(data.frame("HR" = HR.median, "HR.low.CI" = boot.ci.HR$normal[2], "HR.upp.CI" = boot.ci.HR$normal[3], "Method"="Normal bootstrap")) %>%
 rbind(data.frame("HR" = HR.median, "HR.low.CI" = boot.ci.HR.BCA$bca[4], "HR.upp.CI" = boot.ci.HR.BCA$bca[5], "Method"="BCA bootstrap")) %>%
  mutate(HR.95.CI = paste0(sprintf('%.3f', HR), " (", sprintf('%.3f', HR.low.CI), ", ", sprintf('%.3f', HR.upp.CI), ")")) %>%
  select(Method, HR.95.CI)
cox.summ

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Binary summaries ------------------------------------------------

logistic.regr_OR <- suppressWarnings(glm(Binary_event~ARM, family=binomial(link="logit"), data = est_weights$analysis_data, weight = wt))
exp(as.numeric(coef(logistic.regr_OR)[2]))
summary(logistic.regr_OR)
