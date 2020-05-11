# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 2797 exploratory - testing out the code                         #
# Author: SS/ JG (23.04.2020)                                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Load libraries --------------------------------------------------
library(plyr)
library(dplyr)
# library(haven)
# library(ggplot2)
# library(tidyr)
library(boot)
# library(survminer)
library(survival)
library(MAIC)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Directories -----------------------------------------------------
base_dir <- 'G:/Clients/Roche/2797 Development of R Code for MAIC and mixture cure models/Project'
data_path <- file.path(base_dir,'2 Exploratory/Simulated datasets')



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Load data -------------------------------------------------------

#### Intervention data

# Read in ADaM data and rename variables of interest
adsl <- read.csv(file.path(data_path, "adsl.csv")) %>% # subject level data
  dplyr::mutate(SEX=ifelse(SEX=="Male",1,0))

adrs <- read.csv(file.path(data_path, "adrs.csv")) %>% # response data
  dplyr::filter(PARAM=="Response") %>%
  dplyr::select(USUBJID, ARM, Binary_event=AVAL) ## change "Binary_event" to "response"?

adtte <- read.csv(file.path(data_path, "adtte.csv")) %>% # time to event data
  dplyr::filter(PARAMCD=="OS") %>%
  dplyr::mutate(Event=1-CNSR) %>%
  dplyr::select(USUBJID, ARM, Time=AVAL, Event)

# Combine all intervention data
intervention_input <- plyr::join_all(list(adsl, adrs, adtte), type = "full", by=c("USUBJID", "ARM"))

#### Comparator pseudo data

# read in digitised pseudo survival data
comparator_surv <- read.csv(file.path(data_path,"psuedo_IPD.csv"))

# simulate response data based on the known proportion of responders
comparator_n <- nrow(comparator_surv) # total number of patients in the comparator data
comparator_prop_events <- 0.4 # proportion of responders
comparator_binary <- data.frame("Binary_event"=
                                  c(rep(1,comparator_n*comparator_prop_events),
                                    rep(0, comparator_n*(1-comparator_prop_events))))

# join survival and response comparator data note not a 1:1 relationship - the
# rows do not represent the same observation
comparator_input <- cbind(comparator_surv, comparator_binary)

# Baseline aggregate data for the comparator population
target_pop <- read.csv(file.path(data_path,"Aggregate data.csv"))
target_pop


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Center baseline characteristics -----------------

#### center baseline characteristics
# (subtract the aggregate comparator data from the corresponding column of intervention PLD)
names(intervention_input)
intervention_data <- intervention_input %>%
                          dplyr::mutate(Age_centered = AGE - target_pop$age.mean,
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

# Renames target population cols to match match_cov
match_cov
names(target_pop)
target_pop_standard <- target_pop %>%
  #EDIT
  dplyr::rename(N=N,
                Treatment="ARM",
                AGE=age.mean,
                SEX=prop.male,
                SMOKE=prop.smoke,
                ECOG0=prop.ecog0
  ) %>%
  dplyr::select(N, Treatment, match_cov)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  estimate weights, rescaled weights  -----------------

est_weights <- estimate_weights(intervention_data=intervention_data,
                                comparator_data=comparator_input,
                                matching_vars = cent_match_cov)

head(est_weights$analysis_data)
est_weights$matching_vars

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Weight diagnostics ----------------------------------------------
#### ESS
#extract data for the intervention arm and calculate ESS
ESS <- filter(est_weights$analysis_data, ARM == 'Intervention') %>%
  estimate_ess(wt_col = 'wt_rs')
ESS

#### Weights summary
#extract data for the intervention arm and calculate summary statistics for the weights
weight_summ <- filter(est_weights$analysis_data, ARM == 'Intervention') %>%
  summarize_wts()
weight_summ

#extract data for the intervention arm and calculate summary statistics for the rescaled weights
weight_rs_summ <- filter(est_weights$analysis_data, ARM == 'Intervention') %>%
  summarize_wts(wt_col="wt_rs")
weight_rs_summ

##### Weight profiles
wts_profile <- filter(est_weights$analysis_data, ARM == 'Intervention') %>%
  profile_wts(vars = match_cov)
head(wts_profile)

# worth adding something like this?
plot(wts_profile$AGE, wts_profile$wt)
boxplot(wts_profile$SEX , wts_profile$wt)

#### Histograms
# Plot histograms of unscaled and rescaled weights
# bin_width needs to be adapted depending on the sample size in the data set
histogram <- filter(est_weights$analysis_data, ARM == 'Intervention') %>%
  hist_wts(bin = 50)
histogram



#### All weight diagnostics

# Function to produce a set of diagnostics.
# Calls each of the diagnostic functions above except for plotting histograms
diagnostics <- filter(est_weights$analysis_data, ARM == 'Intervention') %>%
  wt_diagnostics(vars = est_weights$matching_vars)

diagnostics$ESS
diagnostics$Summary_of_weights
diagnostics$Weight_profiles


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Baseline summaries ------------------------------------------------

# To be included in the vignette, not as a function in the package

# Set weighted and unweighted intervention data
baseline_analysis_data <- est_weights$analysis_data %>% filter(ARM=="Intervention") %>% dplyr::mutate(Treatment=paste0(ARM, "_matched")) %>%
  rbind(est_weights$analysis_data %>% filter(ARM=="Intervention")%>% mutate(Treatment=paste0(ARM, "_unadjusted"), wt = 1, wt_rs = 1)) %>%
  select(Treatment, match_cov, wt, wt_rs)

# Summerises the baseline characteristics
Baseline_summary <- baseline_analysis_data %>%
  dplyr::group_by(Treatment) %>%
  dplyr::summarise_each(list(~ weighted.mean(., wt)),-c(wt,wt_rs)) %>%
  rbind(target_pop_standard  %>% select(Treatment, match_cov) )

Baseline_summary_n <- baseline_analysis_data %>%
  dplyr::group_by(Treatment) %>%
  dplyr::summarise(
    'N' = n()) %>%
  rbind(target_pop_standard  %>% select(N, Treatment)) %>%
  dplyr::rename(`N/ESS`=N)


Baseline_summary_all <- dplyr::full_join(Baseline_summary_n, Baseline_summary, by="Treatment")

Baseline_summary <- cbind(Baseline_summary_all %>% select(-c(match_cov)),
                          lapply(Baseline_summary_all %>% select(match_cov), sprintf, fmt = "%.2f") %>% as.data.frame())


# replace N with ESS
Baseline_summary$`N/ESS`[Baseline_summary$Treatment == "Intervention_unadjusted"] <- ESS
Baseline_summary


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Bootstrapping ------------------------------------------------
# OR bootstraps-------------------------------------------------------------------------

# Demonstrate functionality of the bootstrap_OR function
# This function returns a single estimate of the odds ratio
# This function is intended to be used in conjunction with the boot function, not called directly by the user
int <- filter(est_weights$analysis_data, ARM == 'Intervention')
comp <- filter(est_weights$analysis_data, ARM == 'Comparator')
OR <- bootstrap_OR(intervention_data=int, comparator_data=comp, matching = est_weights$matching_vars,
                  i=c(1:nrow(intervention_data)), model = 'Binary_event ~ ARM')




# Bootstrap estimates
OR_bootstraps <- boot(data = int, statistic = bootstrap_OR, R=1000, comparator_data=comp,
                      matching = est_weights$matching_vars, model = 'Binary_event ~ ARM'
                      )

# summarise bootstrap estimates
hist(OR_bootstraps$t, main = "", xlab = "Boostrapped OR")
abline(v= quantile(OR_bootstraps$t, probs = c(0.025,0.5,0.975)), lty=2)

OR_median <- median(OR_bootstraps$t)

# Bootstrap CI function - Normal CI
boot_ci_OR <- boot.ci(boot.out = OR_bootstraps, index=1, type="norm") # takes specific values

# Bootstrap CI function - BCA CI
boot_ci_OR_BCA <- boot.ci(boot.out = OR_bootstraps, index=1, type="bca")

# ORs
unweighted_OR <- suppressWarnings(glm(Binary_event~ARM, family=binomial(link="logit"), data = est_weights$analysis_data))
weighted_OR <- suppressWarnings(glm(Binary_event~ARM, family=binomial(link="logit"), data = est_weights$analysis_data, weight = wt))


OR_summ <- rbind(exp(cbind("Odds ratio" = coef(unweighted_OR), confint.default(unweighted_OR, level = 0.95)))[2,],
                 exp(cbind("Odds ratio" = coef(weighted_OR), confint.default(weighted_OR, level = 0.95)))[2,]) %>%
  as.data.frame() %>%
  dplyr::rename(OR = `Odds ratio`, OR_low_CI = `2.5 %`, OR_upp_CI = `97.5 %`) %>%
  mutate(Method = c("Unadjusted", "Weighted")) %>%
  rbind(data.frame("OR" = OR_median, "OR_low_CI" = boot_ci_OR$normal[2], "OR_upp_CI" = boot_ci_OR$normal[3], "Method"="Normal bootstrap")) %>%
  rbind(data.frame("OR" = OR_median, "OR_low_CI" = boot_ci_OR_BCA$bca[4], "OR_upp_CI" = boot_ci_OR_BCA$bca[5], "Method"="BCA bootstrap")) %>%
  dplyr::mutate(OR_95_CI = paste0(sprintf('%.3f', OR), " (", sprintf('%.3f', OR_low_CI), ", ", sprintf('%.3f', OR_upp_CI), ")")) %>%
  dplyr::select(Method, OR_95_CI)
OR_summ



# HR bootstraps -----------------------------------------------------------
# Demonstrate functionality of the bootstrap_HR function
# This function returns a single estimate of the hazard ratio
# This function is intended to be used in conjunction with the boot function, not called directly by the user
HR <- bootstrap_HR(intervention_data=int, comparator_data=comp, matching = est_weights$matching_vars,
                  i=c(1:nrow(intervention_data)), model = Surv(Time, Event==1) ~ ARM
                  )

# Bootstrap estimates
HR_bootstraps <- boot(data = int, statistic = bootstrap_HR, R=1000, comparator_data=comp,
                      matching = est_weights$matching_vars, model = Surv(Time, Event==1) ~ ARM
                      )

# Summarise bootstrap estimates
hist(HR_bootstraps$t, main = "", xlab = "Boostrapped HR")
abline(v= quantile(HR_bootstraps$t, probs = c(0.025, 0.5, 0.975)), lty=2)

HR_median <- median(HR_bootstraps$t)
# Bootstrap CI function - Normal CI
boot_ci_HR <- boot.ci(boot.out = HR_bootstraps, index=1, type="norm") # takes specific values

# Bootstrap CI function - BCA CI
boot_ci_HR_BCA <- boot.ci(boot.out = HR_bootstraps, index=1, type="bca")


## HRs
unweighted_cox <- coxph(Surv(Time, Event==1) ~ ARM, data = est_weights$analysis_data)
weighted_cox <- coxph(Surv(Time, Event==1) ~ ARM, data = est_weights$analysis_data, weights = wt)


HR_summ <- rbind(summary(unweighted_cox)$conf.int, summary(weighted_cox)$conf.int) %>%
  as.data.frame() %>%
  dplyr::select(-`exp(-coef)`) %>% #drop unnecessary column
  dplyr::rename(HR = `exp(coef)`, HR_low_CI = `lower .95`, HR_upp_CI = `upper .95`) %>%
  mutate(Method = c("Unadjusted", "Cox weighted")) %>%
  rbind(data.frame("HR" = HR_median, "HR_low_CI" = boot_ci_HR$normal[2], "HR_upp_CI" = boot_ci_HR$normal[3], "Method"="Normal bootstrap")) %>%
  rbind(data.frame("HR" = HR_median, "HR_low_CI" = boot_ci_HR_BCA$bca[4], "HR_upp_CI" = boot_ci_HR_BCA$bca[5], "Method"="BCA bootstrap")) %>%
  dplyr::mutate(HR_95_CI = paste0(sprintf('%.3f', HR), " (", sprintf('%.3f', HR_low_CI), ", ", sprintf('%.3f', HR_upp_CI), ")")) %>%
  dplyr::select(Method, HR_95_CI)
HR_summ





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Binary summaries ------------------------------------------------



