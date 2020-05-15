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
library(ggplot2)
library(survminer)



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

intervention_input <- full_join(full_join(adsl,adrs, by=c("USUBJID", "ARM")), adtte, by=c("USUBJID", "ARM"))

head(intervention_input)

# List out matching covariates
match_cov <- c("AGE",
               "SEX",
               "SMOKE",
               "ECOG0")

# read in digitised pseudo survival data
comparator_surv <- read.csv(file.path(data_path,"psuedo_IPD.csv"))

# Simulate response data based on the known proportion of responders
comparator_n <- nrow(comparator_surv) # total number of patients in the comparator data
comparator_prop_events <- 0.4 # proportion of responders
comparator_binary <- data.frame("Binary_event"=
                                  c(rep(1,comparator_n*comparator_prop_events),
                                    rep(0, comparator_n*(1-comparator_prop_events))))

# Join survival and response comparator data
# (note the rows do not represent observations from a particular patient)
comparator_input <- cbind(comparator_surv, comparator_binary)

# Baseline aggregate data for the comparator population
target_pop <- read.csv(file.path(data_path,"Aggregate data.csv"))

# Renames target population cols to be consistent with match_cov
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
  dplyr::select(N, Treatment, all_of(match_cov))

names(intervention_input)
intervention_data <-
  intervention_input %>%
  dplyr::mutate(Age_centered = AGE - target_pop$age.mean,
                # matching on both mean and standard deviation for continuous variables:
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

est_weights <- estimate_weights(intervention_data=intervention_data,
                                comparator_data=comparator_input,
                                matching_vars = cent_match_cov)

head(est_weights$analysis_data)

est_weights$matching_vars

# Function to produce a set of diagnostics.
# Calls each of the diagnostic functions above except for plotting histograms
diagnostics <- filter(est_weights$analysis_data, ARM == 'Intervention') %>%
  wt_diagnostics(vars = est_weights$matching_vars)

ESS <- diagnostics$ESS
diagnostics$Summary_of_weights
head(diagnostics$Weight_profiles)

# Plot histograms of unscaled and rescaled weights
# bin_width needs to be adapted depending on the sample size in the data set
histogram <- est_weights$analysis_data %>%
  filter(ARM == 'Intervention') %>%
  hist_wts(bin = 50)
histogram

baseline_summary <- list('Intervention_weighted' = NA, 'Intervention' = NA, 'Comparator' = NA)

# Summarise matching variables for weighted intervention data
baseline_summary$Intervention_weighted <- est_weights$analysis_data %>%
  filter(ARM=="Intervention") %>%
  dplyr::select(all_of(match_cov), wt) %>%
  dplyr::summarise_at(match_cov, list(~ weighted.mean(., wt)))

# Summarise matching variables for unweighted intervention data
baseline_summary$Intervention <- est_weights$analysis_data %>%
  filter(ARM=="Intervention") %>%
  dplyr::select(all_of(match_cov), wt) %>%
  dplyr::summarise_at(match_cov, list(~ mean(.)))

# baseline data for the comparator study
baseline_summary$Comparator <- select(target_pop_standard, all_of(match_cov))

# Combine the three summaries
# Takes a list of data frames and binds these together
trt <- names(baseline_summary)
baseline_summary <- dplyr::bind_rows(baseline_summary) %>%
  transmute_all(sprintf, fmt = "%.2f") %>% #apply rounding for presentation
  mutate(ARM = as.character(trt)) %>% #Add treatment labels
  select(ARM, all_of(match_cov))

# Count the number of patients in the unweighted data
summary_n <- est_weights$analysis_data %>%
  dplyr::group_by(ARM) %>%
  dplyr::summarise('N' = n()) %>%
  dplyr::mutate(ARM = as.character(ARM))

# Combine sample sizes with summary data
baseline_summary <- full_join(baseline_summary, summary_n) %>%
  select(ARM, `N/ESS` = N, all_of(match_cov))

# Insert the ESS as the sample size for the weighted data
# This is calculated above but can also be obtained using the estimate_ess function as shown below
baseline_summary$`N/ESS`[baseline_summary$ARM == "Intervention_weighted"] <- est_weights$analysis_data %>%
  filter(ARM == 'Intervention') %>%
  estimate_ess(wt_col = 'wt')
baseline_summary

# Fit a Cox model without weights to estimate the unweighted HR
unweighted_cox <- coxph(Surv(Time, Event==1) ~ ARM, data = est_weights$analysis_data)

HR_CI_cox <- summary(unweighted_cox)$conf.int %>%
  as.data.frame() %>%
  dplyr::select(-`exp(-coef)`)
HR_CI_cox

# Fit a Cox model with weights to estimate the weighted HR
weighted_cox <- coxph(Surv(Time, Event==1) ~ ARM, data = est_weights$analysis_data, weights = wt)

HR_CI_cox_wtd <- summary(weighted_cox)$conf.int %>%
  as.data.frame() %>%
  dplyr::select(-`exp(-coef)`)
HR_CI_cox_wtd

## Bootstrapping

# Separate data into intervention data and comparator data
# The boostrap_HR function below makes use of the estimate_weights fucntion
# which requires separate datasets
int <- filter(est_weights$analysis_data, ARM == 'Intervention')
comp <- filter(est_weights$analysis_data, ARM == 'Comparator')

# Bootstrap 1000 HRs
HR_bootstraps <- boot(data = int, # intervention data
                      statistic = bootstrap_HR, # bootstrap the HR (defined in the MAIC package)
                      R=1000, # number of bootstrap samples
                      comparator_data = comp, # comparator pseudo data
                      matching = est_weights$matching_vars, # matching variables
                      model = Surv(Time, Event==1) ~ ARM # model to fit
)

# Summarize bootstrap estimates in a histogram
# Vertical lines indicate the median and upper and lower CIs
hist(HR_bootstraps$t, main = "", xlab = "Boostrapped HR")
abline(v= quantile(HR_bootstraps$t, probs = c(0.025, 0.5, 0.975)), lty=2)

# Median of the bootstrap samples
HR_median <- median(HR_bootstraps$t)

# Bootstrap CI - Percentile CI
boot_ci_HR <- boot.ci(boot.out = HR_bootstraps, index=1, type="norm")

# Bootstrap CI - BCa CI
boot_ci_HR_BCA <- boot.ci(boot.out = HR_bootstraps, index=1, type="bca")

## Summary

# Produce a summary of HRs and CIs
HR_summ <-  rbind(HR_CI_cox, HR_CI_cox_wtd) %>% # Unweighted and weights HRs and CIs from Cox models
  dplyr::rename(HR = `exp(coef)`,
                HR_low_CI = `lower .95`,
                HR_upp_CI = `upper .95`) %>%
  mutate(Method = c("HR (95% CI) from unadjusted Cox model",
                    "HR (95% CI) from weighted Cox model")) %>%

  # Median bootstrapped HR and 95% percentile CI
  rbind(data.frame("HR" = HR_median,
                   "HR_low_CI" = boot_ci_HR$normal[2],
                   "HR_upp_CI" = boot_ci_HR$normal[3],
                   "Method"="Bootstrap median HR (95% percentile CI)")) %>%

  # Median bootstrapped HR and 95% bias-corrected and accelerated bootstrap CI
  rbind(data.frame("HR" = HR_median,
                   "HR_low_CI" = boot_ci_HR_BCA$bca[4],
                   "HR_upp_CI" = boot_ci_HR_BCA$bca[5],
                   "Method"="Bootstrap median HR (95% BCa CI)")) %>%

  # Format HR and CI in one variable, rounded to 3 decimal places
  dplyr::mutate(HR_95_CI = paste0(sprintf('%.3f', HR),
                                  " (",
                                  sprintf('%.3f', HR_low_CI),
                                  ", ",
                                  sprintf('%.3f', HR_upp_CI),
                                  ")")
  ) %>%
  dplyr::select(Method, 'HR (95% CI)' = HR_95_CI)
HR_summ

HR_table <- HR_summ %>%
  regulartable() %>% #make it a flextable object
  set_header_labels(Method = "Method",  HR_95_CI = "Hazard ratio (95% CI)")  %>%
  font(font = 'Arial', part = 'all') %>%
  fontsize(size = 10, part = 'all') %>%
  bold(part = 'header') %>%
  align(align = 'center', part = 'all') %>%
  align(j = 1, align = 'left', part = 'all') %>%
  border_outer(border = fp_border()) %>%
  border_inner_h(border = fp_border()) %>%
  border_inner_v(border = fp_border()) %>%
  autofit(add_w = 0.2)
HR_table

# Unweighted intervention data
KM_int <- survfit(formula = Surv(Time, Event==1) ~ 1 ,
                  data = int,
                  type="kaplan-meier")
# Weighted intervention data
KM_int_wtd <- survfit(formula = Surv(Time, Event==1) ~ 1 ,
                      data = int,
                      weights = wt,
                      type="kaplan-meier")
# Comparator data
KM_comp <- survfit(formula = Surv(Time, Event==1) ~ 1 ,
                   data = comp,
                   type="kaplan-meier")

# Combine the survfit objects ready for ggsurvplot
KM_list <- list(Intervention = KM_int,
                Intervention_weighted = KM_int_wtd,
                Comparator = KM_comp)

# Produce the Kaplan-Meier plot


KM_plot <- ggsurvplot(KM_list,
                      combine = TRUE,
                      risk.table=T)

KM_plot <- ggsurvplot(KM_int)


KM_plot <- ggsurvplot(KM_list,
                      combine = TRUE,
                      risk.table=T, # numbers at risk displayed on the plot
                      break.x.by=50,
                      xlab="Time (days)",
                      legend.title = "Treatment",
                      title = "Kaplan-Meier plot of overall survival",
                      legend.labs=c("Intervention", "Intervention weighted", "Comparator"),
                      font.legend = list(size = 10)) +
  guides(colour=guide_legend(nrow = 2))
KM_plot






km.plot <- ggsurvplot(KM_list, combine = TRUE, risk.table = TRUE,
                      break.time.by = 50,
                      conf.int = FALSE,
                      censor=FALSE,
                      legend.title = '',
                      xlab = paste0('Time (days)'),
                       font.x = 16,
                      font.y = 18,
                      font.legend = 16,
                      font.xtickslab = 16,
                      font.ytickslab = 16,
                      fontsize = 6)
km.plot$table <- ggpar(
  km.plot$table,
  font.x        = c(16),
  font.xtickslab = c(16),
  font.ytickslab = c(16)
)




## Calculate ORs

# Fit a logistic regression model without weights to estimate the unweighted OR
unweighted_OR <- glm(formula = Binary_event~ARM,
                     family = binomial(link="logit"),
                     data = est_weights$analysis_data)

# Log odds ratio
log_OR_CI_logit <- cbind("Log odds ratio" = coef(unweighted_OR),
                         confint.default(unweighted_OR, level = 0.95))[2,]

# Odds ratio
OR_CI_logit <- exp(cbind("Odds ratio" = coef(unweighted_OR),
                         confint.default(unweighted_OR, level = 0.95)))[2,]
OR_CI_logit

# Fit a logistic regression model with weights to estimate the weighted OR
weighted_OR <- suppressWarnings(glm(formula = Binary_event~ARM,
                                    family = binomial(link="logit"),
                                    data = est_weights$analysis_data,
                                    weight = wt))

# Weighted log odds ratio
log_OR_CI_logit_wtd <- cbind("Log odds ratio" = coef(weighted_OR),
                             confint.default(weighted_OR, level = 0.95))[2,]

# Weighted odds ratio
OR_CI_logit_wtd <- exp(cbind("Odds ratio" = coef(weighted_OR),
                             confint.default(weighted_OR, level = 0.95)))[2,]
OR_CI_logit_wtd


## Bootstrapping

# Separate data into intervention data and comparator data
# The boostrap_OR function below makes use of the estimate_weights fucntion
# which requires separate datasets
int <- filter(est_weights$analysis_data, ARM == 'Intervention')
comp <- filter(est_weights$analysis_data, ARM == 'Comparator')

# Bootstrap 1000 ORs
OR_bootstraps <- boot(data = int, # intervention data
                      statistic = bootstrap_OR, # bootstrap the OR
                      R = 1000, # number of bootstrap samples
                      comparator_data = comp, # comparator pseudo data
                      matching = est_weights$matching_vars, # matching variables
                      model = 'Binary_event ~ ARM' # model to fit
)

# Summarize bootstrap estimates in a histogram
# Vertical lines indicate the median and upper and lower CIs
hist(OR_bootstraps$t, main = "", xlab = "Boostrapped OR")
abline(v= quantile(OR_bootstraps$t, probs = c(0.025,0.5,0.975)), lty=2)

# Median of the bootstrap samples
OR_median <- median(OR_bootstraps$t)

# Bootstrap CI - Percentile CI
boot_ci_OR <- boot.ci(boot.out = OR_bootstraps, index=1, type="norm")

# Bootstrap CI - BCa CI
boot_ci_OR_BCA <- boot.ci(boot.out = OR_bootstraps, index=1, type="bca")


## Summary

# Produce summary of ORs and CIs
OR_summ <-  rbind(OR_CI_logit, OR_CI_logit_wtd) %>% # Unweighted and weighted ORs and CIs from logistic regression models
  as.data.frame() %>%
  dplyr::rename(OR = `Odds ratio`, OR_low_CI = `2.5 %`, OR_upp_CI = `97.5 %`) %>%
  mutate(Method = c("OR (95% CI) from unadjusted logistic regression model",
                    "OR (95% CI) from weighted logistic regression model")) %>%

  # Median bootstrapped HR and 95% percentile CI
  rbind(data.frame("OR" = OR_median,
                   "OR_low_CI" = boot_ci_OR$normal[2],
                   "OR_upp_CI" = boot_ci_OR$normal[3],
                   "Method"="Bootstrap median HR (95% percentile CI)")) %>%

  # Median bootstrapped HR and 95% bias-corrected and accelerated bootstrap CI
  rbind(data.frame("OR" = OR_median,
                   "OR_low_CI" = boot_ci_OR_BCA$bca[4],
                   "OR_upp_CI" = boot_ci_OR_BCA$bca[5],
                   "Method"="Bootstrap median HR (95% BCa CI)")) %>%

  # Format OR and CI in one variable, rounded to 3 decimal places
  dplyr::mutate(OR_95_CI = paste0(sprintf('%.3f', OR),
                                  " (",
                                  sprintf('%.3f', OR_low_CI),
                                  ", ",
                                  sprintf('%.3f', OR_upp_CI),
                                  ")")
  ) %>%

  dplyr::select(Method, OR_95_CI )

OR_summ2 <- OR_summ %>%
    regulartable() %>% #make it a flextable object
  set_header_labels(Method = "Method",  HR.95.CI = "Hazard ratio (95% CI)")  %>%
    font(font = 'Arial', part = 'all') %>%
  fontsize(size = 10, part = 'all') %>%
  bold(part = 'header') %>%
  align(align = 'center', part = 'all') %>%
  align(j = 1, align = 'left', part = 'all') %>%
  border_outer(border = fp_border()) %>%
  border_inner_h(border = fp_border()) %>%
  border_inner_v(border = fp_border()) %>%
  autofit(add_w = 0.2)
OR_summ2

2

