## ---- include = FALSE-----------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 6
)


## ---- warning = FALSE, message = FALSE------------------------------------------------------------------------------------------------------------------------
library(dplyr)
library(boot)
library(survival)
library(MAIC)
library(ggplot2)
library(survminer)
library(officer)
library(flextable)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
#### Intervention data

# Read in ADaM data and rename variables of interest

adsl <- read.csv(system.file("extdata", "adsl.csv", package = "MAIC", mustWork = TRUE))
adrs <- read.csv(system.file("extdata", "adrs.csv", package = "MAIC", mustWork = TRUE))
adtte <- read.csv(system.file("extdata", "adtte.csv", package = "MAIC", mustWork = TRUE))


adsl <- adsl %>% # Data containing the matching variables
  dplyr::mutate(SEX=ifelse(SEX=="Male", 1, 0)) # Coded 1 for males and 0 for females

adrs <- adrs %>% # Response data
  dplyr::filter(PARAM=="Response") %>%
  dplyr::transmute(USUBJID, ARM, response=AVAL)

adtte <- adtte %>% # Time to event data (overall survival)
  dplyr::filter(PARAMCD=="OS") %>%
  dplyr::mutate(Event=1-CNSR) %>%
  dplyr::transmute(USUBJID, ARM, Time=AVAL, Event)

# Combine all intervention data
intervention_input <- full_join(full_join(adsl,adrs, by=c("USUBJID", "ARM")), adtte,
                                by=c("USUBJID", "ARM"))
head(intervention_input)

# List out matching covariates
match_cov <- c("SEX")



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
#### Comparator pseudo data

# Read in digitised pseudo survival data, col names must match intervention_input
comparator_surv <- read.csv(system.file("extdata", "psuedo_IPD.csv",
                                        package = "MAIC", mustWork = TRUE)) %>%
  rename(Time=Time, Event=Event)


# Simulate response data based on the known proportion of responders
comparator_n <- nrow(comparator_surv) # total number of patients in the comparator data
comparator_prop_events <- 0.4 # proportion of responders
comparator_binary <- data.frame("response"=
                                  c(rep(1,comparator_n*comparator_prop_events),
                                    rep(0, comparator_n*(1-comparator_prop_events))))

# Join survival and response comparator data
# (note the rows do not represent observations from a particular patient)
comparator_input <- cbind(comparator_surv, comparator_binary)




## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Baseline aggregate data for the comparator population
target_pop <- read.csv(system.file("extdata", "Aggregate data.csv",
                                   package = "MAIC", mustWork = TRUE))

# Renames target population cols to be consistent with match_cov
match_cov
names(target_pop)
target_pop_standard <- target_pop %>%
  #EDIT
  dplyr::rename(N=N,
                Treatment=ARM,
                SEX=prop.male
  ) %>%
  dplyr::transmute(N, Treatment, SEX)



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
#### center baseline characteristics
# (subtract the aggregate comparator data from the corresponding column of intervention PLD)
names(intervention_input)
intervention_data <- intervention_input %>%
        dplyr::mutate(
                      Sex_centered = SEX - target_pop$prop.male)
head(intervention_data)

# Set matching covariates
cent_match_cov <- c("Sex_centered")




## -------------------------------------------------------------------------------------------------------------------------------------------------------------
est_weights <- estimate_weights(intervention_data=intervention_data,
                                comparator_data=comparator_input,
                                matching_vars = cent_match_cov)

head(est_weights$analysis_data)

est_weights$matching_vars




## -------------------------------------------------------------------------------------------------------------------------------------------------------------
ESS <- filter(est_weights$analysis_data, ARM == 'Intervention') %>%
  estimate_ess()
ESS



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot histograms of unscaled and rescaled weights
# bin_width needs to be adapted depending on the sample size in the data set
histogram <- est_weights$analysis_data %>%
              filter(ARM == 'Intervention') %>%
              hist_wts(bin = 50)
histogram



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
weight_summ <- filter(est_weights$analysis_data, ARM == 'Intervention') %>%
                summarize_wts()
weight_summ



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
wts_profile <- filter(est_weights$analysis_data, ARM == 'Intervention') %>%
                profile_wts(vars = match_cov)
head(wts_profile)



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function to produce a set of diagnostics.
# Calls each of the diagnostic functions above except for plotting histograms
diagnostics <- filter(est_weights$analysis_data, ARM == 'Intervention') %>%
                wt_diagnostics(vars = est_weights$matching_vars)

diagnostics$ESS

diagnostics$Summary_of_weights
head(diagnostics$Weight_profiles)




## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create an object to hold the output
baseline_summary <- list('Intervention_weighted' = NA, 'Intervention' = NA, 'Comparator' = NA)

# Summarise matching variables for weighted intervention data
baseline_summary$Intervention_weighted <- est_weights$analysis_data %>%
  filter(ARM=="Intervention") %>%
  dplyr::transmute(SEX, wt) %>%
  dplyr::summarise_at(match_cov, list(~ weighted.mean(., wt)))

# Summarise matching variables for unweighted intervention data
baseline_summary$Intervention <- est_weights$analysis_data %>%
  filter(ARM=="Intervention") %>%
  dplyr::transmute(SEX, wt) %>%
  dplyr::summarise_at(match_cov, list(~ mean(.)))

# baseline data for the comparator study
baseline_summary$Comparator <- transmute(target_pop_standard, SEX)

# Combine the three summaries
# Takes a list of data frames and binds these together
trt <- names(baseline_summary)
baseline_summary <- dplyr::bind_rows(baseline_summary) %>%
  transmute_all(sprintf, fmt = "%.2f") %>% #apply rounding for presentation
  transmute(ARM = as.character(trt), SEX)

# Count the number of patients in the unweighted data
summary_n <- est_weights$analysis_data %>%
  dplyr::group_by(ARM) %>%
  dplyr::summarise('N' = n()) %>%
  dplyr::mutate(ARM = as.character(ARM))

# Combine sample sizes with summary data
baseline_summary <- full_join(baseline_summary, summary_n) %>%
  transmute(ARM, `N/ESS` = N, SEX)

# Insert the ESS as the sample size for the weighted data
# This is calculated above but can also be obtained using the estimate_ess function as shown below
baseline_summary$`N/ESS`[baseline_summary$ARM == "Intervention_weighted"] <- est_weights$analysis_data %>%
  filter(ARM == 'Intervention') %>%
  estimate_ess(wt_col = 'wt')

baseline_summary



## ---- warning=FALSE-------------------------------------------------------------------------------------------------------------------------------------------
# Separate data into intervention data and comparator data
# The boostrap_HR function below makes use of the estimate_weights function
# which requires separate datasets
int <- filter(est_weights$analysis_data, ARM == 'Intervention')
comp <- filter(est_weights$analysis_data, ARM == 'Comparator')

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

#Produce the Kaplan-Meier plot
KM_plot <- ggsurvplot(KM_list,
                      combine = TRUE,
                      risk.table=T, # numbers at risk displayed on the plot
                      break.x.by=50,
                      xlab="Time (days)",
                      censor=FALSE,
                      legend.title = "Treatment",
                      title = "Kaplan-Meier plot of overall survival",
                      legend.labs=c("Intervention", "Intervention weighted", "Comparator"),
                      font.legend = list(size = 10)) +
                      guides(colour=guide_legend(nrow = 2))
KM_plot




## -------------------------------------------------------------------------------------------------------------------------------------------------------------
## Calculate HRs

# Fit a Cox model without weights to estimate the unweighted HR
unweighted_cox <- coxph(Surv(Time, Event==1) ~ ARM, data = est_weights$analysis_data)

HR_CI_cox <- summary(unweighted_cox)$conf.int %>%
              as.data.frame() %>%
              dplyr::transmute(`exp(coef)`,`lower .95`,`upper .95`)
HR_CI_cox

# Fit a Cox model with weights to estimate the weighted HR
weighted_cox <- coxph(Surv(Time, Event==1) ~ ARM, data = est_weights$analysis_data, weights = wt)

HR_CI_cox_wtd <- summary(weighted_cox)$conf.int %>%
                  as.data.frame() %>%
                  dplyr::transmute(`exp(coef)`,`lower .95`,`upper .95`)
HR_CI_cox_wtd

## Bootstrapping

# Bootstrap 1000 HRs
HR_bootstraps <- boot(data = int, # intervention data
                      statistic = bootstrap_HR, # bootstrap the HR (defined in the MAIC package)
                      R=1000, # number of bootstrap samples
                      comparator_data = comp, # comparator pseudo data
                      matching = est_weights$matching_vars, # matching variables
                      model = Surv(Time, Event==1) ~ ARM # model to fit
                      )

# Median of the bootstrap samples
HR_median <- median(HR_bootstraps$t)

# Bootstrap CI - Percentile CI
boot_ci_HR <- boot.ci(boot.out = HR_bootstraps, index=1, type="perc")

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
                             "HR_low_CI" = boot_ci_HR$percent[4],
                             "HR_upp_CI" = boot_ci_HR$percent[5],
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
            dplyr::transmute(Method, HR_95_CI)

# turns the results to a table suitable for word/ powerpoint
HR_table <- HR_summ %>%
    regulartable() %>% #make it a flextable object
  set_header_labels(Method = "Method",  HR_95_CI = "Hazard ratio (95% CI)")  %>%
    font(font = 'Arial', part = 'all') %>%
  fontsize(size = 14, part = 'all') %>%
  bold(part = 'header') %>%
  align(align = 'center', part = 'all') %>%
  align(j = 1, align = 'left', part = 'all') %>%
  border_outer(border = fp_border()) %>%
  border_inner_h(border = fp_border()) %>%
  border_inner_v(border = fp_border()) %>%
  autofit(add_w = 0.2, add_h = 2)
HR_table





## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Summarize bootstrap estimates in a histogram
# Vertical lines indicate the median and upper and lower CIs
hist(HR_bootstraps$t, main = "", xlab = "Boostrapped HR")
abline(v= quantile(HR_bootstraps$t, probs = c(0.025, 0.5, 0.975)), lty=2)




## -------------------------------------------------------------------------------------------------------------------------------------------------------------
## Calculate ORs

# Fit a logistic regression model without weights to estimate the unweighted OR
unweighted_OR <- glm(formula = response~ARM,
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
weighted_OR <- suppressWarnings(glm(formula = response~ARM,
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
                      model = 'response ~ ARM' # model to fit
                      )


# Median of the bootstrap samples
OR_median <- median(OR_bootstraps$t)

# Bootstrap CI - Percentile CI
boot_ci_OR <- boot.ci(boot.out = OR_bootstraps, index=1, type="perc")

# Bootstrap CI - BCa CI
boot_ci_OR_BCA <- boot.ci(boot.out = OR_bootstraps, index=1, type="bca")


## Summary

# Produce summary of ORs and CIs
OR_summ <-  rbind(OR_CI_logit, OR_CI_logit_wtd) %>% # Unweighted and weighted ORs and CIs
            as.data.frame() %>%
            dplyr::rename(OR = `Odds ratio`, OR_low_CI = `2.5 %`, OR_upp_CI = `97.5 %`) %>%
            mutate(Method = c("OR (95% CI) from unadjusted logistic regression model",
                              "OR (95% CI) from weighted logistic regression model")) %>%

            # Median bootstrapped HR and 95% percentile CI
            rbind(data.frame("OR" = OR_median,
                             "OR_low_CI" = boot_ci_OR$percent[4],
                             "OR_upp_CI" = boot_ci_OR$percent[5],
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

            dplyr::transmute(Method, OR_95_CI)

# turns the results to a table suitable for word/ powerpoint
OR_table <- OR_summ %>%
    regulartable() %>% #make it a flextable object
  set_header_labels(Method = "Method",  OR_95_CI = "Odds ratio (95% CI)")  %>%
    font(font = 'Arial', part = 'all') %>%
  fontsize(size = 14, part = 'all') %>%
  bold(part = 'header') %>%
  align(align = 'center', part = 'all') %>%
  align(j = 1, align = 'left', part = 'all') %>%
  border_outer(border = fp_border()) %>%
  border_inner_h(border = fp_border()) %>%
  border_inner_v(border = fp_border()) %>%
  autofit(add_w = 0.2)
OR_table



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Summarize bootstrap estimates in a histogram
# Vertical lines indicate the median and upper and lower CIs
hist(OR_bootstraps$t, main = "", xlab = "Boostrapped OR")
abline(v= quantile(OR_bootstraps$t, probs = c(0.025,0.5,0.975)), lty=2)


