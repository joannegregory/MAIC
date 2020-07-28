
# This example code estimates weights for individual patient data from a single
# arm study of 'intervention' based on aggregate baseline characteristics from
# the comparator trial, and then performs diagnostics on the weights

library(dplyr)
library(MAIC)


#### Prepare the data ----------------------------------------------------------

### Intervention data

# Read in relevant ADaM data and rename variables of interest
adsl <- read.csv(system.file("extdata", "adsl.csv", package = "MAIC",
                             mustWork = TRUE))
adrs <- read.csv(system.file("extdata", "adrs.csv", package = "MAIC",
                             mustWork = TRUE))
adtte <- read.csv(system.file("extdata", "adtte.csv", package = "MAIC",
                              mustWork = TRUE))

adsl <- adsl %>% # Data containing the matching variables
  mutate(SEX=ifelse(SEX=="Male", 1, 0)) # Coded 1 for males and 0 for females

adrs <- adrs %>% # Response data
  filter(PARAM=="Response") %>%
  transmute(USUBJID, ARM, response=AVAL)

adtte <- adtte %>% # Time to event data (overall survival)
  filter(PARAMCD=="OS") %>%
  mutate(Event=1-CNSR) %>% #Set up coding as Event = 1, Censor = 0
  transmute(USUBJID, ARM, Time=AVAL, Event)

# Combine all intervention data
intervention_input <- adsl %>%
  full_join(adrs, by=c("USUBJID", "ARM")) %>%
  full_join(adtte, by=c("USUBJID", "ARM"))

# List out the variables in the intervention data that have been identified as
# prognostic factors or treatment effect modifiers and will be used in the
# matching
match_cov <- c("AGE",
               "SEX",
               "SMOKE",
               "ECOG0")

## Baseline data from the comparator trial
# Baseline aggregate data for the comparator population
target_pop <- read.csv(system.file("extdata", "Aggregate data.csv",
                                     package = "MAIC", mustWork = TRUE))

# Rename target population cols to be consistent with match_cov
target_pop_standard <- target_pop %>%
        rename(
           N=N,
           Treatment=ARM,
           AGE=age.mean,
           SEX=prop.male,
           SMOKE=prop.smoke,
           ECOG0=prop.ecog0
              ) %>%
        transmute(N, Treatment, AGE, SEX, SMOKE, ECOG0)


#### Estimate weights ----------------------------------------------------------

### Center baseline characteristics
# (subtract the aggregate comparator data from the corresponding column of
# intervention PLD)
intervention_data <- intervention_input %>%
    mutate(
     Age_centered = AGE - target_pop$age.mean,
     # matching on both mean and standard deviation for continuous variables (optional)
     Age_squared_centered = (AGE^2) - (target_pop$age.mean^2 + target_pop$age.sd^2),
     Sex_centered = SEX - target_pop$prop.male,
     Smoke_centered = SMOKE - target_pop$prop.smoke,
     ECOG0_centered = ECOG0 - target_pop$prop.ecog0)

## Define the matching covariates
cent_match_cov <- c("Age_centered",
                    "Age_squared_centered",
                    "Sex_centered",
                    "Smoke_centered",
                    "ECOG0_centered")

## Optimization procedure
# Following the centering of the baseline characteristics of the intervention
# study, patient weights can be estimated using estimate_weights
# The function output is a list containing (1) a data set of the individual
# patient data with the assigned weights "analysis_data" and (2) a vector
# containing the matching variables "matching_vars"
est_weights <- estimate_weights(intervention_data = intervention_data,
                                matching_vars = cent_match_cov)

#### Weight diagnostics --------------------------------------------------------

### Are the weights sensible?

# The wt_diagnostics function requires the outputs from the est_weights function
# and will output:
# - the effective sample size (ESS)
# - a summary of the weights and rescaled weights (mean, standard deviation,
#   median, minimum and maximum)
# - a unique set of weights with the corresponding patient profile based on the
#   matching variables

diagnostics <- wt_diagnostics(est_weights$analysis_data,
                              vars = est_weights$matching_vars)

diagnostics$ESS
diagnostics$Summary_of_weights
diagnostics$Weight_profiles

# Each of the wt_diagnostics outputs can also be estimated individually
ESS <- estimate_ess(est_weights$analysis_data)
weight_summ <- summarize_wts(est_weights$analysis_data)
wts_profile <- profile_wts(est_weights$analysis_data, vars = match_cov)

# Plot histograms of unscaled and rescaled weights
# bin_width needs to be adapted depending on the sample size in the data set
histogram <- hist_wts(est_weights$analysis_data, bin = 50)
histogram


### Has the optimization worked?

# The following code produces a summary table of the intervention baseline
# characteristics before and after matching compared with the comparator
# baseline characteristics:

# Create an object to hold the output
baseline_summary <- list('Intervention' = NA,
                         'Intervention_weighted' = NA,
                         'Comparator' = NA)

# Summarise matching variables for weighted intervention data
baseline_summary$Intervention_weighted <- est_weights$analysis_data %>%
  transmute(AGE, SEX, SMOKE, ECOG0, wt) %>%
  summarise_at(match_cov, list(~ weighted.mean(., wt)))

# Summarise matching variables for unweighted intervention data
baseline_summary$Intervention <- est_weights$analysis_data %>%
  transmute(AGE, SEX, SMOKE, ECOG0, wt) %>%
  summarise_at(match_cov, list(~ mean(.)))

# baseline data for the comparator study
baseline_summary$Comparator <- transmute(target_pop_standard,
                                         AGE,
                                         SEX,
                                         SMOKE,
                                         ECOG0)

# Combine the three summaries
# Takes a list of data frames and binds these together
trt <- names(baseline_summary)
baseline_summary <-  bind_rows(baseline_summary) %>%
  transmute_all(sprintf, fmt = "%.2f") %>% #apply rounding for presentation
  transmute(ARM = as.character(trt), AGE, SEX, SMOKE, ECOG0)

# Insert N of intervention  as number of patients
baseline_summary$`N/ESS`[baseline_summary$ARM == "Intervention"] <- nrow(est_weights$analysis_data)

# Insert N for comparator from target_pop_standard
baseline_summary$`N/ESS`[baseline_summary$ARM == "Comparator"] <- target_pop_standard$N

# Insert the ESS as the sample size for the weighted data
# This is calculated above but can also be obtained using the estimate_ess function as shown below
baseline_summary$`N/ESS`[baseline_summary$ARM == "Intervention_weighted"] <- est_weights$analysis_data %>%
  estimate_ess(wt_col = 'wt')

baseline_summary <- baseline_summary %>%
  transmute(ARM, `N/ESS`=round(`N/ESS`,1), AGE, SEX, SMOKE, ECOG0)

