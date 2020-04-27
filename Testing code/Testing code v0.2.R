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
base_dir <- 'G:/Clients/Roche/2797 Development of R Code for MAIC and mixture cure models/Project/2 Exploratory'
data_path <- file.path(base_dir,'Simulated datasets')



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

cent_match_cov <- names(intervention_data)[11:15]
cent_match_cov

match_cov <- names(intervention_data)[4:7]
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
  #       - A charcacter vector with the name sof the matching variables
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


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Weight diagnostics ----------------------------------------------

# Weight diagnostics with the intervention data only:
data_for_diagnostics <- est_weights$analysis_data %>%
                        filter(ARM=="A")

#### ESS

# Function to caclulate ESS:
estimate_ess <- function(data, wt=wt){
  ess <- sum(data$wt)^2/sum(data$wt^2)
  return(ess)
}

ESS <- estimate_ess(data_for_diagnostics)
ESS


#### Weights summary

# Function to caclulate min, max, median, mean and sd of wts:
summarize_wts <- function(data, wt=wt){
  summary <- data %>%
    summarise(
      min = min(data$wt),
      max = max(data$wt),
      median = median(data$wt),
      mean = mean(data$wt),
      sd = sd(data$wt)
    )
  return(summary)
}

weight_summ <- summarize_wts(data_for_diagnostics)
weight_summ


##### Weight profiles

profile_wts <- function(data, wt=wt, wt_rs=wt_rs, vars){
  profile_data <-  data %>%
                  select_if(names(data) %in% vars ==TRUE)

  wts <- data %>%
          select(wt, wt_rs)

  profile_wts <- cbind(wts, profile_data) %>%
                 distinct()

  return(profile_wts)
}

profile_wts(data=data_for_diagnostics, vars = est_weights$matching_vars)


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

histogram <- hist_wts(data_for_diagnostics)
histogram

histogram_bw1 <- hist_wts(data_for_diagnostics, bin_width=1)
histogram_bw1



#### All weight diagnostics

##### FUNCTION NOT WORKING
# Function to combine weight diagnostic functions above:
all_wt_diagnostics <- function(data, # analysis data from estimate_weights
                               #arm,
                               vars,
                               wt=wt){

  #intervention_data <- data %>% filter(ARM==arm)

  #ESS
  ESS <- estimate_ess(data)

  #Summary
  summ_wts <- summarize_wts(data)

  #Weight profiles
  #profile <- profile_wts(data, vars)

  output <- list("ESS" = ESS,
                 "Summary_of_weights" = summ_wts#,
                 #"Histogram_of_weights" = hist_wts(intervention_data),
                 # "Weight_profiles" = profile
                 )
  return(output)
}


all_wt_diagnostics(data_for_diagnostics, vars = est_weights$matching_vars)

# all_wt_diagnostics(data = est.weights$analysis_data,
#                    arm = "A",
#                    vars = est.weights$matching_vars)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Bootstrapping ------------------------------------------------

boostrap_OR <- function(intervention_data, i,  cent_vars, comparator_data){

  # Samples the data
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

  #Combines intervention dataset with comparator.data
  analysis_data <- bootstrap_dat_wt %>%
                   rbind.fill(comparator_data)

  # binary data
  prop_A<-sum(bootstrap.dat.wt$weights*bootstrap.dat.wt[,resp])/sum(bootstrap.dat.wt$weights)
  prop_B<-sum(comparator.data[,resp])/nrow(comparator.data)
  logistic.regr <- glm(Binary_event~trt, family=binomial(link="logit"), data = analysis.data, weights = weights)
  OR <- exp(as.numeric(coef(logistic.regr)["trtB"]))
  log_OR <- as.numeric(coef(logistic.regr)["trtB"])

  c("OR" = OR, "Log_OR" = log_OR, "prop_A" = prop_A, "prop_B" = prop_B)
}

# Bootstrap function

boostrap.func <- function(intervention_data, i, vars, comparator.data, resp){
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

  bootstrap.dat.wt <- cbind(bootstrap.data, weights)

  #Combines intervention dataset with comparator.data
  analysis.data <- bootstrap.dat.wt %>%
    #select(trt, Binary_event, weights) %>%
    rbind.fill(comparator.data)

  # binary data
  prop_A<-sum(bootstrap.dat.wt$weights*bootstrap.dat.wt[,resp])/sum(bootstrap.dat.wt$weights)
  prop_B<-sum(comparator.data[,resp])/nrow(comparator.data)
  logistic.regr <- glm(Binary_event~trt, family=binomial(link="logit"), data = analysis.data, weights = weights)
  OR <- exp(as.numeric(coef(logistic.regr)["trtB"]))

  # survival data
  fit.A <- surv_fit(Surv(Time, Event) ~ 1, data = bootstrap.dat.wt, weights=weights)
  median_A<-surv_median(fit.A)$median
  #
  fit.B <- surv_fit(Surv(Time, Event) ~ 1, data = comparator.data)
  median_B<-surv_median(fit.B)$median

  cox_model <- coxph(Surv(Time, Event==1) ~ trt, data = analysis.data, weights = weights)
  HR <- exp(cox_model$coefficients)

  c("OR"=OR, "prop_A"=prop_A, "prop_B"=prop_B,"HR"=HR, "median_A"=median_A, "median_B"=median_B)
}





bootstraps <- boot(intervention, boostrap.func, R=1000, vars = c("Smoke.centered", "ECOG0.centered"), comparator.data=comparator.data, resp="Binary_event")
boot_data <- as.data.frame(bootstraps$t)
colnames(boot_data) <- colnames(t(as.data.frame(bootstraps$t0)))
head(boot_data)
names(boot_data)
bootstraps$t0
summary(bootstraps$t)


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


###############################
# JG attempt - this bootstraps both comparator and int data

JG_boot<-function(data,indices,cova_vec,AD_frame,resp){

# data=usedata2
# cova_vec=c("Smoke.centered","ECOG0.centered")
# AD_frame=target_pop2
# resp="Binary_event"
# wt="wt"
# indices=c(1:nrow(data))


dat<-data[indices,]
#dat<-data

dt <- subset(dat, ARM=="A")

cova<-cova_vec
# Objective function
objfn <- function(a1, X){
  sum(exp(X %*% a1))
}
# Gradient function
gradfn <- function(a1, X){
  colSums(sweep(X, 1, exp(X %*% a1), "*"))
}
# Centred EMs

X.EM.0 <- sweep(as.matrix(dt[,cova]), 2,
                as.numeric(AD_frame[1,cova]), '-')
# Estimate weights
opt1 <- optim(par = rep(0,length(cova)), fn = objfn, gr = gradfn, X = X.EM.0, method = "BFGS")

a1 <- opt1$par
wt <- exp(X.EM.0 %*% a1)
wts <- (wt / sum(wt)) *(sum(wt)^2/sum(wt^2))   # rescaled weights
ess<-(sum(wt)^2/sum(wt^2))
dt<-cbind(dt,wts)
dt_comp <- subset(dat, ARM=="B") %>%
  mutate(wts=1)

data_all <- rbind.fill(dt,dt_comp)

cr<-sum(dt$wts*dt[,resp])/sum(dt$wts)
cr_comp<-sum(dt_comp$wts*dt_comp[,resp])/nrow(dt_comp)
logistic.regr <- glm(Binary_event~ARM, family=binomial(link="logit"), data = data_all)
log_OR <- as.numeric(coef(logistic.regr)["trtB"])

fit.A <- surv_fit(Surv(Time, Event) ~ 1, data = dt,weights=wts)
median_A<-surv_median(fit.A)$median
#
fit.B <- surv_fit(Surv(Time, Event) ~ 1, data = dt_comp)
median_B<-surv_median(fit.B)$median

cox_model <- coxph(Surv(Time, Event==1) ~ trt, data = data_all, weights = wts)
HR <- exp(cox_model$coefficients)

c("prop_A"=cr, "prop_b"=cr_comp, "logOR"=log_OR, "median_A"=median_A, "median_B"=median_B,HR=HR)

}



usedata2 <-rbind.fill(intervention,  comparator_ipd)
target_pop2 <- target_pop %>% dplyr::rename(Smoke.centered=prop.smoke, ECOG0.centered=prop.ecog0)

myBootstrap <- boot(
  usedata2, JG_boot, R=10,strata=usedata2$trt,
  cova_vec=c("Smoke.centered","ECOG0.centered"), AD_frame=target_pop2, resp="Binary_event")



myBootstrap$t
summary(myBootstrap$t)
myBootstrap$t0

plot(myBootstrap, index=1)
plot(myBootstrap, index=4)
boot.ci(boot.out = myBootstrap,index=1)
boot.ci(boot.out = myBootstrap,index=2)
boot.ci(boot.out = myBootstrap,index=3)
boot.ci(boot.out = myBootstrap,index=4)


#

# summaries ---------------------------------------------------------------

##############################

Baseline_summary_bin <- MAIC_analysis_dataset %>%
  filter(Treatment!=Treatment_comparator) %>%
  select(Treatment, matchingvars_binary, Weights) %>%
  dplyr::group_by(Treatment) %>%
  summarise_each(funs(weighted.mean(.*100, Weights)), -Weights)  %>%
  rbind.fill(target_pop_standard  %>% select(matchingvars_binary) %>% transmute_all(funs(.*100)) %>% mutate(Treatment=Treatment_comparator))

Baseline_summary_cont <- MAIC_analysis_dataset %>%
  filter(Treatment!=Treatment_comparator) %>%
  select(Treatment, matchingvars_cont, Weights) %>%
  dplyr::group_by(Treatment) %>%
  summarise_each(funs(weighted.mean(., Weights)),-Weights) %>%
  rbind.fill(target_pop_standard  %>% select(matchingvars_cont) %>% mutate(Treatment=Treatment_comparator))

Baseline_summary_n <- MAIC_analysis_dataset %>%
  filter(Treatment!=Treatment_comparator) %>%
  #select(Treatment, matchingvars) %>%
  dplyr::group_by(Treatment) %>%
  dplyr::summarise(
    'N' = n()) %>%
  rbind.fill(target_pop_standard  %>% select(N) %>% mutate(Treatment=Treatment_comparator)) %>%
  rename(`N/ESS`=N)



Baseline_summary_all <- plyr::join_all(list(Baseline_summary_n, Baseline_summary_bin, Baseline_summary_cont),by="Treatment")

Baseline_summary_all2 <- cbind(Baseline_summary_all %>% select(-matchingvars),
                               lapply(Baseline_summary_all %>% select(matchingvars), sprintf, fmt = "%.1f") %>% as.data.frame())

Baseline_summary_all2$Study <- if_else(Baseline_summary_all2$Treatment == intervention_matched, intervention_study,
                                       if_else(Baseline_summary_all2$Treatment == intervention_unadjusted, intervention_study, comparator_study))

# replace N with ESS
Baseline_summary_all2$`N/ESS`[Baseline_summary_all2$Treatment == intervention_matched] <- ESS

###################################
## Shouldn't need any changes from this point
######################################

# KMs ----------------------------------------------------------------------------------------
message('KMs')

# Km estimate
km.est <- survfit(Surv(Time, Event) ~ Treatment, data = MAIC_analysis_dataset, conf.type = 'plain', weights = Weights)
names(km.est$strata) <- gsub('Treatment=', '', names(km.est$strata))

# Plot KM

km.plot <- ggsurvplot(km.est, data=MAIC_analysis_dataset, risk.table = TRUE,
                      break.time.by = time.breaks,
                      conf.int = FALSE,
                      censor=FALSE,
                      legend.title = '',
                      xlab = paste0('Time (', time.unit, ')'),
                      palette = c(BM.blue, BM.red, BM.Dyellow),
                      font.x = 16,
                      font.y = 18,
                      font.legend = 16,
                      font.xtickslab = 16,
                      font.ytickslab = 16,
                      fontsize = 6,
                      xlim = c(0, max(MAIC_analysis_dataset$Time) + 1))
km.plot$table <- ggpar(
  km.plot$table,
  font.x        = c(16),
  font.xtickslab = c(16),
  font.ytickslab = c(16)
)

jpeg(file = file.path(output_file, "KM.jpeg"), width = 25, height = 20, units = 'cm', res = 300, quality = 100)
km.plot
graphics.off()


km.summ <- summary(km.est)$table %>%
  as.data.frame()  %>%
  rownames_to_column(var = 'Treatment') %>%
  select(Treatment, "N/ ESS" = records, Events = events, Median = median, LowerCI = `0.95LCL`, UpperCI = `0.95UCL`)

km.summ$Study <- if_else(km.summ$Treatment == intervention_matched, intervention_study,
                         if_else(km.summ$Treatment == intervention_unadjusted, intervention_study, comparator_study))

km.summ$`N/ ESS`[km.summ$Treatment == intervention_matched] <- ESS


# Hazard ratios ----------------------------------------------------------------------------------------
message('Hazard ratios unweighted')
unweighted_data.ref1 <- MAIC_analysis_dataset %>%
  filter(Treatment!=intervention_matched) %>%
  mutate(Treatment = factor(Treatment, levels = c(Treatment_comparator, intervention_unadjusted))) # change ref

unweighted_data.ref2 <- MAIC_analysis_dataset %>%
  filter(Treatment!=intervention_matched) %>%
  mutate(Treatment = factor(Treatment, levels = c(intervention_unadjusted,Treatment_comparator))) # change ref

unweighted.cox.ref1 <- coxph(Surv(Time, Event==1) ~ Treatment, data = unweighted_data.ref1)
unweighted.cox.ref2 <- coxph(Surv(Time, Event==1) ~ Treatment, data = unweighted_data.ref2)

#extract HRs
cox.summ.unweighted <- rbind(summary(unweighted.cox.ref1)$conf.int, summary(unweighted.cox.ref2)$conf.int) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'Treatment') %>%
  mutate(Treatment = sub('Treatment', '', Treatment),
         Comparison = c(paste0(Treatment_intervention, " vs ", Treatment_comparator),
                        paste0(Treatment_comparator, " vs ", Treatment_intervention)),
         Method = "Unadjusted") %>%
  select(-`exp(-coef)`) %>% #drop unnecessary column
  rename(HR = `exp(coef)`, HR.low.CI = `lower .95`, HR.upp.CI = `upper .95`)


# Weighted survival analysis -----------------------------------------------------------------------------------------------

message('Weighted survival analysis')
weighted_data.ref1 <- MAIC_analysis_dataset %>%
  filter(Treatment!=intervention_unadjusted) %>%
  mutate(Treatment = factor(Treatment, levels = c(Treatment_comparator, intervention_matched))) # change ref

weighted_data.ref2 <- MAIC_analysis_dataset %>%
  filter(Treatment!=intervention_unadjusted) %>%
  mutate(Treatment = factor(Treatment, levels = c(intervention_matched,Treatment_comparator))) # change ref


weighted.cox.ref1 <- coxph(Surv(Time, Event==1) ~ Treatment, data = weighted_data.ref1, weights = Weights)
weighted.cox.ref2 <- coxph(Surv(Time, Event==1) ~ Treatment, data = weighted_data.ref2, weights = Weights)

#extract HRs
cox.summ.weighted <- rbind(summary(weighted.cox.ref1)$conf.int, summary(weighted.cox.ref2)$conf.int) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'Treatment') %>%
  mutate(Treatment = sub('Treatment', '', Treatment),
         Comparison = c(paste0(Treatment_intervention, " vs ", Treatment_comparator),
                        paste0(Treatment_comparator, " vs ", Treatment_intervention)),
         Method = "Weighted standard CI") %>%
  select(-`exp(-coef)`) %>% #drop unnecessary column
  rename(HR = `exp(coef)`, HR.low.CI = `lower .95`, HR.upp.CI = `upper .95`)

# Bootstrap HRs ------------------------------------------------------------------------------------------------
HR.ref1.Bootstrap <- rep(NA,n.sim)
HR.ref2.Bootstrap <- rep(NA,n.sim)
log.HR.ref1.Bootstrap <- rep(NA,n.sim)
HR.Bootstrap.scaled.weights <- rep(NA,n.sim)
