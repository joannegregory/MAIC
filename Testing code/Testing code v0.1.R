# # # # # #  # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 2797 exploratory -testing out the code                              #
# Author: SS/ JG (23.04.2020)                                   #
# # # # # #  # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Load libraries --------------------------------------------------
library(plyr)
library(dplyr)
library(haven)
library(ggplot2)
library(tidyr)
library(boot)
library(survminer)
library(survminer)


# Load data ---------------------------------------------------------------



# Intervention data
base.dir <- 'G:/Clients/Roche/2797 Development of R Code for MAIC and mixture cure models/Project/2 Exploratory'
data.path <- file.path(base.dir,'Simulated datasets')

# Read in ADaM data and rename
adsl <- read.csv(file.path(data.path, "adsl.csv")) %>%
  mutate(SEX=ifelse(SEX=="Male",1,0))

adrs <- read.csv(file.path(data.path, "adrs.csv")) %>%
  filter(PARAM=="Response") %>%
  select(USUBJID, ARM, Binary_event=AVAL)

adtte <- read.csv(file.path(data.path, "adtte.csv")) %>%
  filter(PARAMCD=="OS") %>%
  mutate(Event=1-CNSR) %>%
  select(USUBJID, ARM, Time=AVAL, Event)


intervention_input <- join_all(list(adsl, adrs, adtte), type = "full", by=c("USUBJID", "ARM"))



comparator.surv <- read.csv(file.path(data.path,"psuedo_IPD.csv"))
comparator.n <- 300
comparator.prop.events <- 0.4
comparator.binary <- data.frame("Binary_event"=c(rep(1,comparator.n*comparator.prop.events),rep(0, comparator.n*(1-comparator.prop.events))))

# join comparator data note not a 1:1 relationship - the rows do not represent the same observation
comparator_input <- cbind(comparator.surv, comparator.binary)

# Baseline agregate data
target_pop <- read.csv(file.path(data.path,"Aggregate data.csv"))
target_pop

##### Calculate weights ------------------------------------------------------------

# center baseline characteristics
intervention_data <- intervention_input %>%
                          mutate(Age.centered = AGE - target_pop$age.mean,
                                 Sex.centered = SEX - target_pop$prop.male,
                                Smoke.centered = SMOKE - target_pop$prop.smoke,
                                 ECOG0.centered = ECOG0 - target_pop$prop.ecog0)
head(intervention)
match.cov <- c("Age.centered", "Sex.centered", "Smoke.centered", "ECOG0.centered")

#### optimisation procedure and calculation of weights

# This is a fuction for Q(b)
objfn <- function(a1, X){
  sum(exp(X %*% a1))
}


# Gradient function => Derivative of Q(b).
gradfn <- function(a1, X){
  colSums(sweep(X, 1, exp(X %*% a1), "*"))
}


## The following is the in built R function used to optimise Q(b) ##
## using Newton-Raphson techniques ##

estimate_weights <- function(intervention_data, vars, comparator_data){
  intervention_data=intervention_data
  vars=match.cov
  comparator_data=comparator_input

  print(opt1 <- optim(par = rep(0,dim(intervention_data[, vars])[2]),
                      fn = objfn,
                      gr = gradfn,
                      X = as.matrix(intervention_data[, vars]),
                      method = "BFGS"))

  a1 <- opt1$par


  # Calculation of weights.
  wt <- as.vector(exp(as.matrix(intervention_data[, vars]) %*% a1))

  # rescaled weights
  wt.rs <- (wt / sum(wt)) * dim(intervention_data)[1]

  # combine data with weights
  data.with.wts <- cbind(intervention_data, wt, wt.rs)
# add weight=1 to comparator data
  comparator_data.wts <- comparator_data %>% mutate(wt=1, wt.rs=1)
 # Add on comparator data
   all_data <- rbind.fill(data.with.wts, comparator_data.wts)
  output <- list(analysis_data=all_data, matching_vars=vars)

  return(output)

}

est.weights <- estimate_weights(intervention_data=intervention_data, vars=match.cov, comparator_data=comparator_input)
head(est.weights$analysis_data)

data_for_diagnostics <- est.weights$analysis_data %>% filter(ARM=="A")

#### Weight diagnostics ----------------------------------------------------------------------

# ESS
estimate.ess <- function(data, wt=wt){
  ess <- sum(data$wt)^2/sum(data$wt^2)
  return(ess)
}
ESS <- estimate.ess(data_for_diagnostics)
ESS


# Weights summary
summarize.wts <- function(intervention_wts, wt=wt){
  summary <- intervention_wts %>%
    summarise(
      min = min(intervention_wts$wt),
      max = max(intervention_wts$wt),
      median = median(intervention_wts$wt),
      mean = mean(intervention_wts$wt),
      sd = sd(intervention_wts$wt)
    )
  return(summary)
}

weight.summ <- summarize.wts(data_for_diagnostics)
weight.summ


# Histograms


hist_wts <- function(data, wt_col="wt", rs_wt_col="wt.rs") {

wt.data <- data[,c(wt_col, rs_wt_col)] %>%
            rename("Weights" = wt_col, "Rescaled weights" = rs_wt_col) %>%
            gather()

hist.plot <- qplot(data = wt.data,
                    value,
                    geom="histogram",
                    xlab = "Histograms of weights and rescaled weights",
                    binwidth=0.05) +
              facet_wrap(~key,  ncol=1) +
              theme_bw()+
              theme(axis.title = element_text(size = 16),
                    axis.text = element_text(size = 16))+
              ylab("Frequency")

return(hist.plot)
}

histogram <- hist_wts(data_for_diagnostics)
histogram


# Weight profiles

profile_wts <- function(data, wt_col="wt", vars){
             profile_data <-  data[, c(wt_col, vars)] %>%
             distinct()
             return(profile_data)
              }

profile_wts(data_for_diagnostics, vars = est.weights$matching_vars)


# all weights - psuedo code
all_diagnostics <- function(est.weights, arm="A"){ # object from estimate_weights
  x <- est.weights$analysis_data %>% filter(ARM==arm)
  vars <- est.weights$matching_vars
  output <- list("ESS"=estimate.ess(x),
                 "summarize.wts"=summarize.wts(x)
                # "hist_wts"=hist_wts(data),"profile_wts"=profile_wts(data, vars)
                 )
  return(output)
}

all_diagnostics(est.weights, arm="A")
#### Bootstrapping --------------------------------------------------------------------------

# Bootstrap function

boostrap.func <- function(data, i, vars, comparator.data, resp){
  # Samples the data
  bootstrap.data <- data[i,]

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



# End of bootstrapping ---------------------------------------------------------------------------------------------------------------





# CI for weighted fedratinib
hist(fed.TSS.Bootstrap,main = "",xlab = "Boostrapped Percentage")
abline(v= quantile( fed.TSS.Bootstrap,probs = c(0.025,0.5,0.975)), lty=2)
quantile(fed.TSS.Bootstrap, probs = c(0.025,0.5,0.975))

reweighted.TSS.percentage.LCI <- data.frame(LCI = quantile(fed.TSS.Bootstrap, probs = c(0.025)), row.names=NULL)
reweighted.TSS.percentage.UCI <- data.frame(UCI = quantile(fed.TSS.Bootstrap, probs = c(0.975)), row.names=NULL)


# CI for RD
hist(RD.TSS.Bootstrap, main = "",xlab = "Boostrapped Percentage")
abline(v= quantile(RD.TSS.Bootstrap,probs = c(0.025,0.5,0.975)), lty=2)
quantile(RD.TSS.Bootstrap, probs = c(0.025,0.5,0.975))

LCI.TSS.BS <- data.frame(LCI = quantile(RD.TSS.Bootstrap, probs = c(0.025)), row.names=NULL)
UCI.TSS.BS <- data.frame(UCI = quantile(RD.TSS.Bootstrap, probs = c(0.975)), row.names=NULL)

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

dt <- subset(dat, trt=="A")

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
dt_comp <- subset(dat, trt=="B") %>%
  mutate(wts=1)

data_all <- rbind.fill(dt,dt_comp)

cr<-sum(dt$wts*dt[,resp])/sum(dt$wts)
cr_comp<-sum(dt_comp$wts*dt_comp[,resp])/nrow(dt_comp)
logistic.regr <- glm(Binary_event~trt, family=binomial(link="logit"), data = data_all)
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
