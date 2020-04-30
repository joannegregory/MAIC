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
