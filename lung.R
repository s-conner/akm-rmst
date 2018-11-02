library(survival)
data(lung)
head(lung)

# Drop individual with missing covariates, recode variables
lung2 <- lung[complete.cases(lung[ , c(2:5,9,10)]),]
lung2$male <- 2-lung2$sex
lung2$status2 <- lung2$status-1

# Obtain weights with logistic model
logit <- glm(male ~ age + meal.cal + wt.loss, data=lung2, family=binomial(link='logit'))
pred <- predict(logit, type='response')
lung2$weight <- lung2$male*pred + (1-lung2$male)*pred

# AKM RMST adjusted for age
akm_rmst(time=lung2$time, status=lung2$status2, group=as.factor(lung2$male), weight=lung2$weight, tau=600)

# Unadjusted
akm_rmst(time=lung2$time, status=lung2$status2, group=as.factor(lung2$male), weight=NULL, tau=600)
