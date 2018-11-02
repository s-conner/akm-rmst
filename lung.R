library(survival)
data(lung)

# Drop individual with missing covariates, recode variables.
# Karnofsky performance scale index dichomotized according to clinical definitions.
lung2 <- lung[complete.cases(lung[ , c(2:9)]),]
lung2$male <- 2-lung2$sex
lung2$status2 <- lung2$status-1
lung2$ph.karno.low <- 0
lung2$ph.karno.low[lung2$ph.karno<=70] <- 1

# Obtain weights with logistic model
# dichotomize by <=70, adjust for age sex PH.ecog and meal.cal
logit <- glm(ph.karno.low ~ male + age + meal.cal + ph.ecog, data=lung2, family=binomial(link='logit'))
pred <- predict(logit, type='response')
lung2$weight <- lung2$ph.karno.low*pred + (1-lung2$ph.karno.low)*(1-pred)

# AKM RMST adjusted for age
akm_rmst(time=lung2$time, status=lung2$status2, group=as.factor(lung2$ph.karno.low), weight=lung2$weight, tau=600)

# Unadjusted (null weights)
akm_rmst(time=lung2$time, status=lung2$status2, group=as.factor(lung2$ph.karno.low), weight=NULL, tau=600)
