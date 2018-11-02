# IPW-adjusted Kaplan-Meier restricted mean survival times
This project contains an R function for inverse probability weighted (IPW) adjusted restricted mean survival times in observational studies. We also provide a working example using the dataset 'lung' from the {survival} package.

We acknowledge the source code of F. Le Borgne and Y. Foucher, authors of the 'adjusted.KM' function from {IPWsurvival} package, and Hajime Uno, Lu Tian, Angel Cronin, Chakib Battioui, and Miki Horiguchi, authors of the 'rmst1 function' from {survRM2} package.

akm_RMST.R
```
This function will produce IPW-adjusted RMSTs, differences and ratios in RMSTs for pair-wise comparisons, and adjusted Kaplan-Meier survival curves. The function requires the following arguments:

Time: time to event
Status: 0 if censored, 1 if event
Group: factor variable for the exposure of interest
Weights: to be obtained previously, ie through logistic models
Tau: a user-specified truncation point. 
     If not specified, the default will be the minimum of the each groups' last event time 
```

lung.R
```
We show how to apply our function to the lung dataset, available in the {survival} package. 
We examine the mean survival time (days) between high and low Karnofsky performance score 
rated by a physician, adjusted for sex, age, calories consumed, and ECOG performance score. 
We obtain inverse probability weights using a logistic model, and apply these weights in our 
function to obtain IPW-adjusted Kaplan-Meier curves and RMSTs.
```
