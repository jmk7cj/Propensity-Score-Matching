#----------------------------------------------------------------------------------------------#
# Author: Joseph Kush (jkush1@jhu.edu) (jmk7cj@virginia.edu)
# 
# Title: Selecting the optimal number of matched cases using propensity score matching methods 
#        for quasi-experiments testing intervention scale-ups R code for Monte Carlo 
#        simulations
#
# Date: 4/20/2021
#
# Purpose: Master .R file to set up and run a Monte Carlo simulation study for propensity 
#          score matching with replacement for scenarios with >50% treatment
#          Step 1: Generate facets for Monte Carlo study
#          Step 2: Generate potential outcomes data, estimate PSM, calculate bias, MSE, etc.
#          Step 3: Store and view results
#----------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------#
# Step 1: Generate facets for Monte Carlo study
#----------------------------------------------------------------------------------------------#
# Load necessary packages, remove all objects
library("foreach"); library("doParallel"); library("doSNOW"); library("MatchIt"); 
library("survey"); library("cobalt"); library("Matching"); library("sjstats"); 
rm(list = ls())

# Facets to vary: 7 x 4 x 3 x 2 = 168 conditions
set.seed(123)
sims <- 1:1 # 1 replication as example
knn <- 1:7 # number of nearest neighbors
sample_size <- c(100, 500, 1000, 5000) # sample size
covs <- c(4, 10, 20) # number of covariates
prop_treat <- c(1, 2) # two different proportion of treatment

# Create progress bar for simulation based on all conditions + replications
total_sim_conditions <- length(unique(sims)) * 
length(unique(knn)) * length(unique(sample_size)) *
length(unique(covs)) * length(unique(prop_treat))

processors <- makeCluster(detectCores()[1]-1) 
registerDoSNOW(processors)
pb <- txtProgressBar(max=total_sim_conditions, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
#----------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------#
# Step 2: Generate potential outcomes data, estimate PSM, calculate bias, MSE, etc.
#----------------------------------------------------------------------------------------------#
sim_function =
foreach(i = sims, .combine=rbind) %:% 
foreach(r = knn, .combine=rbind) %:%
foreach(ss = sample_size, .combine=rbind) %:%
foreach(xs = covs, .combine=rbind) %:%
foreach(prop = prop_treat, .combine=rbind, .options.snow=opts) %dopar% {

library("MatchIt"); library("survey"); library("cobalt"); library("Matching")

#####################
# 4 covariates
#####################
if (xs==4) {
  
generateData <- function(n, prob=1) { 

# Define logistic regression model to relate Pr(treatment==1) to c covariates
# medium effect size for covariates 
alpha <- log(1.25)

# standardized difference d between control and treatment
beta <- .2 

# treatment effect
TE <- 1 

# 50% treatment
if (prob == 1) {
  intercept <- 0
} 

# 80% treatment
else if (prob == 2) {
  intercept <- 1.5
} 

# Create variables
for (i in 1:4) {
assign(paste("x",i, sep = ""), rnorm(n,0,1))
}

# Probability of treatment (equation 6)
logit.treat <- intercept + alpha*x1 + alpha*x2 + alpha*x3 + alpha*x4 

# Exponentiate
p.treat <- exp(logit.treat)/(1 + exp(logit.treat))

# Binary treatment indicator with probability of treatment as sampling probability
treat <- rbinom(n, 1, p.treat)

# Potential outcome for control
y0 <- 0 + beta*x1 + beta*x2 + beta*x3 + beta*x4 +  rnorm(n, 0, 1)

# Potential outcome for treatment 
y1 <- 0 + TE*treat + beta*x1 + beta*x2 + beta*x3 + beta*x4 + rnorm(n, 0, 1)

# Observed outcome
y_obs <- ifelse(treat==1, y1, y0)

return(data.frame(treat, y0,y1,y_obs,x1,x2,x3,x4))
}

# Return dataframe
data <- generateData(n=ss, prob=prop)

# Define population ATT
pop_att <- by(data$y1, data$treat, mean)[[2]] - by(data$y0, data$treat, mean)[[2]]

# Estimate propensity scores and match accordingly 
match = matchit(treat ~ x1+x2+x3+x4, 
data=data, distance="logit", method="nearest", ratio=r, replace=T, verbose=T) 

# Retain matched units only  
matched_data = match.data(match)

# Use weights as matching with replacement
est_att <- svydesign(ids=~1, weights=~weights, data=matched_data)

# Estimate ATT
est_ATT1 <- svyglm(y_obs~treat, est_att, family=gaussian())

# Assess balance
cov_diff <- ((bal.tab(match, s.d.denom="pooled", m.threshold=0.1)$Balanced.mean.diffs[[1]][1])-1)/4

# Store estimated ATT
treat_est <- est_ATT1$coefficients[2]

# Calculate bias, variance, and MSE
bias <- round(((est_ATT1$coefficients[2]-pop_att)/pop_att), digits=3)
var <- (summary(est_ATT1)[13]$coefficients[2,2])^2
mse <- bias^2 + var
}

  
#####################  
# 10 covariates
#####################
else if (xs==10) {
  
generateData <- function(n, prob=1) { 
  
# Define logistic regression model to relate Pr(treatment==1) to c covariates
# medium effect size for covariates 
alpha <- log(1.25)

# standardized difference d between control and treatment
beta <- .2 

# treatment effect
TE <- 1 

# 50% treatment
if (prob == 1) {
  intercept <- 0
} 

# 80% treatment
else if (prob == 2) {
  intercept <- 1.5
} 

# Create variables
for (i in 1:10) {
assign(paste("x",i, sep = ""), rnorm(n,0,1))
}

# Probability of treatment (equation 6)
logit.treat <- intercept + alpha*x1 + alpha*x2 + alpha*x3 + alpha*x4 +
  alpha*x5 + alpha*x6 + alpha*x7 + alpha*x8 + alpha*x9 + alpha*x10

# Exponentiate 
p.treat <- exp(logit.treat)/(1 + exp(logit.treat))

# Binary treatment indicator with probability of treatment as sampling probability
treat <- rbinom(n, 1, p.treat)

# Potential outcome for control
y0 <- 0 + beta*x1 + beta*x2 + beta*x3 + beta*x4 + 
  beta*x5 + beta*x6 + beta*x7 + beta*x8 + beta*x9 + beta*x10 + rnorm(n, 0, 1)

# Potential outcome for treatment
y1 <- 0 + TE*treat + beta*x1 + beta*x2 + beta*x3 + beta*x4 + 
  beta*x5 + beta*x6 + beta*x7 + beta*x8 + beta*x9 + beta*x10 + rnorm(n, 0, 1)

# Observed outcome
y_obs <- ifelse(treat==1, y1, y0)

return(data.frame(treat, y0,y1,y_obs,x1,x2,x3,x4,x5,
                  x6,x7,x8,x9,x10))
}

# Return dataframe
data <- generateData(n=ss, prob=prop)

# Define population ATT
pop_att <- by(data$y1, data$treat, mean)[[2]] - by(data$y0, data$treat, mean)[[2]]

# Estimate propensity scores and match accordingly 
match = matchit(treat ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, 
data=data, distance="logit", method="nearest", ratio=r, replace=T, verbose=T) 

# Retain matched units only
matched_data = match.data(match)

# Use weights as matching with replacement
est_att <- svydesign(ids=~1, weights=~weights, data=matched_data)

# Estimate ATT
est_ATT1 <- svyglm(y_obs~treat, est_att, family=gaussian())

# Assess balance
cov_diff <- ((bal.tab(match, s.d.denom="pooled", m.threshold=0.1)$Balanced.mean.diffs[[1]][1])-1)/10

# Store estimated ATT
treat_est <- est_ATT1$coefficients[2]

# Calculate bias, variance, and MSE
bias <- round(((est_ATT1$coefficients[2]-pop_att)/pop_att), digits=3)
var <- (summary(est_ATT1)[13]$coefficients[2,2])^2
mse <- bias^2 + var
}


#####################
# 20 covariates
#####################
else if (covs==20) {
  
generateData <- function(n, prob=1) { 
  
# Define logistic regression model to relate Pr(treatment==1) to c covariates
# medium effect size for covariates 
alpha <- log(1.25)

# standardized difference d between control and treatment
beta <- .2 

# treatment effect
TE <- 1 

# 50% treatment
if (prob == 1) {
  intercept <- 0
} 

# 80% treatment
else if (prob == 2) {
  intercept <- 1.5
} 

# Create variables
for (i in 1:20) {
assign(paste("x",i, sep = ""), rnorm(n,0,1))
}

# Probability of treatment (equation 6)
logit.treat <- intercept + alpha*x1 + alpha*x2 + alpha*x3 + alpha*x4 +
  alpha*x5 + alpha*x6 + alpha*x7 + alpha*x8 + alpha*x9 + alpha*x10 + 
  alpha*x11 + alpha*x12 + alpha*x13 + alpha*x14 + alpha*x15 + 
  alpha*x16 + alpha*x17 + alpha*x18 + alpha*x19 + alpha*x20
  
# Exponentiate 
p.treat <- exp(logit.treat)/(1 + exp(logit.treat))

# Binary treatment indicator with probability of treatment as sampling probability
treat <- rbinom(n, 1, p.treat)

# Potential outcome for control
y0 <- 0 + beta*x1 + beta*x2 + beta*x3 + beta*x4 + 
  beta*x5 + beta*x6 + beta*x7 + beta*x8 + beta*x9 + beta*x10 +
  beta*x11 + beta*x12 + beta*x13 + beta*x14 + beta*x15 +
  beta*x16 + beta*x17 + beta*x18 + beta*x19 + beta*x20 +
  rnorm(n, 0, 1)

# Potential outcome for treatment
y1 <- 0 + TE*treat + beta*x1 + beta*x2 + beta*x3 + beta*x4 + 
  beta*x5 + beta*x6 + beta*x7 + beta*x8 + beta*x9 + beta*x10 +
    beta*x11 + beta*x12 + beta*x13 + beta*x14 + beta*x15 +
  beta*x16 + beta*x17 + beta*x18 + beta*x19 + beta*x20 +
  rnorm(n, 0, 1)

# Observed outcome
y_obs <- ifelse(treat==1, y1, y0)

return(data.frame(treat, y0,y1,y_obs,x1,x2,x3,x4,x5,
                  x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,
                  x16,x17,x18,x19,x20))
}

# Return dataframe
data <- generateData(n=ss, prob=prop)

# Define population ATT
pop_att <- by(data$y1, data$treat, mean)[[2]] - by(data$y0, data$treat, mean)[[2]]

# Estimate propensity scores and match accordingly 
match = matchit(treat ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+
                  x11+x12+x13+x14+x15+x16+x17+x18+x19+x20, 
data=data, distance="logit", method="nearest", ratio=r, replace=T, verbose=T) 

# Retain matched units only
matched_data = match.data(match)

# Use weights as matching with replacement
est_att <- svydesign(ids=~1, weights=~weights, data=matched_data)

# Estimate ATT
est_ATT1 <- svyglm(y_obs~treat, est_att, family=gaussian())

# Assess balance
cov_diff <- ((bal.tab(match, s.d.denom="pooled", m.threshold=0.1)$Balanced.mean.diffs[[1]][1])-1)/20

# Store estimated ATT
treat_est <- est_ATT1$coefficients[2]

# Calculate bias, variance, and MSE
bias <- round(((est_ATT1$coefficients[2]-pop_att)/pop_att), digits=3)
var <- (summary(est_ATT1)[13]$coefficients[2,2])^2
mse <- bias^2 + var
}
  
# Return list
list(i, r, ss, xs, prop, cov_diff, treat_est, pop_att, bias, var, mse)  

}
#----------------------------------------------------------------------------------------------#  


#----------------------------------------------------------------------------------------------#
# Step 3: Store and view results
#----------------------------------------------------------------------------------------------#
results <- as.data.frame(sim_function)
list_cols <- sapply(results, is.list)
results <- cbind(results[!list_cols], t(apply(results[list_cols], 1, unlist)))
colnames(results) <- c("iteration", "knn", "sample_size", "covs", "prop_treat",
                     "perc_balance", "est_att", "pop_att","bias", "var", "mse")
# Collapse across all replications for average values
results <- aggregate(. ~ knn + sample_size + covs + prop_treat, data=results, FUN=mean)
results$iteration <- length(sims)
summary(results)
#----------------------------------------------------------------------------------------------#  
