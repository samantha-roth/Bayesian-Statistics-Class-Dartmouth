################################################################################
# Filename: extraConvergenceChecks.R
# A toy example for implementing a workflow to
# - fit a model to data
# - assess convergence
################################################################################
# List of Contributors:
# Atieh Alipour (atieh.alipour@dartmouth.edu), 
# Klaus Keller (Klaus@dartmouth.edu),
# Haochen Ye (hxy46@psu.edu) 
# Prabhat Hegde (Prabhat.Hegde.TH@dartmouth.edu)
# Sitara Baboolal (Sitara.D.Baboolal@dartmouth.edu)
# Samantha Roth (samantha.m.roth@dartmouth.edu)
################################################################################
# Copyright & License:
# Copyright by the Trustees of Dartmouth College
# Distributed under the GNU general public license
# no warranty
################################################################################
# Sources:
# codes snippets from:
# Ruckert, K. L., Wong, T. E., Guan, Y., Haran, M., & Applegate, P. J. (n.d.). 
# A Calibration Problem and Markov Chain Monte Carlo. 
# In V. Srikrishnan & K. Keller (Eds.), 
# Advanced Risk Analysis in the Earth Sciences.
#
# Ruckert, K. L., Wong, T. E., Guan, Y., & Haran, M. (n.d.). 
# Applying Markov Chain Monte Carlo to Sea-Level Data. 
# In V. Srikrishnan & K. Keller (Eds.), Advanced Risk Analysis in the Earth Sciences.
#
# Ruckert, K. L., Wong, T. E., Lee, B. S., Guan, Y., & Haran, M. 
# (n.d.). Bayesian Inference and Markov Chain Monte Carlo Basics. In 
# V. Srikrishnan &  Applegate, P. J., & Keller, K. (Eds.). (2016). 
# Risk analysis in the Earth Sciences: 
# A Lab manual. 2nd edition. Leanpub. Retrieved from https://leanpub.com/raes 
#
# R manual (accessed via R studio help)
################################################################################
# Versions
# Version / last changes / by whom / what
##########################################
# 1 / Jan 20 2023 / Contributors from above / write initial version
# 2 / Dec 27 2023 / Klaus Keller / update header, clean up code, rename file
# 3 / Jan 25 2024 / Klaus Keller / clean up comments and format
# 4 / Feb 05 2024 / Klaus Keller / add comments
################################################################################
# How to run:
# - save the file in a directory
# - open R-Studio and navigate to the folder with the file
# - make this directory the work directory
# - open the file
# - source the file
################################################################################
# Clear any existing variables and plots. 

rm(list = ls())
graphics.off()
################################################################################

#set your working directory
setwd("/Users/f007f8t/Documents/GitHub/Bayesian-Statistics-Class-Dartmouth")

################################################################################
# Install and load the required libraries

# Install necessary Packages

list.of.packages <- c("BASS","lhs","caTools","Metrics","sensobol","plotrix",
                      "MASS","fields")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()
                                   [,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Loading required libraries

library(BASS)
library(lhs)
library(caTools)
library(Metrics)
library(sensobol)
library(plotrix)
library(MASS)
library(fields)



################################################################################
# Define a seed for reproducibility
set.seed(314)

################################################################################
# Define the toy example, here is a Quadratic function

Y <- function (theta1, theta2,t) {
  return (theta1*(t^2)+theta2*t)}

################################################################################
# First assume we have a true model with known theta1 and theta2
True_theta1 <- -10
True_theta2 <- 2

################################################################################
# Create some observation data with random error from Gaussian distribution
N_obs <- 20
t_obs <- runif(N_obs, -5, 5)
sigma_obs <- 10

Observation<- c()
for (i in 1:length(t_obs)){
  Observation[i] <- Y(True_theta1,True_theta2,t_obs[i])+rnorm(1,0,sigma_obs)
}


################################################################################

# Generate many true values
X <- seq(-5,5,by=0.01)

True_Y<- c()
for (i in 1:length(X)){
  True_Y[i] <- Y(True_theta1,True_theta2,X[i])
}

###############################################################################
# Perform a Bayesian analysis using MCMC

# Define the unnormalized log posterior
logp = function(theta){
  N = N_obs
  
  # Calculate model simulations.
  theta1 = theta[1]
  theta2 = theta[2]
  model = Y(theta1,theta2,t_obs)
  
  # Estimate the residuals (i.e., the deviation of the observations from the 
  # model simulation).
  resid = 	Observation - model
  
  # Get the log of the likelihood function.
  # note this assumes we know sigma
  log.likelihood = -N/2*log(2*pi) - N/2*log(sigma_obs^2) - 1/2*sum(resid^2)/sigma_obs^2
  
  # Use an improper uniform "relatively uninformative" prior.
  log.prior = 0 # log(1)
  
  # Bayesian updating: update the probability estimates based on the observations.
  log.posterior = log.likelihood + log.prior
  
  # Return the unnormalized log posterior value. 
  return(log.posterior)
}


# Set the number of MCMC iterations.
NI = 130000

#Perform the analysis for three different seeds
seeds <-c(311,312,313)


chain1<-mat.or.vec(NI, length(seeds)) 
chain2<-mat.or.vec(NI, length(seeds)) 

for (s in 1:length(seeds)) {
  
  set.seed(seeds[s])
  
  
  # Start with some initial state of parameter estimates, theta^initial
  # Arbitrary choices.
  theta1.init = runif(1, -20, 0) 
  theta2.init = runif(1, -20, 20)  
  theta = c(theta1.init, theta2.init)
  
  # Evaluate the unnormalized posterior of the parameter values
  # P(theta^initial | y)
  lp = logp(theta)
  
  # Setup some variables and arrays to keep track of:
  theta.best = theta         # the best parameter estimates
  lp.max = lp                # the maximum of the log posterior
  theta.new = rep(NA,2)      # proposed new parameters (theta^new)
  accepts = 0                # how many times the proposed new parameters are accepted
  mcmc.chains = array(dim=c(NI,2)) # and a chain of accepted parameters
  
  # Set the standard deviation of the proposal distribution,
  # i.e. the variability of the step size
  prop_sd = c(0.1, 0.1)
  
  # Metropolis-Hastings algorithm MCMC; the proposal distribution proposes the next 
  # point to which the random walk might move. 
  for(i in 1:NI) {
    # Propose a new state (theta^new) based on the current parameter values 
    # theta and the transition probability
    theta.new = rnorm(2, theta, sd = prop_sd)
    
    # Evaluate the new unnormalized posterior of the parameter values
    # and compare the proposed value to the current state
    lp.new = logp(theta.new)
    lq = lp.new - lp
    # Metropolis test; compute the acceptance ratio
    # Draw some uniformly distributed random number 'lr' from [0,1]; 
    lr = log(runif(1))
    
    # If lr < the new proposed value, then accept the parameters setting the 
    # proposed new theta (theta^new) to the current state (theta).
    if(lr < lq) {
      # Update the current theta and log posterior to the new state.
      theta = theta.new
      lp = lp.new
      
      # If this proposed new parameter value is “better” (higher posterior probability)
      # than the current state, then accept with a probability of 1. Hence, increase
      # the number of acceptations by 1.
      accepts = accepts + 1
      
      # Check if the current state is the best, thus far and save if so.
      if(lp > lp.max) {
        theta.best = theta
        lp.max = lp	
      }
    }
    # Append the parameter estimate to the chain. This will generate a series of parameter
    # values (theta_0, theta_1, ...). 
    mcmc.chains[i,] = theta
  }
  
  
  #Checking the acceptance rate
  #SR: Calculate the parameter acceptance rate; it should be between 20 and 50%
  accept.rate <- (accepts/NI) * 100
  print(accept.rate)
  
  chain1[,s] <- mcmc.chains[,1]
  chain2[,s] <- mcmc.chains[,2]
}

#Check the convergence
for(j in 1:NI){
  theta1_seeds_CI <- 2*qt(0.975,length(seeds))*sd(chain1[j,])/sqrt(length(seeds))
  if (theta1_seeds_CI<=0.05){break}}

for(k in 1:NI){
  theta2_seeds_CI <- 2*qt(0.975,length(seeds))*sd(chain2[k,])/sqrt(length(seeds))
  if (theta2_seeds_CI<=0.05){break}}

################################################################################
####################### ADDED BY SAMANTHA ROTH #################################
################################################################################

#Extra tip to evaluate convergence

# Plot the parameters densities from the whole chain vs only the second half

#does it look like the 2 pdfs are very close, i.e. there are enough steps?
pdf(file="allVShalf2_parameter.pdf",8,10)  
par( mfrow= c(2,1))
par(fig=c(0,1,0.5,1),mar=c(4,4,4,4))

plot(density(chain1[(NI/2+1):NI,1]),main="",xlab ="theta1",xlim=c(-12,-9),
     ylim=c(0,1.5),col = "turquoise",lwd=3)
lines(density(chain1[1:NI,1]),col = "purple",lwd=3)
legend("topright", c("second half","full chain"),
       lty=1, lwd = 3, col = c("turquoise","purple"))

plot(density(chain2[(NI/2+1):NI,1]),main="",xlab ="theta2",xlim=c(-3,5),
     ylim=c(0,0.6),col = "turquoise",lwd=3)
lines(density(chain2[1:NI,1]),col = "purple",lwd=3)
legend("topright", c("second half","full chain"),
       lty=1, lwd = 3, col = c("turquoise","purple"))

dev.off()


#does it look like 20k is enough steps?
pdf(file="second10k_vs_first20k_parameter.pdf",8,10)  
par( mfrow= c(2,1))
par(fig=c(0,1,0.5,1),mar=c(4,4,4,4))

plot(density(chain1[1e4:2e4,1]),main="",xlab ="theta1",xlim=c(-12,-9),
     ylim=c(0,1.5),col = "turquoise",lwd=3)
lines(density(chain1[1:2e4,1]),col = "purple",lwd=3)
legend("topright", c("second half","full chain"),
       lty=1, lwd = 3, col = c("turquoise","purple"))

plot(density(chain2[1e4:2e4,1]),main="",xlab ="theta2",xlim=c(-3,5),
     ylim=c(0,0.6),col = "turquoise",lwd=3)
lines(density(chain2[1:2e4,1]),col = "purple",lwd=3)
legend("topright", c("second half","full chain"),
       lty=1, lwd = 3, col = c("turquoise","purple"))

dev.off()

################################################################################
################################################################################
#what if you're interested in estimating a population proportion with mcmc?
#how should you evaluate the reliability of your estimate?
#does monte carlo standard error change with the proportion you're estimating?

#estimating a population proportion with mcmc
library(batchmeans)

mcmc.chains_noburn<- mcmc.chains[-(1:30000),]

thresholds<- seq(-200,-300,by=-10)

bm_est<- rep(NA,length(thresholds))
bm_se<- rep(NA,length(thresholds))

p_hat<- rep(NA,length(thresholds))
sd_p_hat<- rep(NA,length(thresholds))

badXinds<- which(X< -4.7)

pctBelowThresholds<- matrix(NA,nrow=nrow(mcmc.chains_noburn),
                            ncol=length(thresholds))

for(th in 1:length(thresholds)){
  
  Fit_model_MCMC<- mat.or.vec(nrow(mcmc.chains_noburn), length(badXinds))
  for(j in 1:nrow(mcmc.chains_noburn)){
    for (i in badXinds){
      Fit_model_MCMC[j,i] <- Y(mcmc.chains_noburn[j,1],
                               mcmc.chains_noburn[j,2],X[i])+
        rnorm(1,mean=0,sd=sigma_obs)
    }
  }
  
  ind_mat<- matrix(NA,nrow=nrow(Fit_model_MCMC),ncol=ncol(Fit_model_MCMC))
  for(j in 1:ncol(Fit_model_MCMC)){
    ind_mat[,j]<- ifelse(Fit_model_MCMC[,j]< thresholds[th],1,0)
  }
  
  pctBelowThreshold<-rowSums(ind_mat)/ncol(ind_mat)
  pctBelowThresholds[,th]<- pctBelowThreshold
  
  #print(plot(1:length(pctBelowThreshold),pctBelowThreshold,type="l",ylab=paste0("P(Y<",thresholds[th],"|X<-4.7)"),xlab="step"))
  
  #estimate P(Y < -290 | X< -4.7)
  p_hat[th]<- mean(pctBelowThreshold)
  
  #what's the standard deviation of this estimate using the traditional approach
  #for calculating the standard devation of a sample proportion, where
  # sqrt(n)*(p_hat- p)/ sqrt(p(1-p)) ~ N( 0 , 1 )
  sd_p_hat[th]<- sqrt(p_hat[th]*(1-p_hat[th])/length(pctBelowThreshold))
  
  #now use the batchmeans package to estimate the 
  #mean and standard error
  bm_est[th]<- bm(pctBelowThreshold)$est
  bm_se[th]<- bm(pctBelowThreshold)$se
  
}

print(bm_est) #monte carlo estimate of P(Y<threshold | X< -4.7) for each threshold
print(bm_se) #monte carlo standard errors for all thresholds

#how different do the markov chains look for estimating each population proportion?
pdf(file="PropEstimatesConvergence.pdf",width=8,height=6)
par(mar=c(5.1, 4.1, 6.1, 2.1), xpd=TRUE)
plot(1:length(pctBelowThresholds[,1]),pctBelowThresholds[,1],
           type="l", col= "orange", ylim=c(0,1),
           ylab=paste0("P(Y<y|X<-4.7)"),xlab="step")
#lines(1:length(pctBelowThresholds[,3]),pctBelowThresholds[,3],type="l",col="red")
lines(1:length(pctBelowThresholds[,6]),pctBelowThresholds[,6],type="l",col="purple")
#lines(1:length(pctBelowThresholds[,7]),pctBelowThresholds[,7],type="l",col="blue")
lines(1:length(pctBelowThresholds[,10]),pctBelowThresholds[,10],type="l",col="turquoise")
legend("topright", c("y=-200","y=-250","y=-290"),
       lty=1, lwd = 3, col = c("orange","purple","turquoise"),inset=c(0,-0.25))
dev.off()

# plot the observational data, true model and best fit curve
pdf(file="Threshold_vs_PctBelowThreshold.pdf",width=6,height=5) 
plot(thresholds,bm_est,
     ylab="Estimated P(Y< y|X < -4.7)",
     xlab="y")
dev.off()   

# plot the observational data, true model and best fit curve
pdf(file="PctBelowThreshold_vs_BMSE.pdf",width=6,height=5) 
plot(bm_est,bm_se,
     ylab="Monte Carlo Standard Error",
     xlab="Estimated P(Y< y|X < -4.7)")
dev.off()   

pdf(file="Threshold_vs_MCSE.pdf",width=6,height=5) 
plot(thresholds,bm_se,
     ylab="Monte Carlo Standard Error",
     xlab="y")
dev.off()   



