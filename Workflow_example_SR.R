################################################################################
# Filename: Workflow_Example.R
# A toy example for implementing a workflow to
# - fit a model to data
# - perform a positive control test
# - perform a method comparison (Frequentist Bootstrap vs Bayesian MCMC 
# - assess convergence and 
# - perform two methods for sensitivity analysis
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
setwd("...")

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

if(!dir.exists("figures")) dir.create("figures")
if(!dir.exists("figures/Workflow_example")) dir.create("figures/Workflow_example")

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
# Plot the real values vs the observation data

# Generate many true values
X <- seq(-5,5,by=0.01)

True_Y<- c()
for (i in 1:length(X)){
  True_Y[i] <- Y(True_theta1,True_theta2,X[i])
}

# Save the plot in a pdf file
pdf(file="figures/Workflow_example/Obs_VS_Model.pdf",10,6.17) 
plot(t_obs,Observation, type="p",col="black", lwd=2,ylim=c(-250,10), xlab="t", 
     ylab="y")
lines(X,True_Y, type="l", col="blue", lwd=3)
legend("bottomright",c("Observation","True model"), lty=c(NA,1),lwd=c(2,2), 
       pch=c(1, NA) ,col=c("black","blue"))
dev.off()   

################################################################################
# Fit a quadratic curve to the observation data using min RMSE

# RMSE
rmse <- function(data,par){
  return(sqrt(mean((data$y -Y(par[1],par[2],data$x))^2)))
}


# Find the parameters
fit_minRMSE <- optim(par=c(-2, 5), fn=rmse, data=data.frame(x=t_obs, 
                                                            y=Observation))


# Predict y values at observational point
Fit_Y_obs<- c()
for (i in 1:length(t_obs)){
  Fit_Y_obs[i] <- Y(fit_minRMSE$par[1],fit_minRMSE$par[2],t_obs[i])
}

# Generate the best fit curve on the obtained parameters
Fit_model_Y<- c()
for (i in 1:length(X)){
  Fit_model_Y[i] <- Y(fit_minRMSE$par[1],fit_minRMSE$par[2],X[i])
}

# plot the observational data, true model and best fit curve
pdf(file="figures/Workflow_example/Best_fit_RMSE.pdf",10,6.17) 
plot(t_obs,Observation, type="p",col="black", lwd=2,ylim=c(-300,30),
     xlim=c(-4.5,4.5), xlab="t", ylab="y")
lines(X,Fit_model_Y, type="l", lty=2, col="red", lwd=3)
lines(X,True_Y, type="l", col="blue", lwd=3)
legend("bottomright",c("Observation","Best fit minimizing RMSE","True model"), 
       col=c("black","red","blue"), lty=c(NA,2,1),lwd=c(2,2,2), pch=c(1,NA,NA) )
dev.off()   

################################################################################
# Perform a Sum of Squared Errors (SSE) Scan

# SSE
SSE_f <- function(data,par){
  return(sum((data$y -Y(par[1],par[2],data$x))^2))
}

# Generate x,y,z for the contour plot
V.theta1 <- seq(-20,0,by=0.1)
V.theta2 <- seq(-20,20,by=0.2)

obs=data.frame(x=t_obs, y=Observation)

#generate a matrix of SSE values for different combinations of theta1 and theta2
V.SSE=mat.or.vec(length(V.theta1), length(V.theta2))
for(i in 1:length(V.theta1 )){
  for (j in 1:length(V.theta2 )){
    V.SSE[i,j] <- SSE_f(obs,c(V.theta1[i],V.theta2[j])) }}

# plot the SSE surface 

pdf(file="figures/Workflow_example/SSE.pdf",12,8) 
par( mfrow= c(2,2))
par(fig=c(0,1,0.4,1),mar=c(4,4,4,4))
contour(V.theta1, V.theta2, 
        V.SSE,levels = ceiling(seq(min(V.SSE),min(V.SSE)+sd(V.SSE),
                                   sd(V.SSE)/10)) ,
        col = hcl.colors(11,"heat") ,lwd=2,xlab="theta1", ylab="theta2")
abline(h = V.theta2[100], lwd = 3,lty=2,col="darkgrey");
abline(v = V.theta1[100], lwd = 3,lty=2,col="darkgrey");
points(True_theta1,True_theta2,pch=4, lwd = 4,col="black")
points(fit_minRMSE$par[1],fit_minRMSE$par[2],pch=4, lwd = 4,col="red")
legend("bottomright", c("SSE-Contours", "theta1-fixed","theta2-fixed",
                        "True-value","Best estimate minimizing RMSE"),
       pch=c(NA,NA,NA,4,4), lty=c(1,2,2,NA,NA),lwd=c(2,2,2,4,4), 
       col = c("darkred","darkgrey","darkgrey","black","red"),)

par(fig=c(0,0.5,0.05,0.47), new=TRUE)
plot(V.theta2,V.SSE[101,], type="l",col="darkred", lwd=2, xlab="theta2", 
     ylab="SSE")
abline(v = True_theta2, lwd = 3,lty=2,col="black");
legend("topright", c("SSE for fixed theta1","True value"),lty=c(1,2),
       lwd=c(2,2), 
       col = c("darkred","black"))

par(fig=c(0.5,1,0.05,0.47), new=TRUE)
plot(V.theta1,V.SSE[,111], type="l",col="darkred", lwd=2, xlab="theta1", 
     ylab="SSE")
abline(v = True_theta1, lwd = 3,lty=2,col="black");
legend("topleft", c("SSE for fixed theta2","True value"),lty=c(1,2),lwd=c(2,2), 
       col = c("darkred","black"))

dev.off()   

#SR changed V.SSE[100,] to V.SSE[101,] and V.SSE[,100] to V.SSE[,111]
################################################################################
# Perform a bootstrap analysis as a "simple" methods check 

# Calculate residuals during observed time series (observation - best fit)
res <- Observation - Fit_Y_obs

#SR: Define number of desired bootstrap samples 
N=10000  

#SR: Get a bootstrap sample from the residuals
white.boot = mat.or.vec(N, N_obs) 
white.boot_sd = rep(NA,N)

for(i in 1:N) {
  #SR: sample with replacement from the residuals
  white.boot[i,] = sample(res,size=N_obs,replace=TRUE)
  white.boot_sd[i] = sd(white.boot[i,])
  }


# Fit  quadratic curves to the bootstrapping samples using min RMSE
boot.par=mat.or.vec(N, 2)
boot.fit=mat.or.vec(N, length(X))

for(i in 1:N) {
  #SR: estimate model parameters again but defining the observation as the 
  #y values predicted with original theta1 and theta2 estimates 
  #plus bootstrapped errors
  fit.new <- optim(par=c(-2, 5), fn=rmse, data=data.frame(x=t_obs, 
                                                          y=Fit_Y_obs+white.boot[i,]))
  
  boot.par[i,1]=fit.new$par[1]
  boot.par[i,2]=fit.new$par[2]
  #SR: generate observations using the new distribution you fit 
  #and noise with the same standard deviation as the bootstrapped residuals 
  for(j in 1:length(X)){
    boot.fit[i,j] <-Y(fit.new$par[1],fit.new$par[2],X[j])+rnorm(1,mean=0,
                                                                sd=white.boot_sd[i])}
}

# Apply residuals to best fit
Fit.boot=mat.or.vec(N, length(X)) 

for(i in 1:N) { 
  #SR: add random noise to the estimates of Y from the original model fit
  Fit.boot[i,]=Fit_model_Y+rnorm(1,mean=0,sd=white.boot_sd[i])
}

# Estimate the confidence intervals

Fit.boot.CI = mat.or.vec(length(X), 2)
Fit.boot.mean = c()

#SR: get the mean and 90% prediction interval for each Y corresponding to each X
for(i in 1:length(X)){  
  Fit.boot.CI[i,] <- quantile(Fit.boot[,i], probs = c(.05, .95))
  Fit.boot.mean[i] <- mean(Fit.boot[,i])
  
}

# plot observations and best fit with superimposed bootstraps
pdf(file="figures/Workflow_example/Bootstrap_replicates.pdf",10,6.17)  
plot(X, Fit_model_Y, col="red", type="l",lty=2,lwd=2,ylim=c(-300,30),
     xlim=c(-4.5,4.5) , xlab="t",ylab="y")

lines(X, Fit.boot.CI[,1], col="grey", lty=2, lwd=2)
lines(X, Fit.boot.CI[,2], col="grey", lty=2,lwd=2)
lines(X, True_Y, col="blue", lwd=2)

points(t_obs,Observation, type="p",xlab="t",
       ylab="y",col="black",lwd=3)
legend("bottomright", c("Observations","True model","Best Fit minimizing RMSE",
                        "5-95% ( CI-Bootstrap fits + Noise )"), 
       lty=c(NA,1,2,2),lwd=c(2,2,2,2), pch=c(1,NA,NA,NA), 
       col=c("black","blue","red","grey"))
dev.off()


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
  # i.e. the step size variability
  step = c(0.1, 0.1)
  
  # Metropolis-Hastings algorithm MCMC; the proposal distribution proposes the next 
  # point to which the random walk might move. 
  for(i in 1:NI) {
    # Propose a new state (theta^new) based on the current parameter values 
    # theta and the transition probability
    theta.new = rnorm(2, theta, sd = step)
    
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

#Check the convergence based on the length of CI for each theta parameter
#note: 0.05 is chosen here as an exmample
#the appropriate length will be problem-specific
for(j in 1:NI){
  theta1_seeds_CI <- 2*qt(0.975,length(seeds))*sd(chain1[j,])/sqrt(length(seeds))
  if (theta1_seeds_CI<=0.05){break}}

for(k in 1:NI){
  theta2_seeds_CI <- 2*qt(0.975,length(seeds))*sd(chain2[k,])/sqrt(length(seeds))
  if (theta2_seeds_CI<=0.05){break}}



# Plot chain members convergence

pdf(file="figures/Workflow_example/chain_parameter.pdf",10,8)  
par( mfrow= c(2,1))

par(fig=c(0,1,0.5,1),mar=c(4,4,4,4))

#generating red color bar 
colors_red <- mat.or.vec(length(seeds)+1, 3) 
for (i in 1:length(seeds)){
  colors_red [i,] <- c(1/i,0,0) }

#Add color orange for the true model
colors_red [i+1,] <- c(1,0.5,0)

# theta1 label's name
theta1_names =c()
for(i in 1:length(seeds)){
  theta1_names[i] <- paste("theta1 estimation using MCMC seed",as.character(i)) 
} 

#Add theta1 true value
theta1_names[i+1] <- "theta1 true value"


plot(1:(min(2*j,j+500)),chain1[1:min(2*j,j+500),1],type="l",xlab="iterations",
     ylab="value",col = rgb(1, 0, 0),ylim=c(-20,0))
for (i in 2:length(seeds)){
  lines(1:min(2*j,j+500),chain1[1:min(2*j,j+500),i],type="l",
        col = rgb(matrix(colors_red[i,], ncol = 3)))}
abline(h = True_theta1, lwd = 3,lty = 1,col="orange");
legend("bottomright", theta1_names,
       lty=1, lwd = 3, col = rgb(colors_red))


par(fig=c(0,1,0.05,0.55), new=TRUE)

#generating blue color bar 
colors_blue<- mat.or.vec(length(seeds)+1, 3) 
for (i in 1:length(seeds)){
  colors_blue [i,] <- c(0,0,1/i) }

#Add color orange for the true model
colors_blue[i+1,] <- c(1,0.5,0)

# theta2 label's name
theta2_names =c()
for(i in 1:length(seeds)){
  theta2_names[i] <- paste("theta2 estimation using MCMC seed",as.character(i)) 
} 

#Add theta2 true value
theta2_names[i+1] <- "theta 2 true value"


plot(1:(min(2*k,k+500)),chain2[1:min(2*k,k+500),1],type="l",
     xlab="iterations",ylab="value",col = rgb(0, 0, 1),ylim=c(-20,20))
for (i in 2:length(seeds)){
  lines(1:min(2*k,k+500),chain2[1:min(2*k,k+500),i],type="l",
        col=rgb(matrix(colors_blue[i,], ncol = 3)))}
abline(h = True_theta2, lwd = 2, lty = 1,col="orange");
legend("bottomright", theta2_names,
       lty=1, lwd = 3, col = rgb(colors_blue))
dev.off() 


# Posterior scan

# Generate x,y,z for the contour plot

V.lp=mat.or.vec(length(V.theta1), length(V.theta2))
for(i in 1:length(V.theta1 )){
  for (j in 1:length(V.theta2 )){
    V.lp[i,j] <- logp(c(V.theta1[i],V.theta2[j])) }}

# Plot the Log likelihood surface 

pdf(file="figures/Workflow_example/posterior.pdf",13,10) 

par( mfrow= c(2,2))

par(fig=c(0,1,0.45,1),mar=c(4,4,4,4))
contour(V.theta1, V.theta2, 
        V.lp,levels = ceiling(seq(max(V.lp)-sd(V.lp),max(V.lp),sd(V.lp)/10)),
        col = hcl.colors(11,"heat") ,lwd=2,xlab="theta1", ylab="theta2")
abline(h = V.theta2[100], lwd = 3,lty=2,col="darkgrey");
abline(v = V.theta1[100], lwd = 3,lty=2,col="darkgrey");
points(True_theta1,True_theta2,pch=4, lwd = 4,col="black")
points(mean(chain1[50000:130000,1]),mean(chain2[50000:130000,1]),
       pch=4, lwd = 4,col="red")
legend("bottomright", c("log posterior contours", "theta1-fixed",
                        "theta2-fixed","True value","Best estimate using mcmc"),
       pch=c(NA, NA,NA,4,4) ,lty=c(1,2,2,NA,NA),lwd=c(2,2,2,4,4), 
       col = c("darkred","darkgrey","darkgrey","black","red"),)

par(fig=c(0,0.5,0.05,0.47), new=TRUE)
plot(V.theta2,V.lp[101,], type="l",col="darkred", lwd=2, xlab="theta2", 
     ylab="log posterior")
abline(v = True_theta2, lwd = 3,lty=2,col="black");
legend("bottomright", c("log posterior for fixed theta1","True value"),
       lty=c(1,2),lwd=c(2,2), col = c("darkred","black"))

par(fig=c(0.5,1,0.05,0.47), new=TRUE)
plot(V.theta1,V.lp[,111], type="l",col="darkred", lwd=2, xlab="theta1", 
     ylab="log posterior")
abline(v = True_theta1, lwd = 3,lty=2,col="black");

legend("bottomleft", c("log posterior for fixed theta2","True value"),
       lty=c(1,2),lwd=c(2,2), col = c("darkred","black"))

dev.off() 
#SR: changed 100 to 101 in line 527, 100 to 111 in line 534

#SR: Generate random realizations from the fit curves 
#using last 10000 MCMC simulations
Fit_model_MCMC<- mat.or.vec(10000, length(X))
for(j in 1:10000){
  for (i in 1:length(X)){
    Fit_model_MCMC[j,i] <- Y(mcmc.chains[length(mcmc.chains)/2-j,1],
                             mcmc.chains[length(mcmc.chains)/2-j,2],X[i])+
      rnorm(1,mean=0,sd=sigma_obs)
  }
}



# Generate the fit curves using mean and 5-95%CI of the last 10000 MCMC simulations

Fit.mcmc.CI = mat.or.vec(length(X), 2)
Fit.mcmc.mean = c()

for(i in 1:length(X)){  
  Fit.mcmc.CI[i,] <- quantile(Fit_model_MCMC[,i], probs = c(.05, .95))
  Fit.mcmc.mean[i] <- mean(Fit_model_MCMC[,i])
  
}


# plot the observational data, true model and best fit curve using MCMC
pdf(file="figures/Workflow_example/Best_fit_MCMC.pdf",10,6.17) 

plot(X, Fit.mcmc.mean, col="red", type="l",lty=2,lwd=2,ylim=c(-300,30),
     xlim=c(-4.5,4.5) , xlab="t",
     ylab="y")

lines(X, Fit.mcmc.CI[,1], col="grey", lty=2, lwd=2)
lines(X, Fit.mcmc.CI[,2], col="grey", lty=2,lwd=2)
lines(X, True_Y, col="blue", lwd=2)

points(t_obs,Observation, type="p",xlab="t",
       ylab="y",col="black",lwd=3)
legend("bottomright", c("Observations","True model",
                        "Best Fit-mean mcmc","5-95% ( CI-mcmc + Noise )"), 
       lty=c(NA,1,2,2),lwd=c(2,2,2,2), pch=c(1,NA,NA,NA), 
       col=c("black","blue","red","grey"))
dev.off()


# Plot the parameters densities for different seeds

pdf(file="figures/Workflow_example/seeds_parameter.pdf",8,10)  
par( mfrow= c(2,1))
par(fig=c(0,1,0.5,1),mar=c(4,4,4,4))

# label's name
names =c()
for(i in 1:length(seeds)){
  names[i] <- paste("seed",as.character(i)) 
} 

#Add true value
names[i+1] <- "True value"

# Plot
plot(density(chain1[50000:130000,1]),main="",xlab ="theta1",xlim=c(-12,-9),
     ylim=c(0,1.5),col = rgb(1, 0, 0))
for(i in 2:length(seeds)) {
  lines(density(chain1[50000:130000,i]),col = rgb(matrix(colors_red[i,], 
                                                         ncol = 3)))}
abline(v = True_theta1, lwd = 2, lty = 1,col="orange");
legend("topright", names,
       lty=1, lwd = 3, col = rgb(colors_red))

par(fig=c(0,1,0.05,0.55), new=TRUE)

plot(density(chain2[50000:130000,1]),main="",xlab ="theta2",xlim=c(-3,5),
     ylim=c(0,0.6),col = rgb(0, 0, 1))
for(i in 2:length(seeds)) {
  lines(density(chain2[50000:130000,i]),col = rgb(matrix(colors_blue[i,], ncol = 3)))}
abline(v = True_theta2, lwd = 2, lty = 1,col="orange");
legend("topright", names,
       lty=1, lwd = 3, col = rgb(colors_blue))

dev.off()

# Compare and plot the parameters densities using bootstrapping and MCMC 
# for the same sample number

#sample from an example MCMC chain

set.seed(314)
chain1_sample=c()
chain2_sample=c()

for(i in 1:length(boot.par[,1])){
  ind <- sample(50000:NI,1) #SR: only sample steps after burn-in
  
  chain1_sample[i] <- chain1[ind,1]
  chain2_sample[i] <- chain2[ind,1]
  
}


pdf(file="figures/Workflow_example/boots_vs_mcmc_parameter.pdf",8,10)  
par( mfrow= c(2,1))
par(fig=c(0,1,0.5,1),mar=c(4,4,4,4))


plot(density(chain1_sample),main="",xlab ="theta1",xlim=c(-12,-9),
     ylim=c(0,1.5),col = "red")
lines(density(boot.par[,1]),col = "orange")
abline(v = True_theta1, lwd = 2, lty = 1,col="black");
legend("topleft", c("theta1 estimation using mcmc",
                    "theta1 estimation using bootstrapping",
                    "theta1 true value"),
       lty=1, lwd = 3, col = c("red","orange","black"))

par(fig=c(0,1,0.05,0.55), new=TRUE)

plot(density(chain2_sample),main="",xlab ="theta2",xlim=c(-3,5),ylim=c(0,0.6),
     col = "blue")
lines(density(boot.par[,2]),col = "cyan")
abline(v = True_theta2, lwd = 2, lty = 1,col="black");
legend("topleft", c("theta2 estimation using mcmc",
                    "theta2 estimation using bootstrapping","theta2 true value"),
       lty=1, lwd = 3, col = c("blue","cyan","black"))

dev.off()



# Plot and compare theta1 and theta2 values for MCMC and bootstrapping for 
# the same sample number
pdf(file="figures/Workflow_example/mcmc_parameters.pdf",8,10)  
par( mfrow= c(2,1))
par(fig=c(0,1,0.5,1),mar=c(4,4,4,4))

plot(chain1_sample,chain2_sample,lwd=1,ylab ="theta2",xlab ="theta1",
     col="blue",pch=20)
points(True_theta1,True_theta2,pch=4, lwd = 4,col="black")
legend("topleft", c("Estimation using mcmc","True value"),
       pch=c(20,4), lty=c(NA,NA),lwd = c(1,3), col = c("blue","black"))

par(fig=c(0,1,0.05,0.55), new=TRUE)
plot(boot.par[,1],boot.par[,2],lwd=1,ylab ="theta2",xlab ="theta1",col="blue",pch=20)
points(True_theta1,True_theta2,pch=4, lwd = 4,col="black")
legend("topleft", c("Estimation using bootstrapping","True value"),
       pch=c(20,4), lty=c(NA,NA),lwd = c(1,3), col = c("blue","black"))
dev.off()



# Compare 2D density plots of parameters estimation using MCMC and bootstrapping
pdf(file="figures/Workflow_example/2Ddensities.pdf",8,10)  
par( mfrow= c(2,1))
par(fig=c(0,1,0.5,1),mar=c(4,4,4,4))
t1 <-kde2d(chain1_sample,chain2_sample,n=200)
image.plot (t1,ylab ="theta2",xlab ="theta1",
            legend.args = list( text = "Density\nMCMC\n",cex = .8),
            col=rev(heat.colors(15)))

par(fig=c(0,1,0.05,0.55), new=TRUE)
t2 <-kde2d(boot.par[,1],boot.par[,2],n=200)
image.plot (t2,ylab ="theta2",xlab ="theta1",  
            legend.args = list( text = "Density\nBootstrapping\n",cex = .8),
            col=rev(heat.colors(15)))
dev.off()

################################################################################
####################### ADDED BY SAMANTHA ROTH #################################
################################################################################

#Extra tip to evaluate convergence

# Plot the parameters densities from the whole chain vs only the second half
# after getting rid of burn in 

#does it look like the 2 pdfs are very close, i.e. there are enough steps?

pdf(file="figures/Workflow_example/overlaid_densities_allVShalf2_parameter.pdf",8,10)  
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


#does it look like 10k is enough steps?
pdf(file="figures/Workflow_example/overlaid_densities_second5k_vs_first10k_parameter.pdf",8,10)  
par( mfrow= c(2,1))
par(fig=c(0,1,0.5,1),mar=c(4,4,4,4))

plot(density(chain1[5e3:1e4,1]),main="",xlab ="theta1",xlim=c(-12,-9),
     ylim=c(0,1.5),col = "turquoise",lwd=3)
lines(density(chain1[1:1e4,1]),col = "purple",lwd=3)
legend("topright", c("second 5k","first 10k"),
       lty=1, lwd = 3, col = c("turquoise","purple"))

plot(density(chain2[5e3:1e4,1]),main="",xlab ="theta2",xlim=c(-3,5),
     ylim=c(0,0.6),col = "turquoise",lwd=3)
lines(density(chain2[1:1e4,1]),col = "purple",lwd=3)
legend("topright", c("second 5k","first 10k"),
       lty=1, lwd = 3, col = c("turquoise","purple"))

dev.off()

################################################################################
################################################################################
#what if you're interested in estimating a population proportion with mcmc?
#how should you evaluate the reliability of your estimate?
#does monte carlo standard error change with the proportion you're estimating?

#estimating a population proportion with mcmc
library(batchmeans)

#mcmc.chains<- mcmc.chains[-(1:30000),]

thresholds<- seq(-200,-300,by=-10)

bm_est<- rep(NA,length(thresholds))
bm_se<- rep(NA,length(thresholds))

p_hat<- rep(NA,length(thresholds))
sd_p_hat<- rep(NA,length(thresholds))

badXinds<- which(X< -4.7)

pctBelowThresholds<- matrix(NA,nrow=nrow(mcmc.chains),
                            ncol=length(thresholds))

for(th in 1:length(thresholds)){
  
  Fit_model_MCMC<- mat.or.vec(nrow(mcmc.chains), length(badXinds))
  for(j in 1:nrow(mcmc.chains)){
    for (i in badXinds){
      Fit_model_MCMC[j,i] <- Y(mcmc.chains[j,1],
                               mcmc.chains[j,2],X[i])+
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
pdf(file="figures/Workflow_example/PropEstimatesConvergence.pdf",width=8,height=6)
par(mar=c(5.1, 4.1, 6.1, 2.1), xpd=TRUE)
plot(1:length(pctBelowThresholds[,2]),pctBelowThresholds[,2],
     type="l", col= "orange", ylim=c(0,1),
     ylab=paste0("P(Y<y|X<-4.7)"),xlab="step")
#lines(1:length(pctBelowThresholds[,3]),pctBelowThresholds[,3],type="l",col="red")
lines(1:length(pctBelowThresholds[,6]),pctBelowThresholds[,6],type="l",col="purple")
#lines(1:length(pctBelowThresholds[,7]),pctBelowThresholds[,7],type="l",col="blue")
lines(1:length(pctBelowThresholds[,10]),pctBelowThresholds[,10],type="l",col="turquoise")
legend("topright", c(paste0("y=",thresholds[2]),
                     paste0("y=",thresholds[6]),
                     paste0("y=",thresholds[10])),
       lty=1, lwd = 3, col = c("orange","purple","turquoise"),inset=c(0,-0.25))
dev.off()

#how early do the markov chains for estimating each population proportion converge?
pdf(file="figures/Workflow_example/PropEstimatesConvergenceFirst5k.pdf",width=8,height=6)
par(mar=c(5.1, 4.1, 6.1, 2.1), xpd=TRUE)
plot(1:length(pctBelowThresholds[1:5e3,2]),pctBelowThresholds[1:5e3,2],
     type="l", col= "orange", ylim=c(0,1),
     ylab=paste0("P(Y<y|X<-4.7)"),xlab="step")
#lines(1:length(pctBelowThresholds[,3]),pctBelowThresholds[,3],type="l",col="red")
lines(1:length(pctBelowThresholds[1:5e3,6]),pctBelowThresholds[1:5e3,6],type="l",col="purple")
#lines(1:length(pctBelowThresholds[,7]),pctBelowThresholds[,7],type="l",col="blue")
lines(1:length(pctBelowThresholds[1:5e3,10]),pctBelowThresholds[1:5e3,10],type="l",col="turquoise")
legend("topright", c(paste0("y=",thresholds[2]),
                     paste0("y=",thresholds[6]),
                     paste0("y=",thresholds[10])),
       lty=1, lwd = 3, col = c("orange","purple","turquoise"),inset=c(0,-0.25))
dev.off()

# How does the population proportion estimate change with the threshold?
pdf(file="figures/Workflow_example/Threshold_vs_PctBelowThreshold.pdf",width=6,height=5) 
plot(thresholds,bm_est,
     ylab="Estimated P(Y< y|X < -4.7)",
     xlab="y")
dev.off()   

# How does the monte carlo standard error change with the population proportion estimate?
pdf(file="figures/Workflow_example/PctBelowThreshold_vs_MCSE.pdf",width=6,height=5) 
plot(bm_est,bm_se,
     ylab="Monte Carlo Standard Error",
     xlab="Estimated P(Y< y|X < -4.7)")
dev.off()   

# How does the monte carlo standard error change with the threshold?
pdf(file="figures/Workflow_example/Threshold_vs_MCSE.pdf",width=6,height=5) 
plot(thresholds,bm_se,
     ylab="Monte Carlo Standard Error",
     xlab="y")
dev.off()   


pctBelowThresholds[1:5e3,6]

#does it look like 10k is enough steps?
pdf(file="figures/Workflow_example/overlaid_densities_second5k_vs_first10k_p2.pdf")  

plot(density(pctBelowThresholds[5e3:1e4,2]),main="",
     xlab =paste0("P(Y<",thresholds[2],"|X<-4.7)"),
     col = "turquoise",lwd=3)
lines(density(pctBelowThresholds[1:1e4,2]),col = "purple",lwd=3)
legend("topright", c("second 5k","first 10k"),
       lty=1, lwd = 3, col = c("turquoise","purple"))

dev.off()

pdf(file="figures/Workflow_example/overlaid_densities_second5k_vs_first10k_p6.pdf")  

plot(density(pctBelowThresholds[5e3:1e4,6]),main="",
     xlab =paste0("P(Y<",thresholds[6],"|X<-4.7)"),
     col = "turquoise",lwd=3)
lines(density(pctBelowThresholds[1:1e4,6]),col = "purple",lwd=3)
legend("topright", c("second 5k","first 10k"),
       lty=1, lwd = 3, col = c("turquoise","purple"))

dev.off()


pdf(file="figures/Workflow_example/overlaid_densities_second5k_vs_first10k_p10.pdf")  

plot(density(pctBelowThresholds[5e3:1e4,10]),main="",
     xlab =paste0("P(Y<",thresholds[10],"|X<-4.7)"),
     col = "turquoise",lwd=3)
lines(density(pctBelowThresholds[1:1e4,10]),col = "purple",lwd=3)
legend("topright", c("second 5k","first 10k"),
       lty=1, lwd = 3, col = c("turquoise","purple"))

dev.off()

#does it look like 130k is enough steps?
pdf(file="figures/Workflow_example/overlaid_densities_half1_vs_all_p2.pdf")  

plot(density(pctBelowThresholds[(NI/2+1):NI,2]),main="",
     xlab =paste0("P(Y<",thresholds[2],"|X<-4.7)"),
     col = "turquoise",lwd=3)
lines(density(pctBelowThresholds[1:NI,2]),col = "purple",lwd=3)
legend("topright", c("second half","all"),
       lty=1, lwd = 3, col = c("turquoise","purple"))

dev.off()

pdf(file="figures/Workflow_example/overlaid_densities_half1_vs_all_p6.pdf")  

plot(density(pctBelowThresholds[(NI/2+1):NI,6]),main="",
     xlab =paste0("P(Y<",thresholds[6],"|X<-4.7)"),
     col = "turquoise",lwd=3)
lines(density(pctBelowThresholds[1:NI,6]),col = "purple",lwd=3)
legend("topright", c("second half","all"),
       lty=1, lwd = 3, col = c("turquoise","purple"))

dev.off()


pdf(file="figures/Workflow_example/overlaid_densities_half1_vs_all_p10.pdf")  

plot(density(pctBelowThresholds[(NI/2+1):NI,10]),main="",
     xlab =paste0("P(Y<",thresholds[10],"|X<-4.7)"),
     col = "turquoise",lwd=3)
lines(density(pctBelowThresholds[1:NI,10]),col = "purple",lwd=3)
legend("topright", c("second half","all"),
       lty=1, lwd = 3, col = c("turquoise","purple"))

dev.off()

#compare the rounded uncertainty bounds to the estimate for each significant figure

z <- half_width <- rep(NA, length(thresholds))
interval <- matrix(data = NA, nrow = length(thresholds), ncol = 2)
rownames(interval)<- thresholds
colnames(interval)<- c("lower_bound", "upper_bound")
for(i in 1:length(thresholds)){
  z[i] <- (mean(pctBelowThresholds[,i]) - bm_est[i])/bm_se[i]
  half_width[i] <- z[i] * bm_se[i ]
  interval[i,1] <- bm_est[i] - half_width[i]
  interval[i,2] <- bm_est[i] + half_width[i]
}
print(interval)
print(bm_est)


################################################################################
# Sobol Sensitivity Analysis on Peak value using the model

# "sensobol" package has a special design for a sobol matrix used for the analysis,
#     N here is the length of the matrix.
# The relationship between N and the sample size 
# (when calculating up to 2nd order indices) is:
#     sample size = N*(d+2+d*(d-1)/2), where d is the number of parameters
# Note we take the nearest integer, hence the actual sample size is an 
# approximation of the list of sample size we consider. 

# parameter range

theta_1_max <- 0
theta_1_min <- -20

theta_2_max <- 20
theta_2_min <- -20


# peak estimation
peak <- function(parameter){
  
  parameter[1]=parameter[1]*(theta_1_max-theta_1_min)+theta_1_min 
  parameter[2]=parameter[2]*(theta_2_max-theta_2_min)+theta_2_min
  output=c()
  tt<- seq(-4.5,4.5,0.1)
  for(i in 1:length(tt)){ 
    output[i]<-Y(parameter[1],parameter[2],tt[i])}
  return(max(output)) }



# Function dimension
d <- 2

# Sample size
set.seed(314)
sample_size <-130000
N <- floor(sample_size/(d+2+d*(d-1)/2))

# The input matrix is generated by the Sobol sequence algorithm, distributed in (0, 1)
mat <- sobol_matrices(N = N, params = c("theta1","theta2"), order = "second")
X_sobol <- split(t(mat), rep(1:dim(mat)[1], each = d))

# The output matrix
Y_sobol <- sapply(X_sobol, peak)

# Sensitivity analysis
S <- sobol_indices(Y=Y_sobol,N=N,params = c("theta1","theta2"),
                   boot=TRUE,R=100,order="second")

# Convergence check: if all total-order sensitivity indices have a 
# 95%CI width within 0.05, then it is converged
# If it was not converged increase sample size
Range <- S$results$high.ci[(d+1):(2*d)]-S$results$low.ci[(d+1):(2*d)]
if (all(Range<=0.05)){print("converged")}


# save the results
sobol_result <- matrix(c(S$results$original[3],S$results$original[4],
                         S$results$original[1],S$results$original[2],
                         S$results$original[5]),nrow=1,ncol=5,byrow = FALSE)


#################################################################################
# Sobol Sensitivity Analysis on Peak value using the BASS emulator
# Note: the Sobol analysis in this section directly uses the function in 
# "BASS" package instead of sensobol package

# Choose sample size
sample_size <- 200

# Sample parameters using Latin hypercube sampling

samples <-  randomLHS(sample_size, 2)

# Define the train input and output
train_input <- samples

train_output <- c()
for(i in 1:sample_size){
  train_output[i]<- peak(train_input[i,])}

#Fit the BASS emulator using the train data
#Users may set related parameters for the emulator structure and MCMC
# see more detailed information in bass() function documentation.
#The emulator consists of a series MCMC samples

bass_model <- bass(train_input, train_output,nmcmc = 10000)


# BASS analysis: calculate Sobol' indices based on the emulator
#Sensitivity indices for each MCMC sample will be calculated.
BASS_sensitivity <- sobol(bass_model)

# Convergence check: if all total-order sensitivity indices have a 
# 95%CI width within 0.05, then it is converged
# If it was not converged increase sample size

Range_S_bass <-c()
Range_S_bass[1] <- 2*qt(0.975,length(BASS_sensitivity$S[,1]))*sd(BASS_sensitivity$S[,1])/sqrt(length(BASS_sensitivity$S[,1]))
Range_S_bass[2] <- 2*qt(0.975,length(BASS_sensitivity$S[,2]))*sd(BASS_sensitivity$S[,2])/sqrt(length(BASS_sensitivity$S[,2]))
Range_S_bass[3] <- 2*qt(0.975,length(BASS_sensitivity$T[,1]))*sd(BASS_sensitivity$T[,1])/sqrt(length(BASS_sensitivity$T[,1]))
Range_S_bass[4] <- 2*qt(0.975,length(BASS_sensitivity$T[,2]))*sd(BASS_sensitivity$T[,2])/sqrt(length(BASS_sensitivity$T[,2]))
Range_S_bass[5] <- 2*qt(0.975,length(BASS_sensitivity$S[,3]))*sd(BASS_sensitivity$S[,3])/sqrt(length(BASS_sensitivity$S[,3]))
if (all(Range_S_bass<=0.05)){print("converged")}


# Save the BASS mean sensitivity indices
BASS_result <- matrix(c(mean(BASS_sensitivity$T[ ,1]),
                        mean(BASS_sensitivity$T[ ,2]),mean(BASS_sensitivity$S[ ,1]),
                        mean(BASS_sensitivity$S[ ,2]),mean(BASS_sensitivity$S[ ,3]))
                      ,nrow=1,ncol=5,byrow = FALSE)

# Plot and save the results from both analyses in a pdf file 

Result_sensitivity <- matrix(c(sobol_result, BASS_result),nrow=2,ncol=5,byrow=TRUE)
colnames(Result_sensitivity ) <- c( "theta1-total-effect", 
                                    "theta2-total-effect","theta1-first-order", 
                                    "theta2-first-order", "theta1/theta2-second-order")

pdf(file="figures/Workflow_example/Sensitivity_Analysis.pdf",14,7)  
par(fig=c(0,1,0.05,1),mar = c(3, 6, 2, 2))
barplot(Result_sensitivity, beside = TRUE,cex.axis=1,
        cex.names=1,col=c("cornsilk4","orange"),ylim = c(0, 0.8),
        las=1,ylab="Sobol Sensitivity")
legend("topright", 
       legend = c("Full model-130000 model runs","BASS-Emulation- 200 model runs"), 
       fill = c("cornsilk4", "orange"))

dev.off()

