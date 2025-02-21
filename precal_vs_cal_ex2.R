#Samantha Roth 2025
#here we compare estimation of model parameters using precalibration
#VS Bayesian calibration with MCMC.

################################################################################
# Clear any existing variables and plots. 

rm(list = ls())
graphics.off()
################################################################################

#set your working directory
#setwd("...")

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
  return (theta1*theta2*t + 5*theta1*cos(t)+theta2*sin(t))}

################################################################################
# First assume we have a true model with known theta1 and theta2
True_theta1 <- -10
True_theta2 <- 2

################################################################################
# Create some observation data with random error from Gaussian distribution
N_obs <- 20
set.seed(20)
t_obs <- runif(N_obs, -5, 5)
sigma_obs <- 20

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


################################################################################
###############################PRECALIBRATION###################################
################################################################################

#Get a lhs sample
#theta_1~ unif(-20, 0) 
#theta_2 ~ unif(-20, 20)  
theta1_range<- c(-20,0)
theta2_range<- c(-20,20)

nSamples<- 10000
nPars<- 2

##100 samples, 4 variables
set.seed(33)
lhs_0_1<- lhs::maximinLHS(nSamples,nPars)
#A 1000 by 2 Latin Hypercube Sample matrix with values uniformly distributed on [0,1]
#transform into samples from the right uniform distributions
lhs_sample<- matrix(0, nSamples, nPars)
lhs_sample[,1]<- qunif(lhs_0_1[,1], min = theta1_range[1], max = theta1_range[2])
lhs_sample[,2]<- qunif(lhs_0_1[,2], min = theta2_range[1], max = theta2_range[2])

colnames(lhs_sample)<- c("theta1","theta2")

################################################################################

#define an interval which model runs must be within to be considered 'behavioral'

behavioral_func<- function(sigma_obs, obs, model_eval){
  return(model_eval > obs - 1.96*sigma_obs & model_eval < obs + 1.96*sigma_obs)
}

################################################################################

#get evaluate the function Y at all lhs samples of (theta_1, theta_2) and 
#all t_obs values

#create a matrix where the row corresponds to the lhs sample 
#and the column corresponds to the t_obs value, then
#evaluate for which lhs samples all function evaluations are within the 
#uncertainty bounds for all observations

in_bds<- matrix(NA, nrow= nrow(lhs_sample), ncol= length(Observation))
for(i in 1:nrow(lhs_sample)){
  for(j in 1:length(t_obs)){
    in_bds[i,j]<- behavioral_func(sigma_obs, Observation[j], Y(lhs_sample[i,1],lhs_sample[i,2],t_obs[j]))
  }
}

tot_in_bds<- rowSums(in_bds)

good_inds<- which(tot_in_bds==length(Observation))
good_samples<- lhs_sample[good_inds,]

#estimate the expected values of theta_1 and theta_2

e_theta1<- mean(lhs_sample[good_inds,1])
e_theta2<- mean(lhs_sample[good_inds,2])


################################################################################

#perform precalibration using sizes of lhs samples increased by 10k at each step
#tracking convergence of expected values for theta1 and theta2

maxSamples<- 1.3e5
nAdd<- 10000

e_theta1_chain<- e_theta1
e_theta2_chain<- e_theta2

test_ss<- seq(nSamples+nAdd,maxSamples,by=nAdd)

################################################################################
#THE BELOW CODE CAN BE COMMENTED OUT TO SAVE TIME 
#AND THE RESULTS CAN BE LOADED INSTEAD

for(ss in test_ss){
  new_lhs_0_1 <- augmentLHS(lhs_0_1, m = nAdd)  # Add 5 new samples to the existing design
  lhs_new<- matrix(0, nrow(new_lhs_0_1), nPars)
  #convert the samples from U(0,1) distributions to the parameters' priors
  lhs_new[,1]<- qunif(new_lhs_0_1[,1], min = -20, max = 0)
  lhs_new[,2]<- qunif(new_lhs_0_1[,2], min = -20, max = 20)
  
  #create a matrix that indicates whether the predicted value is in the 
  #uncertainty bounds for each lhs sample (row) and each observation (column)
  in_bds<- matrix(NA, nrow= nrow(lhs_new), ncol= length(Observation))
  for(i in 1:nrow(lhs_new)){
    for(j in 1:length(Observation)){
      in_bds[i,j]<- behavioral_func(sigma_obs, Observation[j], Y(lhs_new[i,1],lhs_new[i,2],t_obs[j]))
    }
  }
  #get the total number of observations for which the model prediction is in bounds
  #for each lhs sample
  tot_in_bds<- rowSums(in_bds) 
  
  #identify the lhs sample for which model predictions are within the 
  #uncertainty bounds for each observation
  good_inds<- which(tot_in_bds==length(Observation))
  good_samples<- lhs_new[good_inds,]
  
  #estimate the expected values of theta_1 and theta_2
  e_theta1_new<- mean(lhs_new[good_inds,1])
  e_theta2_new<- mean(lhs_new[good_inds,2])
  
  #create a chain of expected values of theta1 and theta2 
  #using the accepted lhs samples
  e_theta1_chain<- c(e_theta1_chain,e_theta1_new)
  e_theta2_chain<- c(e_theta2_chain,e_theta2_new)
  
  #update the lhs sample that will be added to in the next step 
  lhs_0_1<- new_lhs_0_1
  
  print(ss)
}


#save precalibration data
if(!dir.exists("data")) dir.create("data")
if(!dir.exists("data/precal_ex2")) dir.create("data/precal_ex2")

save(e_theta1_chain,file="data/precal_ex2/e_theta1_chain")
save(e_theta2_chain,file="data/precal_ex2/e_theta2_chain")

save(good_inds,file="data/precal_ex2/good_inds")
save(lhs_new,file="data/precal_ex2/lhs_new")

#LOAD DATA BELOW IF YOU WANT TO SKIP RUNNING PRECALIBRATION
load("data/precal_ex2/e_theta1_chain")
load("data/precal_ex2/e_theta2_chain")

load("data/precal_ex2/good_inds")
load("data/precal_ex2/lhs_new")

################################################################################

#prepare precalibration data for plotting 
#by mapping accepted (theta1,theta2) samples to the LHS sample size considered

all_ss_precal<- c(nSamples,test_ss)
all_good_ss<- vector()
all_good_theta1<- vector()
all_good_theta2<- vector()
for(ss in all_ss_precal){
  in_bds<- matrix(NA, nrow= ss, ncol= length(Observation))
  for(i in 1:ss){
    for(j in 1:length(Observation)){
      in_bds[i,j]<- behavioral_func(sigma_obs, Observation[j], Y(lhs_new[i,1],lhs_new[i,2],t_obs[j]))
    }
  }
  
  tot_in_bds<- rowSums(in_bds)
  
  good_inds<- which(tot_in_bds==length(Observation))
  good_samples<- lhs_new[good_inds,]
  #all_good_samples[[which(all_ss_precal==ss)]]<- good_samples 
  print(ss)
  
  if(is.matrix(good_samples)){
    all_good_theta1<- c(all_good_theta1,good_samples[,1])
    all_good_theta2<- c(all_good_theta2,good_samples[,2])
  }
  if(is.vector(good_samples)){
    all_good_theta1<- c(all_good_theta1,good_samples[1])
    all_good_theta2<- c(all_good_theta2,good_samples[2])
  }

  all_good_ss<- c(all_good_ss,rep(ss,length(good_inds)))
}

################################################################################

#compare convergence of theta estimates from precalibration 
#to convergence of estimates of theta from Bayesian calibration

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
NI = 1.3e5

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

#save calibration data
if(!dir.exists("data")) dir.create("data")
if(!dir.exists("data/cal_ex2")) dir.create("data/cal_ex2")

save(chain1,file="data/cal_ex2/chain1")
save(chain2,file="data/cal_ex2/chain2")


#LOAD DATA BELOW IF YOU WANT TO SKIP RUNNING CALIBRATION
load("data/cal_ex2/chain1")
load("data/cal_ex2/chain2")

################################################################################
#plot showing convergence of estimates for theta1 and theta2 using 
#Bayesian calibration with MCMC vs precalibration

if(!dir.exists("figures")) dir.create("figures")
if(!dir.exists("figures/precal_vs_cal_ex2")) dir.create("figures/precal_vs_cal_ex2")

pdf(file="figures/precal_vs_cal_ex2/precal_vs_cal_theta1.pdf",10,8)  

plot(1:NI,chain1[,1],type="l",main="theta1",xlab="iterations",
     ylab="value",col = "red",ylim=c(-20,0),lwd = 3)
#lines(c(nSamples,test_ss),e_theta1_chain,type="l",col = "blue")
points(all_good_ss,all_good_theta1,col="blue")
abline(h = True_theta1, lwd = 3,lty = 1,col="black");
legend("topright", c("mcmc","precalibration","truth"),
       lty=1, lwd = 3, col = c("red","blue","black"))

dev.off()

pdf(file="figures/precal_vs_cal_ex2/precal_vs_cal_theta2.pdf",10,8)

plot(1:NI,chain2[,1],type="l",main="theta2",xlab="iterations",
     ylab="value",col = "red",ylim=c(-10,10),lwd = 3)
#lines(c(nSamples,test_ss),e_theta1_chain,type="l",col = "blue")
points(all_good_ss,all_good_theta2,col="blue")
abline(h = True_theta2, lwd = 3,lty = 1,col="black");
legend("topright", c("mcmc","precalibration","truth"),
       lty=1, lwd = 3, col = c("red","blue","black"))

dev.off()

################################################################################

#plot the posterior densities of theta1 and theta2 as estimated by 
#mcmc vs precalibration

pdf(file="figures/precal_vs_cal_ex2/overlaid_densities_cal_vs_precal.pdf",8,10)  
par( mfrow= c(2,1))
par(fig=c(0,1,0.5,1),mar=c(4,4,4,4))

plot(density(chain1[,1]),main="",xlab ="theta1",xlim=c(-20,0),
     col = "purple",lwd=3)
abline(v = all_good_theta1[which(all_good_ss==NI)],col = "turquoise",lwd=3)
abline(v = True_theta1, lwd = 3,lty = 1,col="black")
legend("topright", c("mcmc","precalibration","truth"),
       lty=1, lwd = 3, col = c("purple","turquoise","black"))

plot(density(chain2[,1]),main="",xlab ="theta2",xlim=c(-10,10),
     col = "purple",lwd=3)
abline(v = all_good_theta2[which(all_good_ss==NI)],col = "turquoise",lwd=3)
abline(v = True_theta2, lwd = 3,lty = 1,col="black")
legend("topright", c("mcmc","precalibration","truth"),
       lty=1, lwd = 3, col = c("purple","turquoise","black"))

dev.off()

# Bayesian calibration converges first
plot(density(all_good_theta1[which(all_good_ss==2e4)]))
lines(density(all_good_theta1[which(all_good_ss==1e4)]),col="red")

plot(density(chain1[1:2e4,1]))
lines(density(chain1[1000:1e4,1]),col="red")


plot(density(all_good_theta1[which(all_good_ss==6e4)]))
lines(density(all_good_theta1[which(all_good_ss==1.2e5)]),col="red")
