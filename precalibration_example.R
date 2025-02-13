
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


################################################################################
###############################PRECALIBRATION###################################
################################################################################

#Get a lhs sample
#theta_1~ unif(-20, 0) 
#theta_2 ~ unif(-20, 20)  
theta1_range<- c(-20,0)
theta2_range<- c(-20,20)

nSamples<- 1000
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
#and the column corresponds to the t_obs value

Y_evals_lhs<- matrix(NA, nrow= nrow(lhs_sample), ncol= length(t_obs))

for(i in 1:nrow(lhs_sample)){
  for(j in 1:length(t_obs)){
    Y_evals_lhs[i,j]<- Y(lhs_sample[i,1],lhs_sample[i,2],t_obs[j])
  }
}

################################################################################
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

maxSamples<- 5e6
nAdd<- 1000

e_theta1_chain<- e_theta1
e_theta2_chain<- e_theta2

e_theta1_old<- e_theta1
e_theta2_old<- e_theta2

test_ss<- seq(nSamples+nAdd,maxSamples,by=nAdd)

for(ss in test_ss){
  new_lhs_0_1 <- augmentLHS(lhs_0_1, m = nAdd)  # Add 5 new samples to the existing design
  lhs_new<- matrix(0, nrow(new_lhs_0_1), nPars)
  lhs_new[,1]<- qunif(new_lhs_0_1[,1], min = -20, max = 0)
  lhs_new[,2]<- qunif(new_lhs_0_1[,2], min = -20, max = 20)
  
  in_bds<- matrix(NA, nrow= nrow(lhs_new), ncol= length(t_obs))
  for(i in 1:nrow(lhs_new)){
    for(j in 1:length(t_obs)){
      in_bds[i,j]<- behavioral_func(sigma_obs, t_obs[j], Y(lhs_new[i,1],lhs_new[i,2],t_obs[j]))
    }
  }
  
  tot_in_bds<- rowSums(in_bds)
  
  good_inds<- which(tot_in_bds==length(t_obs))
  good_samples<- lhs_new[good_inds,]
  
  #estimate the expected values of theta_1 and theta_2
  
  e_theta1_new<- mean(lhs_new[good_inds,1])
  e_theta2_new<- mean(lhs_new[good_inds,2])
  
  e_theta1_chain<- c(e_theta1_chain,e_theta1_new)
  e_theta2_chain<- c(e_theta2_chain,e_theta2_new)
  
  if(abs(e_theta1_new-e_theta1_old)<(theta1_range[2]-theta1_range[1])/1000 &
     abs(e_theta2_new-e_theta2_old)<(theta2_range[2]-theta2_range[1])/1000){
    break
  } 
  else{
    lhs_0_1<- new_lhs_0_1
    e_theta1_old<- e_theta1_new
    e_theta2_old<- e_theta2_new
  }
}

