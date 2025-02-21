#code adapted from:
#Ruckert, K. L., Wong, T. E., Guan, Y., Haran, M., & Applegate, P. J. (n.d.). 
#A Calibration Problem and Markov Chain Monte Carlo. 
#In V. Srikrishnan & K. Keller (Eds.), Advanced Risk Analysis in the Earth Sciences.

################################################################################
#Note from SR:
#Methods kept the same with the exception of a methodological change to the 
#convergence diagnostic that compares uncertainty bounds
#rounded to each number of significant figures to the corresponding significant figure in the mean estimate.
#The original approach used the wrong value for z.

#I also added code showing how to continue a Markov chain where the last one left off
#and computed convergence metrics and plots for the extended chain.
################################################################################

# Clear away any existing variables or figures.
rm(list = ls())
graphics.off()
# Install and read in packages.
# install.packages("coda")
# install.packages("mcmc")
# install.packages("batchmeans")
library(coda)
library(mcmc)
library(batchmeans)
# Set the seed for random sampling.
# All seeds in this tutorial are arbitrary.

setwd("...")

if(!dir.exists("figures")) dir.create("figures")

if(!dir.exists("figures/reading")) dir.create("figures/reading")

set.seed(1)

# Read in some observations with non-correlated measurement error
# (independent and identically distributed assumption).

#the below file does not exist
data <- read.table("observations.txt", header=TRUE)
t <- data$time
observations <- data$observations


# Plot data.
#par(mfrow = c(1,1))
pdf(file="figures/reading/time_vs_observations.pdf",width=7,height=5)  
plot(t, observations, pch = 20, xlab = "Time", ylab = "Observations")
dev.off()

# Set up a simple linear equation as a physical model.
model <- function(parm,t){ # Inputs are parameters and length of data
  model.p <- length(parm) # number of parameters in the physical model
  alpha <- parm[1]
  beta <- parm[2]
  y.mod <- alpha*t + beta # This linear equation represents a simple physical model
  return(list(mod.obs = y.mod, model.p = model.p))
}

# Sample function for calculating the root mean squared error given a set of
# parameters, a vector of time values, and a vector of observations.
fn <- function(parameters, t, obs){
  alpha <- parameters[1]
  beta <- parameters[2]
  data <- alpha*t + beta
  resid <- obs - data
  rmse <- sqrt(mean(resid^2))
  # return the root mean square error
  return(rmse)
}
# Plug in random values for the parameters.
parameter_guess <- c(0.5, 2)
# Optimize the physical model to find initial starting values for parameters.
# For optim to print more information add the arguments:
# method = "L-BFGS-B", control=list(trace=6))
result <- optim(parameter_guess, fn, gr=NULL, t, observations)
start_alpha <- result$par[1]
start_beta <- result$par[2]
parameter <- c(start_alpha, start_beta)
# Use the optimized parameters to generate a fit to the data and
# calculate the residuals.
y.obs <- model(parameter,t)
res <- observations - y.obs$mod.obs
start_sigma <- sd(res)


#par(mfrow = c(1,1))
pdf("figures/reading/time_vs_residuals.pdf",width=7,height=5)  
plot(res, type = "l", ylab = "Residuals", xlab = "Time")
points(res, pch = 20)
abline(h = 0, lty = 2)
dev.off()

# Set up priors.
bound.lower <- c(-1, -1, 0)
bound.upper <- c( 1, 3, 1)

# Name the parameters and specify the number of physical model parameters (alpha and beta).
# sigma is a statistical parameter and will not be counted in the number.
parnames <- c("alpha", "beta", "sigma")
model.p <- 2
# Load the likelihood model for measurement errors
source("iid_obs_likelihood.R")
# Optimize the likelihood function to estimate initial starting values
p <- c(start_alpha, start_beta, start_sigma)
p0 <- c(0.3, 1, 0.6) # random guesses
p0 <- optim(p0, function(p) -log.post(p))$par
print(round(p0,4))

# Set the step size and number of iterations.
step <- c(0.001, 0.01, 0.01)
NI <- 1E3
# Run MCMC calibration.
mcmc.out <- metrop(log.post, p0, nbatch = NI, scale = step)
prechain <- mcmc.out$batch

# Print the acceptance rate as a percent.
acceptrate <- mcmc.out$accept * 100
cat("Accept rate =", acceptrate, "%\n")

# Identify the burn-in period and subtract it from the chains
burnin <- seq(1, 0.01*NI, 1)
mcmc.chains <- prechain[-burnin, ]

## Check #1: Trace Plots:
pdf(file=paste0("figures/reading/traceplots_",NI,"iterations.pdf"),width=10,height=7)  
par(mfrow = c(2,2))
for(i in 1:3){
  plot(mcmc.chains[ ,i], type="l", main = ""
       ,
       ylab = paste('Parameter = ', parnames[i], sep = ''), xlab = "Number of Runs")
}
dev.off()
## Check #2: Monte Carlo Standard Error:
# est: approximation of the mean
# se: estimates the MCMC standard error
bm_est <- bmmat(mcmc.chains)
print(bm_est)

#is the mcse small relative to the size of the estimate?
for(i in 1:length(parnames)){
  print(paste0("ratio of Monte Carlo standard error to Monte Carlo estimate for ",parnames[i],": ", as.numeric(bm_est[i,2]/bm_est[i,1])))
}

# Evaluate the number of significant figures
# assume we want a 95% confidence interval
z <- 1.96
interval <- matrix(data = NA, nrow = length(parnames), ncol = 2,
                   dimnames = list(c(1:length(parnames)), c("lower_bound", "upper_bound")))
for(i in 1:length(parnames)){
  half_width<- z * bm_est[i ,"se"]
  interval[i,1] <- bm_est[i ,"est"] - half_width
  interval[i,2] <- bm_est[i ,"est"] + half_width
}

for(i in 1:length(parnames)){
  print(paste0(parnames[i]," estimate :",as.numeric(bm_est[i,1])))
  print(paste0(parnames[i]," interval :",as.numeric(interval[i,1]),",",as.numeric(interval[i,2])))
}
# "we are strongly suggesting that an estimate of the Monte Carlo standard error
# should be used to assess simulation error and reported. 
# Without an attached MCSE a point estimate should not be trusted" (Flegal, Haran, and Jones, 2008)

## Check #3: Heidelberger and Welch's convergence diagnostic:
heidel.diag(mcmc.chains, eps = 0.1, pvalue = 0.05)


## Check #4: Gelman and Rubin's convergence diagnostic:
set.seed(111)
p0 <- c(0.05, 1.5, 0.6) # Arbitrary choice.
mcmc.out2 <- metrop(log.post, p0, nbatch=NI, scale=step)
prechain2 <- mcmc.out2$batch
set.seed(1708)
p0 <- c(0.1, 0.9, 0.3) # Arbitrary choice.
mcmc.out3 <- metrop(log.post, p0, nbatch=NI, scale=step)
prechain3 <- mcmc.out3$batch
set.seed(1234)
p0 <- c(0.3, 1.1, 0.5) # Arbitrary choice.
mcmc.out4 <- metrop(log.post, p0, nbatch=NI, scale=step)
prechain4 <- mcmc.out4$batch
# The burn-in has already been subtracted from the first chain.
# Thus, the burn-in only needs to be subtracted from the three other
# chains mcmc1 <- mcmc2 <- mcmc3 <- shrink factor

mcmc1<- as.mcmc(mcmc.chains)
mcmc2<- as.mcmc(prechain2[-burnin,])
mcmc3<- as.mcmc(prechain3[-burnin,])
mcmc4<- as.mcmc(prechain4[-burnin,])

set.seed(1)
mcmc_chain_list<- mcmc.list(list(mcmc1,mcmc2,mcmc3,mcmc4))
gelman.diag(mcmc_chain_list)

pdf(file=paste0("figures/reading/gelman_plot_",NI,"iterations.pdf"),width=10,height=7)  
gelman.plot(mcmc_chain_list)
dev.off()

#all point estimates and upper limits of the scale reduction factor are >1.1, so
#this is evidence that the Markov chains have not converged

# Calculate the 90% highest posterior density CI.
# HPDinterval() requires an mcmc object; this was done in the code block above.
hpdi = HPDinterval(mcmc1, prob = 0.90)
# Create density plot of each parameter.

pdf(file=paste0("figures/reading/HPDI_vs_equaltailCI_",NI,"iterations.pdf"),width=10,height=8) 
par(mfrow = c(2,2))
for(i in 1:3){
  # Create density plot.
  p.dens = density(mcmc.chains[,i])
  plot(p.dens, xlab = paste('Parameter ='
                            ,
                            ' ', parnames[i], sep = ''), main="")
  # Add mean estimate.
  abline(v = bm(mcmc.chains[,i])$est, lwd = 2)
  # Add 90% equal-tail CI.
  CI = quantile(mcmc.chains[,i], prob = c(0.05, 0.95))
  lines(x = CI, y = rep(0, 2), lwd = 2)
  points(x = CI, y = rep(0, 2), pch = 16)
  # Add 90% highest posterior density CI.
  lines(x = hpdi[i, ], y = rep(mean(p.dens$y), 2), lwd = 2, col = "red")
  points(x = hpdi[i, ], y = rep(mean(p.dens$y), 2), pch = 16, col = "red")
}
dev.off()
#black is equal tail. red is high posterior density.

################################################################################
#REASSESS CONVERGENCE AFTER ADDING MORE ITERATIONS

NIextra<- NI*20

set.seed(23)
mcmc.out <- metrop(log.post, mcmc1[nrow(mcmc1),], nbatch = NIextra, scale = step)
prechain1 <- mcmc.out$batch

chain1<- rbind(as.matrix(mcmc1),prechain1)

# No more burn-in since we're starting at the last step of the previous chain
mcmc1_all <- as.mcmc(chain1)

set.seed(24)
mcmc.out <- metrop(log.post, mcmc2[nrow(mcmc2),], nbatch = NIextra, scale = step)
prechain2 <- mcmc.out$batch

chain2<- rbind(as.matrix(mcmc2),prechain2)

# No more burn-in since we're starting at the last step of the previous chain
mcmc2_all <- as.mcmc(chain2)

set.seed(25)
mcmc.out<- metrop(log.post, mcmc3[nrow(mcmc3),], nbatch = NIextra, scale = step)
prechain3 <- mcmc.out$batch

chain3<- rbind(as.matrix(mcmc3),prechain3)

# No more burn-in since we're starting at the last step of the previous chain
mcmc3_all <- as.mcmc(chain3)

set.seed(26)
mcmc.out <- metrop(log.post, mcmc4[nrow(mcmc4),], nbatch = NIextra, scale = step)
prechain4 <- mcmc.out$batch

chain4<- rbind(as.matrix(mcmc4),prechain4)

# No more burn-in since we're starting at the last step of the previous chain
mcmc4_all <- as.mcmc(chain4)

## Check #1: Trace Plots:
pdf(file=paste0("figures/reading/traceplots_",NI+NIextra,"iterations.pdf"),width=10,height=7)  
par(mfrow = c(2,2))
for(i in 1:3){
  plot(chain1[ ,i], type="l", main = ""
       ,
       ylab = paste('Parameter = ', parnames[i], sep = ''), xlab = "Number of Runs")
}
dev.off()


set.seed(1)
mcmc_chain_list<- mcmc.list(list(mcmc1_all,mcmc2_all,mcmc3_all,mcmc4_all))
gelman.diag(mcmc_chain_list)

pdf(file=paste0("figures/gelman_plot_",NI+NIextra,"iterations.pdf"),width=10,height=7)  
gelman.plot(mcmc_chain_list)
dev.off()
