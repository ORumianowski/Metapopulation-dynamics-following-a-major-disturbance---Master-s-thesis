


library(IPMbook); library(nimble)
data(cormorant)
str(cormorant)

marr <- marray(cormorant$ms.ch, unobs=3)

phi <- c(0.4, 0.8)
kappa <- 0.25
rho <- 1.5
A <- matrix(c(
  phi[2] * (1-kappa), phi[1] * rho,
  phi[2] * kappa, phi[2]), byrow=TRUE, ncol=2)
z <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[,z])
matrix(revec / sum(revec)) # Standardized right eigenvecto


jags.data <- list(marr=marr, n.years=ncol(cormorant$ms.ch), rel=rowSums(marr), ns=9,
                  zero=matrix(0, ncol=9, nrow=9), ones=diag(9), C=cormorant$count) 
str(jags.data)

# Write JAGS model file
cat(file = "model1.txt", "
model {
# -------------------------------------------------
# Stages:
# N: not-yet recruited individuals
# B: breeders
# Parameters:
# phi[age, site, time]: survival probability
# eta[departure site, arrival site, time]: natal dispersal
# nu[departure site, arrival site, time]: breeding dispersal
# kappa[site, time]: recruitment probability
# rho[site, time]: productivity
# p[site, time]: recapture probability
# -------------------------------------------------
# Priors and linear models
# Productivity
for (t in 1:n.years){
for (s in 1:3){
rho[s,t] <- mean.rho[s]
} #s
} #t
mean.rho[1] ~ dunif(0, 4)
mean.rho[2] ~ dunif(0, 4)
mean.rho[3] ~ dunif(0, 4)
# Parameters of multistate model
for (t in 1:(n.years-1)){
phi[1,1,t] ~ dunif(0, 1)
logit.phi[1,1,t] <- logit(phi[1,1,t])
logit.phi[2,1,t] <- logit.phi[1,1,t] + beta.phi[2,1]
phi[2,1,t] <- ilogit(logit.phi[2,1,t])
logit.kappa[1,t] <- mu.kappa[1] + mu.kappa[2] * t
kappa[1,t] <- ilogit(logit.kappa[1,t])
for (s in 2:3){
logit.phi[1,s,t] <- logit.phi[1,1,t] + beta.phi[1,s]
phi[1,s,t] <- ilogit(logit.phi[1,s,t])
logit.phi[2,s,t] <- logit.phi[1,1,t] + beta.phi[2,s]
phi[2,s,t] <- ilogit(logit.phi[2,s,t])
logit.kappa[s,t] <- logit.kappa[1,t] + mu.kappa[s+1]
kappa[s,t] <- ilogit(logit.kappa[s,t])
} #s
for (s in 1:3){
p[s,t] ~ dunif(0, 1)
eta[1,s,t] <- mean.eta[1,s]
eta[2,s,t] <- mean.eta[2,s]
eta[3,s,t] <- mean.eta[3,s]
nu[1,s,t] <- mean.nu[1,s]
nu[2,s,t] <- mean.nu[2,s]
nu[3,s,t] <- mean.nu[3,s]
} #s
} #t
# Multinomial logit link to impose the constraint that natal dispersal does 
# only depend on the site of departure
logit.mean.eta[1,2] ~ dnorm(0, 0.01)
logit.mean.eta[1,3] <- logit.mean.eta[1,2]
mean.eta[1,2] <- exp(logit.mean.eta[1,2]) / (1 + exp(logit.mean.eta[1,2]) + 
 exp(logit.mean.eta[1,3]))
mean.eta[1,3] <- exp(logit.mean.eta[1,3]) / (1 + exp(logit.mean.eta[1,2]) + 
 exp(logit.mean.eta[1,3]))
mean.eta[1,1] <- 1 - mean.eta[1,2] - mean.eta[1,3]
logit.mean.eta[2,1] ~ dnorm(0, 0.01)
logit.mean.eta[2,3] <- logit.mean.eta[2,1]
mean.eta[2,1] <- exp(logit.mean.eta[2,1]) / (1 + exp(logit.mean.eta[2,1]) + 
 exp(logit.mean.eta[2,3]))
mean.eta[2,3] <- exp(logit.mean.eta[2,3]) / (1 + exp(logit.mean.eta[2,1]) + 
 exp(logit.mean.eta[2,3]))
mean.eta[2,2] <- 1 - mean.eta[2,1] - mean.eta[2,3]
logit.mean.eta[3,1] ~ dnorm(0, 0.01)
logit.mean.eta[3,2] <- logit.mean.eta[3,1]
mean.eta[3,1] <- exp(logit.mean.eta[3,1]) / (1 + exp(logit.mean.eta[3,1]) + 
 exp(logit.mean.eta[3,2]))
mean.eta[3,2] <- exp(logit.mean.eta[3,2]) / (1 + exp(logit.mean.eta[3,1]) + 
 exp(logit.mean.eta[3,2]))
mean.eta[3,3] <- 1 - mean.eta[3,1] - mean.eta[3,2]
# Multinomial logit link to impose the constraint that breeding dispersal does only depend on the
# site of arrival
logit.mean.nu[1,2] ~ dnorm(0, 0.01)
logit.mean.nu[1,3] ~ dnorm(0, 0.01)
logit.mean.nu[2,1] ~ dnorm(0, 0.01)
mean.nu[1,2] <- exp(logit.mean.nu[1,2]) / (1 + exp(logit.mean.nu[1,2]) + 
 exp(logit.mean.nu[1,3]) + exp(logit.mean.nu[2,1]))
mean.nu[1,3] <- exp(logit.mean.nu[1,3]) / (1 + exp(logit.mean.nu[1,2]) + 
 exp(logit.mean.nu[1,3]) + exp(logit.mean.nu[2,1]))
mean.nu[2,1] <- exp(logit.mean.nu[2,1]) / (1 + exp(logit.mean.nu[1,2]) + 
 exp(logit.mean.nu[1,3]) + exp(logit.mean.nu[2,1]))
mean.nu[1,1] <- 1 - mean.nu[1,2] - mean.nu[1,3]
mean.nu[2,2] <- 1 - mean.nu[2,1] - mean.nu[1,3]
mean.nu[2,3] <- mean.nu[1,3]
mean.nu[3,1] <- mean.nu[2,1]
mean.nu[3,2] <- mean.nu[1,2]
mean.nu[3,3] <- 1 - mean.nu[3,1] - mean.nu[3,2]
for (i in 1:4){
mu.kappa[i] ~ dnorm(0, 0.01)
}
beta.phi[2,1] ~ dnorm(0, 0.01)
beta.phi[1,2] ~ dnorm(0, 0.01)
beta.phi[2,2] ~ dnorm(0, 0.01)
beta.phi[1,3] ~ dnorm(0, 0.01)
beta.phi[2,3] ~ dnorm(0, 0.01)
# Residual (observation) error
for (s in 1:3){
sigma[s] ~ dunif(0.05, 1000)
tau[s] <- pow(sigma[s], -2)
}
# Population count data (state-space model)
# Models for the initial population size: uniform priors
N[1,1] ~ dunif(4270, 7940)
B[1,1] ~ dunif(3500, 6500)
N[2,1] ~ dunif(1710, 3180)
B[2,1] ~ dunif(1400, 2600)
N[3,1] ~ dunif(680, 1270)
B[3,1] ~ dunif(560, 1040)
