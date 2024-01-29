


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

# Process model over time: our model of population dynamics
for (t in 1:(n.years-1)){
N[1,t+1] <- B[1,t] * rho[1,t] * phi[1,1,t] * eta[1,1,t] + B[2,t] * rho[2,t] * phi[1,2,t] * 
 eta[2,1,t] + B[3,t] * rho[3,t] * phi[1,3,t] * eta[3,1,t] + N[1,t] * phi[2,1,t] * 
 (1 - kappa[1,t]) 
N[2,t+1] <- B[1,t] * rho[1,t] * phi[1,1,t] * eta[1,2,t] + B[2,t] * rho[2,t] * phi[1,2,t] * 
 eta[2,2,t] + B[3,t] * rho[3,t] * phi[1,3,t] * eta[3,2,t] + N[2,t] * phi[2,2,t] * 
 (1 - kappa[2,t]) 
N[3,t+1] <- B[1,t] * rho[1,t] * phi[1,1,t] * eta[1,3,t] + B[2,t] * rho[2,t] * phi[1,2,t] * 
 eta[2,3,t] + B[3,t] * rho[3,t] * phi[1,3,t] * eta[3,3,t] + N[3,t] * phi[2,3,t] *
 (1 - kappa[3,t]) 
B[1,t+1] <- B[1,t] * phi[2,1,t] * nu[1,1,t] + B[2,t] * phi[2,2,t] * nu[2,1,t] + B[3,t] * 
 phi[2,3,t] * nu[3,1,t] + N[1,t] * phi[2,1,t] * kappa[1,t] 
B[2,t+1] <- B[1,t] * phi[2,1,t] * nu[1,2,t] + B[2,t] * phi[2,2,t] * nu[2,2,t] + B[3,t] * 
 phi[2,3,t] * nu[3,2,t] + N[2,t] * phi[2,2,t] * kappa[2,t] 
B[3,t+1] <- B[1,t] * phi[2,1,t] * nu[1,3,t] + B[2,t] * phi[2,2,t] * nu[2,3,t] + B[3,t] * 
 phi[2,3,t] * nu[3,3,t] + N[3,t] * phi[2,3,t] * kappa[3,t] 
}
# Observation model
for (t in 1:n.years){
for (s in 1:3){
C[s,t] ~ dnorm(B[s,t], tau[s])
} #s
} #t

# Multistate capture-recapture data (with multinomial likelihood)
# Define state-transition and re-encounter probabilities
for (t in 1:(n.years-1)){
psi[1,t,1] <- 0
psi[1,t,2] <- 0
psi[1,t,3] <- 0
psi[1,t,4] <- 0
psi[1,t,5] <- 0
psi[1,t,6] <- 0
psi[1,t,7] <- phi[1,1,t] * eta[1,1,t]
psi[1,t,8] <- phi[1,1,t] * eta[1,2,t]
psi[1,t,9] <- phi[1,1,t] * eta[1,3,t]
psi[2,t,1] <- 0
psi[2,t,2] <- 0
psi[2,t,3] <- 0
psi[2,t,4] <- 0
psi[2,t,5] <- 0
psi[2,t,6] <- 0
psi[2,t,7] <- phi[1,2,t] * eta[2,1,t]
psi[2,t,8] <- phi[1,2,t] * eta[2,2,t]
psi[2,t,9] <- phi[1,2,t] * eta[2,3,t]
psi[3,t,1] <- 0
psi[3,t,2] <- 0
psi[3,t,3] <- 0
psi[3,t,4] <- 0
psi[3,t,5] <- 0
psi[3,t,6] <- 0
psi[3,t,7] <- phi[1,3,t] * eta[3,1,t]
psi[3,t,8] <- phi[1,3,t] * eta[3,2,t]
psi[3,t,9] <- phi[1,3,t] * eta[3,3,t]
psi[4,t,1] <- 0
psi[4,t,2] <- 0
psi[4,t,3] <- 0
psi[4,t,4] <- phi[2,1,t] * nu[1,1,t]
psi[4,t,5] <- phi[2,1,t] * nu[1,2,t]
psi[4,t,6] <- phi[2,1,t] * nu[1,3,t]
psi[4,t,7] <- 0
psi[4,t,8] <- 0
psi[4,t,9] <- 0
psi[5,t,1] <- 0
psi[5,t,2] <- 0
psi[5,t,3] <- 0
psi[5,t,4] <- phi[2,2,t] * nu[2,1,t]
psi[5,t,5] <- phi[2,2,t] * nu[2,2,t]
psi[5,t,6] <- phi[2,2,t] * nu[2,3,t]
psi[5,t,7] <- 0
psi[5,t,8] <- 0
psi[5,t,9] <- 0
psi[6,t,1] <- 0
psi[6,t,2] <- 0
psi[6,t,3] <- 0
psi[6,t,4] <- phi[2,3,t] * nu[3,1,t]
psi[6,t,5] <- phi[2,3,t] * nu[3,2,t]
psi[6,t,6] <- phi[2,3,t] * nu[3,3,t]
psi[6,t,7] <- 0
psi[6,t,8] <- 0
psi[6,t,9] <- 0
psi[7,t,1] <- 0
psi[7,t,2] <- 0
psi[7,t,3] <- 0
psi[7,t,4] <- phi[2,1,t] * kappa[1,t]
psi[7,t,5] <- 0
psi[7,t,6] <- 0
psi[7,t,7] <- phi[2,1,t] * (1-kappa[1,t])
psi[7,t,8] <- 0
psi[7,t,9] <- 0
psi[8,t,1] <- 0
psi[8,t,2] <- 0
psi[8,t,3] <- 0
psi[8,t,4] <- 0
psi[8,t,5] <- phi[2,2,t] * kappa[2,t]
psi[8,t,6] <- 0
psi[8,t,7] <- 0
psi[8,t,8] <- phi[2,2,t] * (1-kappa[2,t])
psi[8,t,9] <- 0
psi[9,t,1] <- 0
psi[9,t,2] <- 0
psi[9,t,3] <- 0
psi[9,t,4] <- 0
psi[9,t,5] <- 0
psi[9,t,6] <- phi[2,3,t] * kappa[3,t]
psi[9,t,7] <- 0
psi[9,t,8] <- 0
psi[9,t,9] <- phi[2,3,t] * (1-kappa[3,t])

po[1,t] <- 0
po[2,t] <- 0
po[3,t] <- 0
po[4,t] <- p[1,t]
po[5,t] <- p[2,t]
po[6,t] <- p[3,t]
po[7,t] <- 0
po[8,t] <- 0
po[9,t] <- 0

# Calculate probability of non-encounter (dq) and reshape the array for the encounter
# probabilities
for (s in 1:ns){
dp[s,t,s] <- po[s,t]
dq[s,t,s] <- 1-po[s,t]
} #s
for (s in 1:(ns-1)){
for (m in (s+1):ns){
dp[s,t,m] <- 0
dq[s,t,m] <- 0
} #m
} #s
for (s in 2:ns){
for (m in 1:(s-1)){
dp[s,t,m] <- 0
dq[s,t,m] <- 0
} #m
} #s
} #t

# Define the multinomial likelihood
for (t in 1:((n.years-1)*ns)){
marr[t,1:(n.years*ns-(ns-1))] ~ dmulti(pr[t,], rel[t])
}
# Define the cell probabilities of the multistate m-array
for (t in 1:(n.years-2)){
U[(t-1)*ns+(1:ns), (t-1)*ns+(1:ns)] <- ones
for (j in (t+1):(n.years-1)){
U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-2)*ns+(1:ns)] %*% 
 psi[,t,] %*% dq[,t,]
} #j
} #t
U[(n.years-2)*ns+(1:ns), (n.years-2)*ns+(1:ns)] <- ones
# Diagonal
for (t in 1:(n.years-2)){
pr[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] 
 %*% psi[,t,] %*% dp[,t,]
# Above main diagonal
for (j in (t+1):(n.years-1)){
pr[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] %*% 
 psi[,j,] %*% dp[,j,]
} #j
} #t
pr[(n.years-2)*ns+(1:ns), (n.years-2)*ns+(1:ns)] <- psi[,n.years-1,] %*% dp[,n.years-1,] 
# Below main diagonal
for (t in 2:(n.years-1)){
for (j in 1:(t-1)){
pr[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- zero
} #j
} #t
# Last column: probability of non-recapture
for (t in 1:((n.years-1)*ns)){
pr[t,(n.years*ns-(ns-1))] <- 1-sum(pr[t,1:((n.years-1)*ns)])
} #t
}
")

# Initial values
inits <- function(cc=cormorant$count){
  B <- array(NA, dim=c(3, 14))
  B[1,1] <- rpois(1, cc[1,1])
  B[2,1] <- rpois(1, cc[2,1])
  B[3,1] <- rpois(1, cc[3,1])
  N <- B * 0.55/0.45
  return(list(B=B, N=N))
}

# Parameters monitored
parameters <- c("beta.phi", "phi", "mu.kappa", "kappa", "p", "mean.eta", "mean.nu", "mean.rho", "sigma",
                "N", "B")
# MCMC settings
ni <- 150000; nb <- 50000; nc <- 3; nt <- 100; na <- 3000
# Call JAGS from R (ART 143 min) and check convergence
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
             n.thin=nt, n.adapt=na, parallel=TRUE) 
traceplot(out1)