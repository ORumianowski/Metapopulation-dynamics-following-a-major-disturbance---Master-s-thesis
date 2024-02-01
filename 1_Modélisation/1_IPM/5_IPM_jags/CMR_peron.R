

library(IPMbook); library(jagsUI)

source("simul_peron.R")

marr <- marray(y, unobs=5)

n.colony = 5
n.age.class = 3
n.years = ncol(y)
ns=n.colony*n.age.class

jags.data <- list(marr=marr, n.years=n.years, rel=rowSums(marr), ns=ns,
                  zero=matrix(0, ncol=ns, nrow=ns), ones=diag(ns),
                  n.colony = n.colony) 

# Write JAGS model file
cat(file = "model1.txt", "
model {

  phi[1] ~ dunif(0, 1) # survie la première année
  phi[2] ~ dunif(0, 1) # Survie des subadultes et adultes
  
  kappa[1] ~ dunif(0, 1) # probabilité de recruitement de la Ronze
  kappa[2] ~ dunif(0, 1) # probabilité de recruitement des autres colonies
  
  for (dep in 1:n.colony){
    for (arr in 1:n.colony){
      eta[dep, arr] ~ dunif(0, 1) # natal dispersal
      nu[dep, arr] ~ dunif(0, 1) # breeding dispersal
    }
  }
  
  for (t in 1:(n.years-1)){
    for (s in 1:n.colony){
      p[s,t] ~ dunif(0, 1)
    } #s
  } #t
  
  # Define the multinomial likelihood
  for (t in 1:((n.years-1)*ns)){
    marr[t,1:(n.years*ns-(ns-1))] ~ dmulti(pr[t,1:(n.years*ns-(ns-1))], rel[t])
  }
  
  
  # Multistate capture-recapture data (with multinomial likelihood)
  # Define state-transition
  
  # Nestling
  
  psi[1,1] <- 0
  psi[1,2] <- 0
  psi[1,3] <- 0
  psi[1,4] <- 0
  psi[1,5] <- 0
  psi[1,6] <- 0
  psi[1,7] <- 0
  psi[1,8] <- 0
  psi[1,9] <- 0
  psi[1,10] <- 0
  psi[1,11] <- phi[1] * eta[1,1]
  psi[1,12] <- phi[1] * eta[1,2]
  psi[1,13] <- phi[1] * eta[1,3]
  psi[1,14] <- phi[1] * eta[1,4]
  psi[1,15] <- phi[1] * eta[1,5]
  
  psi[2,1] <- 0
  psi[2,2] <- 0
  psi[2,3] <- 0
  psi[2,4] <- 0
  psi[2,5] <- 0
  psi[2,6] <- 0
  psi[2,7] <- 0
  psi[2,8] <- 0
  psi[2,9] <- 0
  psi[2,10] <- 0
  psi[2,11] <- phi[1] * eta[2,1]
  psi[2,12] <- phi[1] * eta[2,2]
  psi[2,13] <- phi[1] * eta[2,3]
  psi[2,14] <- phi[1] * eta[2,4]
  psi[2,15] <- phi[1] * eta[2,5]
  
  psi[3,1] <- 0
  psi[3,2] <- 0
  psi[3,3] <- 0
  psi[3,4] <- 0
  psi[3,5] <- 0
  psi[3,6] <- 0
  psi[3,7] <- 0
  psi[3,8] <- 0
  psi[3,9] <- 0
  psi[3,10] <- 0
  psi[3,11] <- phi[1] * eta[3,1]
  psi[3,12] <- phi[1] * eta[3,2]
  psi[3,13] <- phi[1] * eta[3,3]
  psi[3,14] <- phi[1] * eta[3,4]
  psi[3,15] <- phi[1] * eta[3,5]
  
  psi[4,1] <- 0
  psi[4,2] <- 0
  psi[4,3] <- 0
  psi[4,4] <- 0
  psi[4,5] <- 0
  psi[4,6] <- 0
  psi[4,7] <- 0
  psi[4,8] <- 0
  psi[4,9] <- 0
  psi[4,10] <- 0
  psi[4,11] <- phi[1] * eta[4,1]
  psi[4,12] <- phi[1] * eta[4,2]
  psi[4,13] <- phi[1] * eta[4,3]
  psi[4,14] <- phi[1] * eta[4,4]
  psi[4,15] <- phi[1] * eta[4,5]
  
  psi[5,1] <- 0
  psi[5,2] <- 0
  psi[5,3] <- 0
  psi[5,4] <- 0
  psi[5,5] <- 0
  psi[5,6] <- 0
  psi[5,7] <- 0
  psi[5,8] <- 0
  psi[5,9] <- 0
  psi[5,10] <- 0
  psi[5,11] <- phi[1] * eta[5,1]
  psi[5,12] <- phi[1] * eta[5,2]
  psi[5,13] <- phi[1] * eta[5,3]
  psi[5,14] <- phi[1] * eta[5,4]
  psi[5,15] <- phi[1] * eta[5,5]
  
  # Breeders
  
  psi[6,1] <- 0
  psi[6,2] <- 0
  psi[6,3] <- 0
  psi[6,4] <- 0
  psi[6,5] <- 0
  psi[6,6] <- phi[2] * nu[1,1]
  psi[6,7] <- phi[2] * nu[1,2]
  psi[6,8] <- phi[2] * nu[1,3]
  psi[6,9] <- phi[2] * nu[1,4]
  psi[6,10] <- phi[2] * nu[1,5]
  psi[6,11] <- 0
  psi[6,12] <- 0
  psi[6,13] <- 0
  psi[6,14] <- 0
  psi[6,15] <- 0
  
  psi[7,1] <- 0
  psi[7,2] <- 0
  psi[7,3] <- 0
  psi[7,4] <- 0
  psi[7,5] <- 0
  psi[7,6] <- phi[2] * nu[2,1]
  psi[7,7] <- phi[2] * nu[2,2]
  psi[7,8] <- phi[2] * nu[2,3]
  psi[7,9] <- phi[2] * nu[2,4]
  psi[7,10] <- phi[2] * nu[2,5]
  psi[7,11] <- 0
  psi[7,12] <- 0
  psi[7,13] <- 0
  psi[7,14] <- 0
  psi[7,15] <- 0
  
  psi[8,1] <- 0
  psi[8,2] <- 0
  psi[8,3] <- 0
  psi[8,4] <- 0
  psi[8,5] <- 0
  psi[8,6] <- phi[2] * nu[3,1]
  psi[8,7] <- phi[2] * nu[3,2]
  psi[8,8] <- phi[2] * nu[3,3]
  psi[8,9] <- phi[2] * nu[3,4]
  psi[8,10] <- phi[2] * nu[3,5]
  psi[8,11] <- 0
  psi[8,12] <- 0
  psi[8,13] <- 0
  psi[8,14] <- 0
  psi[8,15] <- 0
  
  psi[9,1] <- 0
  psi[9,2] <- 0
  psi[9,3] <- 0
  psi[9,4] <- 0
  psi[9,5] <- 0
  psi[9,6] <- phi[2] * nu[4,1]
  psi[9,7] <- phi[2] * nu[4,2]
  psi[9,8] <- phi[2] * nu[4,3]
  psi[9,9] <- phi[2] * nu[4,4]
  psi[9,10] <- phi[2] * nu[5,5]
  psi[9,11] <- 0
  psi[9,12] <- 0
  psi[9,13] <- 0
  psi[9,14] <- 0
  psi[9,15] <- 0
  
  psi[10,1] <- 0
  psi[10,2] <- 0
  psi[10,3] <- 0
  psi[10,4] <- 0
  psi[10,5] <- 0
  psi[10,6] <- phi[2] * nu[5,1]
  psi[10,7] <- phi[2] * nu[5,2]
  psi[10,8] <- phi[2] * nu[5,3]
  psi[10,9] <- phi[2] * nu[5,4]
  psi[10,10] <- phi[2] * nu[5,5]
  psi[10,11] <- 0
  psi[10,12] <- 0
  psi[10,13] <- 0
  psi[10,14] <- 0
  psi[10,15] <- 0
  
  # Pre-Breeders
  
  psi[11,1] <- 0
  psi[11,2] <- 0
  psi[11,3] <- 0
  psi[11,4] <- 0
  psi[11,5] <- 0
  psi[11,6] <- phi[2] * kappa[1]
  psi[11,7] <- 0
  psi[11,8] <- 0
  psi[11,9] <- 0
  psi[11,10] <- 0
  psi[11,11] <- phi[2] * (1 - kappa[1])
  psi[11,12] <- 0
  psi[11,13] <- 0
  psi[11,14] <- 0
  psi[11,15] <- 0
  
  psi[12,1] <- 0
  psi[12,2] <- 0
  psi[12,3] <- 0
  psi[12,4] <- 0
  psi[12,5] <- 0
  psi[12,6] <- 0 
  psi[12,7] <- phi[2] * kappa[2]
  psi[12,8] <- 0
  psi[12,9] <- 0
  psi[12,10] <- 0
  psi[12,11] <- 0
  psi[12,12] <- phi[2] * (1 - kappa[2])
  psi[12,13] <- 0
  psi[12,14] <- 0
  psi[12,15] <- 0
  
  psi[13,1] <- 0
  psi[13,2] <- 0
  psi[13,3] <- 0
  psi[13,4] <- 0
  psi[13,5] <- 0
  psi[13,6] <- 0 
  psi[13,7] <- 0
  psi[13,8] <- phi[2] * kappa[2]
  psi[13,9] <- 0
  psi[13,10] <- 0
  psi[13,11] <- 0
  psi[13,12] <- 0
  psi[13,13] <- phi[2] * (1 - kappa[2])
  psi[13,14] <- 0
  psi[13,15] <- 0
  
  psi[14,1] <- 0
  psi[14,2] <- 0
  psi[14,3] <- 0
  psi[14,4] <- 0
  psi[14,5] <- 0
  psi[14,6] <- 0 
  psi[14,7] <- 0
  psi[14,8] <- 0
  psi[14,9] <- phi[2] * kappa[2]
  psi[14,10] <- 0
  psi[14,11] <- 0
  psi[14,12] <- 0
  psi[14,13] <- 0
  psi[14,14] <- phi[2] * (1 - kappa[2])
  psi[14,15] <- 0
  
  psi[15,1] <- 0
  psi[15,2] <- 0
  psi[15,3] <- 0
  psi[15,4] <- 0
  psi[15,5] <- 0
  psi[15,6] <- 0 
  psi[15,7] <- 0
  psi[15,8] <- 0
  psi[15,9] <- 0
  psi[15,10] <- phi[2] * kappa[2]
  psi[15,11] <- 0
  psi[15,12] <- 0
  psi[15,13] <- 0
  psi[15,14] <- 0
  psi[15,15] <- phi[2] * (1 - kappa[2])
  
  #  Define re-encounter probabilities
  
  for (t in 1:(n.years-1)){
    
    po[1,t] <- 0
    po[2,t] <- 0
    po[3,t] <- 0
    po[4,t] <- 0
    po[5,t] <- 0
    po[6,t] <- p[1,t]
    po[7,t] <- p[2,t]
    po[8,t] <- p[3,t]
    po[9,t] <- p[4,t]
    po[10,t] <- p[5,t]
    po[11,t] <- 0
    po[12,t] <- 0
    po[13,t] <- 0
    po[14,t] <- 0
    po[15,t] <- 0
    
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
  
  # Define the cell probabilities of the multistate m-array
  
  for (t in 1:(n.years-2)){
    U[(t-1)*ns+(1:ns), (t-1)*ns+(1:ns)] <- ones
    
  for (j in (t+1):(n.years-1)){
    U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-2)*ns+(1:ns)] %*% psi[,] %*% dq[,t,]
   } #j
  } #t
  
  U[(n.years-2)*ns+(1:ns), (n.years-2)*ns+(1:ns)] <- ones
  
  
  # Diagonal
  for (t in 1:(n.years-2)){
    pr[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] %*% psi[,] %*% dp[,t,]
    
  # Above main diagonal
  for (j in (t+1):(n.years-1)){
    pr[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] %*% psi[,] %*% dp[,j,]
    } #j
  } #t
  
  pr[(n.years-2)*ns+(1:ns), (n.years-2)*ns+(1:ns)] <- psi[,] %*% dp[,n.years-1,]
  
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

inits <- function(){
  
  phi=c(0.2,0.8)
  kappa=runif(2, 0.6, 0.9)
  eta=eta
  nu=nu
  p=p
  
  return(list(phi=phi, kappa=kappa, nu=nu, eta=eta))
}

# Parameters monitored
parameters <- c("phi", "kappa", "eta", "nu"#, "sigma",  "rho",
               ) # "N", "B")
# MCMC settings
#ni <- 150000; nb <- 50000; nc <- 3; nt <- 100; na <- 3000
ni <- 1000; nb <- 200; nc <- 2; nt <- 1; na <- 200

# Call JAGS from R (ART 143 min) and check convergence
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
             n.thin=nt, n.adapt=na, parallel=TRUE) 
