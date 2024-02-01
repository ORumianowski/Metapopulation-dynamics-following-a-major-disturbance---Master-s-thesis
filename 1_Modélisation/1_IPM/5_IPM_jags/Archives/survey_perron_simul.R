
library(IPMbook)
library(jagsUI)
library(dplyr)


eff=matrix( c(
  2500,900, 01, 466, 3376,
  2500,1000, 01, 480, 2823,
  2500,1200, 01, 434, 2820,
  2500,1000, 01, 685, 2518,
  2500,01, 01, 410, 1552,
  2500,500, 245, 610, 2077,
  2500,01, 01, 413, 2099,
  3000,700, 01, 84, 2077,
  3000,700, 01, 01, 1275,
  3000,500, 01, 469, 1499,
  3000,280, 210, 560, 1000,
  5000,560, 119, 179, 859,
  5000,272, 196, 672, 894,
  5000,295, 150, 530, 1133,
  5000,280, 42, 406, 1087,
  5000,525, 01, 283, 1512,
  5000,01, 10, 213, 866,
  5000,450, 01, 127, 1502,
  5000,420, 01, 01, 1109,
  5000,455, 01, 01, 1356), ncol= 5, byrow = T) %>% 
  t()



# Hypothèses pour l'initialisation ----------------------------------------


phi <- c(0.213, 0.860)
kappa <- 0.7
rho <- 1
A <- matrix(c(
  phi[2] * (1-kappa), phi[1] * rho,
  phi[2] * kappa, phi[2]), byrow=TRUE, ncol=2)
z <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[,z])
matrix(revec / sum(revec)) # Standardized right eigenvector
# 0.2 0.8

pop_init = data.frame(estim = c(2500, 900, 1, 466, 3376)) %>% 
  mutate(B_inf = estim*0.7)%>% 
  mutate(B_sup = estim*1.3)%>% 
  mutate(N_inf = estim*(0.2/0.8)*0.7)%>% 
  mutate(N_sup = estim*(0.2/0.8)*1.3) %>% 
  round() %>% 
  as.matrix()


n.colony = 5
n.age.class = 3
n.years = ncol(eff)
ns=n.colony*n.age.class

jags.data <- list(C=eff, 
                  ns=ns, n.colony = n.colony,
                  n.years = n.years,
                  pop_init=pop_init) 




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
  # Priors 
  
  rho[1] ~ dunif(0, 4) # productivité de la Ronze
  rho[2] ~ dunif(0, 4) # productivité des autres colonies
  
  phi[1] ~ dunif(0, 1) # survie la première année
  phi[2] ~ dunif(0, 1) # Survie des subadultes et adultes
  
  kappa[1] ~ dunif(0, 1) # probabilité de recruitement de la Ronze
  kappa[2] ~ dunif(0, 1) # probabilité de recruitement des autres colonies
  
  for (dep in 1:3){
    for (arr in 1:3){
      eta[dep, arr] ~ dunif(0, 1) # natal dispersal
      nu[dep, arr] ~ dunif(0, 1) # breeding dispersal
    }
  }
  
  
  # Population count data (state-space model)
  # Models for the initial population size: uniform priors
  for (col in 1:n.colony){
    N[col,1] ~ dunif(pop_init[col,4], pop_init[col,5]) 
    B[col,1] ~ dunif(pop_init[col,2], pop_init[col,3]) 
  }
                 
  # Process model over time: our model of population dynamics
  for (t in 1:(n.years-1)){
    # La Ronze
  N[1,t+1] <- B[1,t] * rho[1] * phi[1] * eta[1,1] + 
              B[2,t] * rho[2] * phi[1] * eta[2,1] + 
              B[3,t] * rho[2] * phi[1] * eta[3,1] +
              N[1,t] * phi[2] * (1 - kappa[1]) 

  B[1,t+1] <- B[1,t] * phi[2] * nu[1,1] +
              B[2,t] * phi[2] * nu[2,1] + 
              B[3,t] * phi[2] * nu[3,1] + 
              N[1,t] * phi[2] * kappa[1] 
  
  # Les satellites
  N[2,t+1] <- B[1,t] * rho[1] * phi[1] * eta[1,2] + 
              B[2,t] * rho[2] * phi[1] * eta[2,2] + 
              B[3,t] * rho[2] * phi[1] * eta[3,2] + 
              N[2,t] * phi[2] * (1 - kappa[2]) 
  
  N[3,t+1] <- B[1,t] * rho[1] * phi[1] * eta[1,3] + 
              B[2,t] * rho[2] * phi[1] * eta[2,3] + 
              B[3,t] * rho[2] * phi[1] * eta[3,3] + 
              N[3,t] * phi[2] * (1 - kappa[2]) 
  
  B[2,t+1] <- B[1,t] * phi[2] * nu[1,2] + 
              B[2,t] * phi[2] * nu[2,2] +
              B[3,t] * phi[2] * nu[3,2] + 
              N[2,t] * phi[2] * kappa[2] 
  
  B[3,t+1] <- B[1,t] * phi[2] * nu[1,3] + 
              B[2,t] * phi[2] * nu[2,3] + 
              B[3,t] * phi[2] * nu[3,3] + 
              N[3,t] * phi[2] * kappa[2] 
  }
  
  # Residual (observation) error
  
  sigma[1] ~ dunif(0.05, 1000) # La ronze
  sigma[2] ~ dunif(0.05, 1000) # Les autres
  
  for (t in 1:n.years){
    
    R[1,t] <- sigma[1] * B[1,t] 
    tau[1,t] <- pow(R[1,t], -2)
    
    for (sat in 2:3){
      R[sat,t] <- sigma[2] * B[sat,t] 
      tau[sat,t] <- pow(R[sat,t], -2)
    }
  }
  
  # Observation model
  for (t in 1:n.years){
    for (s in 1:3){
      C[s,t] ~ dnorm(B[s,t], tau[s,t])
    } #s
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
parameters <- c("phi", "kappa", "eta", "nu", "mean.rho", "sigma",
                "N", "B")
# MCMC settings
#ni <- 150000; nb <- 50000; nc <- 3; nt <- 100; na <- 3000
ni <- 30000; nb <- 5000; nc <- 3; nt <- 10; na <- 3000

# Call JAGS from R (ART 143 min) and check convergence
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
             n.thin=nt, n.adapt=na, parallel=TRUE) 
traceplot(out1)
