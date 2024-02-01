


library(IPMbook); library(jagsUI)


# Hypotheses - Population at first year -----------------------------------
# We are calculating the proportion of pre breeders versus breeders at equilibrium
phi <- c(0.213, 0.860)
kappa <- 0.7
rho <- 1
A <- matrix(c(
  phi[2] * (1-kappa), phi[1] * rho,
  phi[2] * kappa, phi[2]), byrow=TRUE, ncol=2)
z <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[,z])
matrix(revec / sum(revec)) 
# 0.2 0.8 

# We determine the numbers of pre breeders and breeders for the priors, the first year
pop_init = data.frame(estim = c(2500, 900, 1, 466, 3376)) %>% 
  mutate(B_inf = estim*0.7+1)%>% 
  mutate(B_sup = estim*1.3+1)%>% 
  mutate(N_inf = estim*(0.2/0.8)*0.7+1)%>% 
  mutate(N_sup = estim*(0.2/0.8)*1.3+1) %>% 
  round() %>% 
  as.matrix()

# The colony 3 is problematic with one pair reported, the first year
pop_init[3,2] = 1
pop_init[3,3] = 10
pop_init[3,4] = 1
pop_init[3,5] = 10


# Survey - dataset --------------------------------------------------------

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
  5000,455, 01, 01, 1356) %>% as.integer(), ncol= 5, byrow = T) %>% 
  t() 

# Bundle data 
survey_data = eff
n.colony = nrow(eff)
n.years = ncol(survey_data)
n.age.class = 3
ns=n.colony*n.age.class

jags.data <- list(C=survey_data,
                  n.colony=n.colony,
                  n.years=n.years,
                  pop_init=pop_init) 

# Write JAGS model file
cat(file = "model1.txt", "
model {
 # -------------------------------------------------
  # Stages:
  # N: not-yet recruited individuals
  # B: breeders
  # Parameters:
  # phi[age]: survival probability
  # eta[departure site, arrival site, time]: natal dispersal
  # nu[departure site, arrival site, time]: breeding dispersal
  # kappa[site]: recruitment probability - two classes (LaRonze or not)
  # rho[site]: productivity - two classes (LaRonze or not)
  # p[site, time]: recapture probability # for CMR only
  # -------------------------------------------------
  # Priors 
  
  rho[1] ~ dunif(0, 4) # productivité de la Ronze
  rho[2] ~ dunif(0, 4) # productivité des autres colonies
  
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
  
  
  # Population count data (state-space model)
  # Models for the initial population size: uniform priors
  for (col in 1:n.colony){
    B[col,1] ~ dunif(pop_init[col,2], pop_init[col,3]) 
    N[col,1] ~ dunif(pop_init[col,4], pop_init[col,5]) 
  }

                 
  # Process model over time: our model of population dynamics
  for (t in 1:(n.years-1)){
    # La Ronze
  N[1,t+1] <- B[1,t] * rho[1] * phi[1] * eta[1,1] + 
              B[2,t] * rho[2] * phi[1] * eta[2,1] + 
              B[3,t] * rho[2] * phi[1] * eta[3,1] +
              B[4,t] * rho[2] * phi[1] * eta[4,1] + 
              B[5,t] * rho[2] * phi[1] * eta[5,1] +
              N[1,t] * phi[2] * (1 - kappa[1]) 

  B[1,t+1] <- B[1,t] * phi[2] * nu[1,1] +
              B[2,t] * phi[2] * nu[2,1] + 
              B[3,t] * phi[2] * nu[3,1] + 
              B[4,t] * phi[2] * nu[4,1] + 
              B[5,t] * phi[2] * nu[5,1] + 
              N[1,t] * phi[2] * kappa[1] 
  
  # Satellite colonies
  for (s in 2:n.colony){
    
  N[s,t+1] <- B[1,t] * rho[1] * phi[1] * eta[1,s] + 
              B[2,t] * rho[2] * phi[1] * eta[2,s] + 
              B[3,t] * rho[2] * phi[1] * eta[3,s] + 
              B[4,t] * rho[2] * phi[1] * eta[4,s] + 
              B[5,t] * rho[2] * phi[1] * eta[5,s] + 
              N[s,t] * phi[2] * (1 - kappa[2]) 
              
  B[s,t+1] <- B[1,t] * phi[2] * nu[1,s] + 
              B[2,t] * phi[2] * nu[2,s] +
              B[3,t] * phi[2] * nu[3,s] + 
              B[4,t] * phi[2] * nu[4,s] +
              B[5,t] * phi[2] * nu[5,s] + 
              N[s,t] * phi[2] * kappa[2]
    }
  }
  
  # Residual (observation) error
  
  sigma[1] ~ dunif(0.05, 1000) # La ronze
  sigma[2] ~ dunif(0.05, 1000) # Satellites colonies
  
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
inits <- function(){
  
  return(list())
}

# Parameters monitored
parameters <- c("phi", "kappa", "eta", "nu", "rho", "sigma",
                "N", "B")
# MCMC settings
#ni <- 150000; nb <- 50000; nc <- 3; nt <- 100; na <- 3000
ni <- 1000; nb <- 100; nc <- 3; nt <- 10; na <- 300

# Call JAGS from R and check convergence
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
             n.thin=nt, n.adapt=na, parallel=TRUE) 
#traceplot(out1)

whiskerplot(out1,parameters=c('B'))
