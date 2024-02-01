


library(IPMbook); library(nimble)
data(cormorant)
str(cormorant)

# Block survey du modèle
# Dans le modèle de Péron:
# - les paramètres de dispersion et d'effectifs sont fixés à 0 les années sans reproduction
# - les paramètres de dispersion sont changés après un année sans reproduction

# Pour le choix des états initiaux des effectifs des pre-breeding
phi <- c(0.4, 0.8)
kappa <- 0.25
rho <- 1.5
A <- matrix(c(
  phi[2] * (1-kappa), phi[1] * rho,
  phi[2] * kappa, phi[2]), byrow=TRUE, ncol=2)
z <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[,z])
matrix(revec / sum(revec)) # Standardized right eigenvecto


# Bundle data and produce data overview
mydata <- list(C=cormorant$count)

myconsts <- list(n.years=ncol(cormorant$ms.ch))

mycode = nimbleCode(code ={
  
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
  N[1,1] ~ dunif(4270, 7940) # cormorant$count[1,1] * (0.55/0.45) +/- 30%
  B[1,1] ~ dunif(3500, 6500) # cormorant$count[1,1] +/- 30%
  N[2,1] ~ dunif(1710, 3180)
  B[2,1] ~ dunif(1400, 2600)
  N[3,1] ~ dunif(680, 1270)
  B[3,1] ~ dunif(560, 1040)
                 
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
    
})

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
parameters <- c("phi", "kappa", "eta", "nu", "rho", "sigma",
                "N", "B")
# MCMC settings
ni <- 150000; nb <- 50000; nc <- 3; nt <- 100; na <- 3000


nimble_model <- nimbleModel(code = mycode, name = "IPM", constants = myconsts,
                            data = mydata, inits = inits())

Cnimble_model<- compileNimble(nimble_model)

Cnimble_modelConf <- configureMCMC(Cnimble_model, enableWAIC = FALSE, print = FALSE)

MCMC <- buildMCMC(Cnimble_modelConf)
CMCMC <- compileNimble(MCMC)

mcmc.out <- runMCMC(CMCMC,
                    niter = ni, nchains = nc,
                    nburnin = nb, thin = nt,
                    inits = inits, setSeed = TRUE,
                    samples = TRUE, summary = TRUE)

MCMCvis::MCMCtrace(object = mcmc.out[["samples"]],
                   pdf = FALSE, 
                   ind = TRUE)

MCMCvis::MCMCsummary(mcmc.out[["samples"]])








