


library(IPMbook); library(nimble)
data(cormorant)
str(cormorant)

# CMR survey du modèle

# Bundle data and produce data overview

# Transform to m-array

marr <- marray(cormorant$ms.ch, unobs=3)

mydata <- list(marr=marr, rel=rowSums(marr))

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
  
  # Define the multinomial likelihood
  for (t in 1:((n.years-1)*ns)){
    marr[t,1:(n.years*ns-(ns-1))] ~ dmulti(pr[t,], rel[t])
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
parameters <- c("phi", "kappa", "eta", "nu", "sigma", # "rho",
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








