

library(IPMbook); library(nimble)
data(woodchat5)
str(woodchat5)


marr <- marrayAge(woodchat5$ch, woodchat5$age)

#  Bundle data and produce data overview
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2],
                  J=woodchat5$repro[,1],
                  year=woodchat5$repro[,2], age=woodchat5$repro[,3], C=woodchat5$count, pNinit=dUnif(1, 300))
str(jags.data)


myconsts <- list(rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), n.occasions=dim(marr)[2], len_J = length(J))


mycode = nimbleCode(code = { 
  # Priors and linear models
  for (t in 1:(n.occasions-1)){
    
    # effet aléatoire sur les survies annuelles
    logit.sj[t] ~ dnorm(mu.sj, tau.sj)
    sj[t] <- ilogit(logit.sj[t]) # Back-transformation from logit scale
    
    logit.sa[t] ~ dnorm(mu.sa, tau.sa)
    sa[t] <- ilogit(logit.sa[t]) # Back-transformation from logit scale
    
    p[t] <- mean.p
  }
  #log.f = matrix(NaN, nrow = 2, ncol = n.occasions) # a retirer
  #f = matrix(NaN, nrow = 2, ncol = n.occasions) # a retirer
  for (t in 1:n.occasions){
    
    # effet aléatoire sur la fécondité annuelles 
    #log.f[1,t] = 1.5 
      log.f[1,t] ~ dnorm(mu.f[1], tau.f[1])
    f[1,t] <- exp(log.f[1,t]) # Back-transformation from log scale
    
    #log.f[2,t] = 0.5 
    log.f[2,t] ~ dnorm(mu.f[2], tau.f[2])
    f[2,t] <- exp(log.f[2,t]) # Back-transformation from log scale
  }
  mean.sj ~ dunif(0, 1)
  mu.sj <- logit(mean.sj) # Logit transformation
  mean.sa ~ dunif(0, 1)
  mu.sa <- logit(mean.sa) # Logit transformation
  sigma.sj ~ dunif(0, 3)
  tau.sj <- pow(sigma.sj, -2)
  sigma.sa ~ dunif(0, 3)
  tau.sa <- pow(sigma.sa, -2)
  
  for (j in 1:2){
    mean.f[j] ~ dunif(0, 10)
    mu.f[j] <- log(mean.f[j]) # Log transformation
    sigma.f[j] ~ dunif(0, 3)
    tau.f[j] <- pow(sigma.f[j], -2)
  }
  mean.p ~ dunif(0, 1)
  
  # Erreur d'observation
  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)
  
  # Population count data (state-space model)
  # Model for initial stage-spec. population sizes: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)
  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois(N[1,t] * f[1,t] / 2 * sj[t] + N[2,t] * f[2,t] / 2 * sj[t])
    N[2,t+1] ~ dbin(sa[t], N[1,t] + N[2,t])
  }
  
  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }
  # Productivity data (Poisson regression model)
  
  #J = c() # 
  for (i in 1:len_J){
    J[i] ~ dpois(f[age[i], 2])  # J[i] ~ dpois(f[age[i],year[i]]) ## lots of waarnings !!
  }
  # Capture-recapture data (CJS model with multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,1:n.occasions], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,1:n.occasions], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    q[t] <- 1 - p[t] # Probability of non-recapture
    pr.j[t,t] <- sj[t] * p[t]
    pr.a[t,t] <- sa[t] * p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj[t] * prod(sa[(t+1):j]) * prod(q[t:(j-1)]) * p[j]
      pr.a[t,j] <- prod(sa[t:j]) * prod(q[t:(j-1)]) * p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
  }
  # Derived parameters
  # Annual population growth rate (added 0.001 to avoid possible division by 0) 
  for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t] + 0.001)
  }
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
)

# Initial values
inits <- function(){list(mean.sj=runif(1, 0, 0.5))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.f", "mean.p", "sigma.sj", "sigma.sa", "sigma.f",
                "sigma", "sj", "sa", "f", "N", "ann.growth.rate", "Ntot") 

# MCMC settings
ni <- 12000; nb <- 2000; nc <- 3; nt <- 2; na <- 1000

# Run ---------------------------------------------------------------------

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




