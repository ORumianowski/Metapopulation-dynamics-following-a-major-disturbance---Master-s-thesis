<<<<<<< HEAD

# ce script a pour ojectif de se rapprocher de la stochasticité

# il y a un pbm au niveau dela nature discrete de N qui necessite un chg de prior - a faire plus tard


library(IPMbook); library(nimble)
data(woodchat5)
str(woodchat5)


marr <- marrayAge(woodchat5$ch, woodchat5$age)

# Bundle data and produce data overview
mydata <- list(marr.j=marr[,,1], marr.a=marr[,,2],
               J=woodchat5$repro[,1],
               age=woodchat5$repro[,3],
               C=woodchat5$count,
               pNinit=dUnif(1, 300))
str(jags.data)

myconsts <- list(rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), n.occasions=dim(marr)[2],len_J = length(J))


mycode = nimbleCode(code ={
  
  # Priors and linear models
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f[1] ~ dunif(0, 10)
  mean.f[2] ~ dunif(0, 10)
  
  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }
  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)
  # Population count data (state-space model)
  # Model for the initial stage-specific population sizes: uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)
  
  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- N[1,t] * mean.f[1] / 2 * mean.sj + N[2,t] * mean.f[2] / 2 * mean.sj 
    N[2,t+1] ~ dbin(mean.sa, (N[1,t] + N[2,t]))
  }
  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }
  # Productivity data (Poisson regression model)
  for (i in 1:len_J){
    J[i] ~ dpois(mean.f[age[i]])
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
  # Annual population growth rate
  for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t]) 
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
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma", "ann.growth.rate", "Ntot")

# MCMC settings
ni <- 2000; nb <- 400; nc <- 3; nt <- 2; na <- 3000

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




=======

# ce script a pour ojectif de se rapprocher de la stochasticité

# il y a un pbm au niveau dela nature discrete de N qui necessite un chg de prior - a faire plus tard


library(IPMbook); library(nimble)
data(woodchat5)
str(woodchat5)


marr <- marrayAge(woodchat5$ch, woodchat5$age)

# Bundle data and produce data overview
mydata <- list(marr.j=marr[,,1], marr.a=marr[,,2],
               J=woodchat5$repro[,1],
               age=woodchat5$repro[,3],
               C=woodchat5$count,
               pNinit=dUnif(1, 300))
str(jags.data)

myconsts <- list(rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), n.occasions=dim(marr)[2],len_J = length(J))


mycode = nimbleCode(code ={
  
  # Priors and linear models
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f[1] ~ dunif(0, 10)
  mean.f[2] ~ dunif(0, 10)
  
  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }
  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)
  # Population count data (state-space model)
  # Model for the initial stage-specific population sizes: uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)
  
  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- N[1,t] * mean.f[1] / 2 * mean.sj + N[2,t] * mean.f[2] / 2 * mean.sj 
    N[2,t+1] ~ dbin(mean.sa, (N[1,t] + N[2,t]))
  }
  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }
  # Productivity data (Poisson regression model)
  for (i in 1:len_J){
    J[i] ~ dpois(mean.f[age[i]])
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
  # Annual population growth rate
  for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t]) 
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
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma", "ann.growth.rate", "Ntot")

# MCMC settings
ni <- 2000; nb <- 400; nc <- 3; nt <- 2; na <- 3000

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




>>>>>>> 2c33ff2863c6f6671e2c82b275eacb5fbac9c0f0
