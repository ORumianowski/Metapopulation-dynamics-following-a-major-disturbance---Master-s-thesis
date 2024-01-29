
library(nimble)

# Simulated Dataset 

source("simulated_dataset.R")

total.marked.j = 55000
n.occasions = (2023 - 1986)

CH.J = simul_data(total.marked.j = total.marked.j, n.occasions = n.occasions)$CH.J


# Récupération des constantes nécessaires au modèle -----------------------

# Create vector with occasion of marking   
get.first <- function(x) min(which(x!=0))
f.j <- apply(CH.J, 1, get.first)


# Create matrices X indicating age classes
x.j <- matrix(NA, ncol = dim(CH.J)[2]-1, nrow = dim(CH.J)[1])
for (i in 1:dim(CH.J)[1]){
  for (t in f.j[i]:(dim(CH.J)[2]-1)){
    x.j[i,t] <- 2
    x.j[i,f.j[i]] <- 1 # on ne pourrait pas placer cette ligne à l'exrt. de la boucle t ?!!
  } #t
} #i


# We combine the two data sets into a common set.
CH <- rbind(CH.J)
f <- c(f.j)
x <- rbind(x.j)


# Function to create a matrix with information about known latent state z
known.state.cjs <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}


# Model and computational parameters --------------------------------------

CJSCode = nimbleCode(code ={
  
  # Priors and constraints
  for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- beta[x[i,t]]
      p[i,t] <- mean.p
    } #t
  } #i
  
  for (u in 1:2){
    beta[u] ~ dunif(0, 1) # Priors for age-specific survival
  }
  mean.p ~ dunif(0, 1) # Prior for mean recapture
  
  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
    } #t
  } #i
}
)

# Function to create a matrix of initial values for latent state z
# Je ne comprends pa l'opération suivante (p182 du kery 2012):
cjs.init.z <- function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2 <- max(which(ch[i,]==1))
    ch[i,f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), beta = runif(2, 0, 1),
                         mean.p = runif(1, 0, 1))}

# Parameters monitored
parameters <- c("beta", "mean.p")

# MCMC settings
ni <- 500
nt <- 1
nb <- 50
nc <- 2

CJSConsts <- list(f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), x = x)

CJSData <- list(y = CH)

# Run ---------------------------------------------------------------------

# mcmc.out <- nimbleMCMC(code = CJSCode, constants = CJSConsts,
#                        data = CJSData, inits = inits,
#                        monitors = parameters,
#                        niter = ni, nchains = nc, nburnin = nb, thin = nt,
#                        summary = TRUE, WAIC = TRUE)

# Version self-designed

start_time <- proc.time()

CJS <- nimbleModel(code = CJSCode, name = "CJS", constants = CJSConsts,
                   data = CJSData, inits = inits())
CCJS<- compileNimble(CJS)

CJSConf <- configureMCMC(CJS, enableWAIC = FALSE, print = FALSE)

CJSMCMC <- buildMCMC(CJSConf)
CCJSMCMC <- compileNimble(CJSMCMC, project = CJS)

mcmc.out <- runMCMC(CCJSMCMC,
                    niter = ni, nchains = nc,
                    nburnin = nb, thin = nt,
                    inits = inits, setSeed = TRUE,
                    samples = TRUE, summary = TRUE)


# CCJSMCMC$run(niter = ni, time = TRUE)
# CCJSMCMC$getTimes()
# calculateWAIC(CCJSMCMC)


end_time <- proc.time()

elapsed_time <- end_time - start_time
print(elapsed_time)

# Check and plot ----------------------------------------------------------


MCMCvis::MCMCtrace(object = mcmc.out[["samples"]],
                   pdf = FALSE, # no export to PDF
                   ind = TRUE)

MCMCvis::MCMCsummary(mcmc.out[["samples"]])

samples = rbind(mcmc.out$samples$chain1,
                mcmc.out$samples$chain2,
                mcmc.out$samples$chain3)

prob.recap.mean = 0.2
phi.juv <- 0.243  
phi.ad <- 0.844

par(mfrow = c(1, 3), las = 1)

hist(samples[,1], nclass = 30, col = "gray", main = "",
     xlab = "mean.p", ylab = "Frequency")
abline(v = prob.recap.mean, col = "red", lwd = 2)

hist(samples[,3], nclass = 30, col = "gray", main = "",
     xlab = "mean.phijuv", ylab = "Frequency")
abline(v = phi.juv, col = "red", lwd = 2)

hist(samples[,2], nclass = 30, col = "gray", main = "",
     xlab = "mean.phiad", ylab = "Frequency")
abline(v = phi.ad, col = "red", lwd = 2)

