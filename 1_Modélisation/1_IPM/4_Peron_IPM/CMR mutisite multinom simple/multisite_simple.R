

# Choose constants in simulation
nmarked <- 50 # Number of marked individuals each year and site
nyears <- 6 # Number of years
phiA <- 0.8 # Apparent survival probability when at site A
phiB <- 0.7 # Apparent survival probability when at site B
psiAB <- 0.3 # Probability to move from A to B
psiBA <- 0.5 # Probability to move from B to A
pA <- 0.7 # Recapture probability when at site A
pB <- 0.4 # Recapture probability when at site B

# Determine occasion when an individual first captured and marked
f <- rep(rep(1:(nyears-1), each=nmarked), 2)
nind <- length(f) # Total number of marked individuals
# Construct the transition probability matrix (TPM), above called OMEGA
# Includes the dead state as state number 3
# Departure state in rows (time t), arrival state in columns (t+1)
TPM <- matrix(c(
  phiA * (1-psiAB), phiA * psiAB, 1-phiA,
  phiB * psiBA, phiB * (1-psiBA), 1-phiB,
  0, 0, 1 ), nrow=3, byrow=TRUE)
# Construct the observation probability matrix (OPM), above called THETA
# True state is in rows, observed state is in columns
# Includes nondetection as observation event number 3
OPM <- matrix(c(
  pA, 0, 1-pA,
  0, pB, 1-pB,
  0, 0, 1 ), nrow=3, byrow=TRUE)
# State or ecological process
# Simulate true system state
z <- array(NA, dim=c(nind, nyears)) # Empty alive/dead matrix
# Initial conditions: all individuals alive at f(i)
initial.state <- c(rep(1, nind/2), rep(2, nind/2)) 
for (i in 1:nind){
  z[i,f[i]] <- initial.state[i]
}
set.seed(2) # Initialize the RNGs in R
# Propagate alive/dead process forwards via transition rule (=TPM=OMEGA)
for (i in 1:nind){
  for (t in (f[i]+1):nyears){
    departure.state <- z[i,t-1]
    arrival.state <- which(rmultinom(1,1, TPM[departure.state,])==1)
    z[i,t] <- arrival.state
  } #t
} #i
#z # Not shown, but useful if you do look at this
# Observation process: simulate observations using observation matrix OPM (=THETA)
y <- array(3, dim=c(nind, nyears))
for (i in 1:nind){
  y[i,f[i]] <- z[i,f[i]]
  for (t in (f[i]+1):nyears){
    true.state <- z[i,t-1]
    observed.state <- which(rmultinom(1,1, OPM[true.state,])==1)
    y[i,t] <- observed.state
  } #t
} #i


# Analyse -----------------------------------------------------------------

library(IPMbook); library(nimble)

# Transform to m-array
y[y==3] <- 0
marr <- marray(y)

# Bundle data
# Calculate the number of states
ns <- length(unique(as.numeric(y))) - 1

mydata <- list(marr=marr, nyears=ncol(y), rel=rowSums(marr))

myconsts <- list(ns=ns, zero=matrix(0, ns, ns),
                 ones=diag(ns))

mycode = nimbleCode(code ={
  
  # Priors and linear models
  for (t in 1:(nyears-1)){
    phiA[t] <- mean.phi[1]
    phiB[t] <- mean.phi[2]
    psiAB[t] <- mean.psi[1]
    psiBA[t] <- mean.psi[2]
    pA[t] <- mean.p[1]
    pB[t] <- mean.p[2]
  }
  for (u in 1:2){
    mean.phi[u] ~ dunif(0, 1) # Priors for mean state-specific survival
    mean.psi[u] ~ dunif(0, 1) # Priors for mean transitions
    mean.p[u] ~ dunif(0, 1) # Priors for mean state-specific recapture
  }
  # Define state-transition and re-encounter probabilities
  for (t in 1:(nyears-1)){
    psi[1,t,1] <- phiA[t] * (1-psiAB[t])
    psi[1,t,2] <- phiA[t] * psiAB[t]
    psi[2,t,1] <- phiB[t] * psiBA[t]
    psi[2,t,2] <- phiB[t] * (1-psiBA[t])
    po[1,t] <- pA[t]
    po[2,t] <- pB[t]
  }
  # From here onwards, no changes needed regardless of which model is fitted
  # Calculate probability of non-encounter (dq) and reshape the array for the encounter
  # probabilities
  for (t in 1:(nyears-1)){
    for (s in 1:ns){
      dp[s,t,s] <- po[s,t]
      dq[s,t,s] <- 1-po[s,t]
    } #s
    for (s in 1:(ns-1)){
      for (m in (s+1):ns){
        dp[s,t,m] <- 0
        dq[s,t,m] <- 0
      } #s
    } #m
    for (s in 2:ns){
      for (m in 1:(s-1)){
        dp[s,t,m] <- 0
        dq[s,t,m] <- 0
      } #s
    } #m
  } #t
  
  # Define the multinomial likelihood
  for (t in 1:((nyears-1)*ns)){
    marr[t,1:(nyears *ns-(ns-1))] ~ dmulti(pi[t,1:(nyears *ns-(ns-1))], rel[t])
  }
  # Define cell probabilities of the multistate m-array 
  # Matrix U: product of probabilities of state-transition and non-encounter (needed because 
  # there is no product function for matrix multiplication in JAGS)
  for (t in 1:(nyears-2)){
    U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- ones
    for (j in (t+1):(nyears-1)){
      U[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(j-2)*ns+(1:ns)] %*% psi[,t,] %*%
        dq[,t,]
    } #j
  } #t
  U[(nyears-2)*ns+(1:ns), (nyears-2)*ns+(1:ns)] <- ones
  for (t in 1:(nyears-2)){
    # Diagonal
    pi[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] %*% psi[,t,] %*% 
      dp[,t,]
    # Above main diagonal
    for (j in (t+1):(nyears-1)){
      pi[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] %*% psi[,j,] %*% 
        dp[,j,]
    } #j
  } #t
  pi[(nyears-2)*ns+(1:ns),(nyears-2)*ns+(1:ns)] <- psi[,nyears-1,] %*% dp[,nyears-1,] 
  # Below main diagonal
  for (t in 2:(nyears-1)){
    for (j in 1:(t-1)){
      pi[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- zero
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:((nyears-1)*ns)){
    pi[t,(nyears*ns-(ns-1))] <- 1-sum(pi[t,1:((nyears-1)*ns)])
  } #t
}
)

# Initial values
inits <- function(){list(mean.phi=runif(2, 0, 1))}

# Parameters monitored
parameters <- c("mean.phi", "mean.psi", "mean.p")

# MCMC settings
ni <- 3000; nb <- 1000; nc <- 3; nt <- 1; na <- 1000

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



