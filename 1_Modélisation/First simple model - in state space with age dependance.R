
library(nimble)

# Simulated Dataset  ------------------------------------------------------

# Dataset properties : define parameter values

total.marked.j <- 550 # 55000   # Total number of marked juveniles

n.occasions <-  (2023 - 1986) # Number of capture occasions  # p32 thèse Perron
marked.j <- rep( trunc(total.marked.j / n.occasions), n.occasions-1) # Annual number of newly marked juveniles
phi.juv <- 0.243  # Juvenile annual survival # p42 thèse Perron
phi.ad <- 0.844 # Adult annual survival # p42 thèse Perron
p <- rep(0.2, n.occasions-1) # Recapture # p137 Perron

#The estimated detection probabilities (averaged overtime) were 0.08 (0.05; 0.10) and 0.48 (0.41; 0.55) in the low- and high-detectability classes. 

phi.j <- c(phi.juv, rep(phi.ad,n.occasions-2))


# Define function to simulate a capture-history (CH) matrix
simul.cjs <- function(PHI, P, marked){
  n.occasions <- dim(PHI)[2] + 1
  CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
  # Define a vector with the occasion of marking
  mark.occ <- rep(1:length(marked), marked[1:length(marked)]) # marked[1:length(marked)] ==?? marked
  
  # Fill the CH matrix
  for (i in 1:sum(marked)){    # sum(marked) == nombre de marqués
    
    CH[i, mark.occ[i]] <- 1 # Write an 1 at the release occasion
    if (mark.occ[i]==n.occasions) next # normalement pas possible, car marked ne comprend pas les marqués de dernière occasion
    
    for (t in (mark.occ[i]+1):n.occasions){
      # Bernoulli trial: does individual survive occasion?
      sur <- rbinom(1, 1, PHI[i,t-1])
      if (sur==0) break # If dead, move to next individual
      # Bernoulli trial: is individual recaptured?
      rp <- rbinom(1, 1, P[i,t-1])
      if (rp==1) CH[i,t] <- 1
    } #t
  } #i
  return(CH)
}

# Define matrices with survival and recapture probabilities
PHI.J <- matrix(0, ncol = n.occasions-1, nrow = sum(marked.j))
for (i in 1:(length(marked.j)-1)){
  PHI.J[(sum(marked.j[1:i])-
           marked.j[i]+1):sum(marked.j[1:i]),i:(n.occasions-1)] <-
    matrix(rep(phi.j[1:(n.occasions-i)], marked.j[i]),
           ncol = n.occasions-i, byrow = TRUE)
}
P.J <- matrix(rep(p, sum(marked.j)), ncol =
                n.occasions-1, nrow = sum(marked.j), byrow = TRUE) 

# Apply simulation function
CH.J <- simul.cjs(PHI.J, P.J, marked.j)



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
ni <- 600
nt <- 2
nb <- 100
nc <- 2

CJSConsts <- list(f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), x = x)

CJSData <- list(y = CH)

CJSInits = inits()

# -------------------------------------------------------------------------
start_time <- proc.time()
mcmc.out <- nimbleMCMC(code = CJSCode, constants = CJSConsts,
                       data = CJSData, inits = CJSInits,
                       monitors = parameters,
                       niter = ni, nchains = nc, nburnin = nb, thin = nt,
                       summary = TRUE, WAIC = TRUE)
end_time <- proc.time()

elapsed_time <- end_time - start_time
print(elapsed_time)

mcmc.out$summary

# Plot --------------------------------------------------------------------

samples = rbind(mcmc.out$samples$chain1,
                mcmc.out$samples$chain2)

par(mfrow = c(1, 3), las = 1)

hist(samples[,3], nclass = 30, col = "gray", main = "",
     xlab = "mean.p", ylab = "Frequency")
abline(v = p, col = "red", lwd = 2)

hist(samples[,1], nclass = 30, col = "gray", main = "",
     xlab = "mean.phijuv", ylab = "Frequency")
abline(v = phi.juv, col = "red", lwd = 2)

hist(samples[,2], nclass = 30, col = "gray", main = "",
     xlab = "mean.phiad", ylab = "Frequency")
abline(v = phi.ad, col = "red", lwd = 2)


