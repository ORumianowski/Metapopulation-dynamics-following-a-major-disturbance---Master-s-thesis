
library(nimble)

# Simulated Dataset  ------------------------------------------------------

# Dataset properties : define parameter values

total.marked.j <- 55000 # 55000   # Total number of marked juveniles

n.occasions <-  (2023 - 1986) # Number of capture occasions  # p32 thèse Perron
marked.j <- rep( trunc(total.marked.j / n.occasions), n.occasions-1) # Annual number of newly marked juveniles
marked.a <- rep(1, n.occasions-1) # Annual number of newly marked adults
phi.juv <- 0.243  # Juvenile annual survival # p42 thèse Perron
phi.ad <- 0.844 # Adult annual survival # p42 thèse Perron
p <- rep(0.2, n.occasions-1) # Recapture # p137 Perron

#The estimated detection probabilities (averaged overtime) were 0.08 (0.05; 0.10) and 0.48 (0.41; 0.55) in the low- and high-detectability classes. 

phi.j <- c(phi.juv, rep(phi.ad,n.occasions-2))
phi.a <- rep(phi.ad, n.occasions-1)


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
PHI.A <- matrix(rep(phi.a, sum(marked.a)), ncol = n.occasions-1,
                nrow = sum(marked.a), byrow = TRUE)
P.A <- matrix(rep(p, sum(marked.a)), ncol = n.occasions-1,
              nrow = sum(marked.a), byrow = TRUE)

# Apply simulation function
CH.J <- simul.cjs(PHI.J, P.J, marked.j)
CH.A <- simul.cjs(PHI.A, P.A, marked.a)


# Convert CH to m-array ---------------------------------------------------

# Function to create a m-array based on capture-histories (CH)
marray <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow =
                      n.occasions)
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] -
      sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

# Split CH according age-class --------------------------------------------

cap <- apply(CH.J, 1, sum)
ind <- which(cap >= 2)
CH.J.R <- CH.J[ind,] # Juvenile CH recaptured at least once
CH.J.N <- CH.J[-ind,] # Juvenile CH never recaptured

# Remove first capture
first <- numeric()
for (i in 1:dim(CH.J.R)[1]){
  first[i] <- min(which(CH.J.R[i,]==1))
}
CH.J.R1 <- CH.J.R
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R1[i,first[i]] <- 0
}

# Add grown-up juveniles to adults and create m-array
CH.A.m <- rbind(CH.A, CH.J.R1)
CH.A.marray <- marray(CH.A.m)

# Create CH matrix for juveniles, ignoring subsequent recaptures
second <- numeric()
for (i in 1:dim(CH.J.R1)[1]){
  second[i] <- min(which(CH.J.R1[i,]==1))
}
CH.J.R2 <- matrix(0, nrow = dim(CH.J.R)[1], ncol = dim(CH.J.R)[2])
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R2[i,first[i]] <- 1
  CH.J.R2[i,second[i]] <- 1
}

# Create m-array for these
CH.J.R.marray <- marray(CH.J.R2)

# The last column ought to show the number of juveniles not recaptured again and should all be zeros, since all of them are released as adults
CH.J.R.marray[,dim(CH.J)[2]] <- 0

# Create the m-array for juveniles never recaptured and add it to the previous m-array
CH.J.N.marray <- marray(CH.J.N)
CH.J.marray <- CH.J.R.marray + CH.J.N.marray


# Model and computational parameters --------------------------------------

CJSCode = nimbleCode(code =
                       {
                         # Priors and constraints
                         for (t in 1:(n.occasions-1)){
                           phi.juv[t] <- mean.phijuv
                           phi.ad[t] <- mean.phiad
                           p[t] <- mean.p
                         }
                         mean.phijuv ~ dunif(0, 1) # Prior for mean juv. survival
                         mean.phiad ~ dunif(0, 1) # Prior for mean ad. survival
                         mean.p ~ dunif(0, 1) # Prior for mean recapture
                         # Define the multinomial likelihood
                         for (t in 1:(n.occasions-1)){
                           marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,1:n.occasions], r.j[t])
                           marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,1:n.occasions], r.a[t])
                         }
                         # Define the cell probabilities of the m-arrays
                         # Main diagonal
                         for (t in 1:(n.occasions-1)){
                           q[t] <- 1-p[t] # Probability of non-recapture
                           pr.j[t,t] <- phi.juv[t]*p[t]
                           pr.a[t,t] <- phi.ad[t]*p[t]
                           # Above main diagonal
                           for (j in (t+1):(n.occasions-1)){
                             pr.j[t,j] <- phi.juv[t]*prod(phi.ad[(t+1):j])*prod(q[t:
                                                                                    (j-1)])*p[j]
                             pr.a[t,j] <- prod(phi.ad[t:j])*prod(q[t:(j-1)])*p[j]
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
                         } #t 
                       }
)



r.j = c()
r.a = c()
# Calculate the number of birds released each year
for (t in 1:(n.occasions-1)){
  r.j[t] <- sum(CH.J.marray[t,])
  r.a[t] <- sum(CH.A.marray[t,])
}

# Constants
CJSConsts <- list(n.occasions = dim(CH.J.marray)[2], r.j = r.j,  r.a = r.a)

# Bundle data
CJSData <- list(marr.j = CH.J.marray, marr.a = CH.A.marray,
                  n.occasions = dim(CH.J.marray)[2])
# Initial values
inits <- function(){list(mean.phijuv = runif(1, 0, 1), mean.phiad =
                           runif(1, 0, 1), mean.p = runif(1, 0, 1))}
inits = inits()

# Parameters monitored
parameters <- c("mean.phijuv", "mean.phiad", "mean.p")

# MCMC settings
ni <- 3000
nt <- 3
nb <- 1000
nc <- 3


# Run ---------------------------------------------------------------------

mcmc.out <- nimbleMCMC(code = CJSCode, constants = CJSConsts,
                       data = CJSData, inits = inits,
                       monitors = parameters,
                       niter = ni, nchains = nc, nburnin = nb, thin = nt,
                       summary = TRUE, WAIC = TRUE)

mcmc.out$summary


# Plot --------------------------------------------------------------------

samples = rbind(mcmc.out$samples$chain1,
                mcmc.out$samples$chain2,
                mcmc.out$samples$chain3)

par(mfrow = c(1, 3), las = 1)

hist(samples[,1], nclass = 30, col = "gray", main = "",
     xlab = "mean.p", ylab = "Frequency")
abline(v = p, col = "red", lwd = 2)

hist(samples[,3], nclass = 30, col = "gray", main = "",
     xlab = "mean.phijuv", ylab = "Frequency")
abline(v = phi.juv, col = "red", lwd = 2)

hist(samples[,2], nclass = 30, col = "gray", main = "",
     xlab = "mean.phiad", ylab = "Frequency")
abline(v = phi.ad, col = "red", lwd = 2)

