
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


# 7.7 MODELS WITH AGE EFFECTS


# Define parameter values
n.occasions <- 10 # Number of capture occasions
marked.j <- rep(200, n.occasions-1) # Annual number of newly marked juveniles
marked.a <- rep(30, n.occasions-1) # Annual number of newly marked adults
phi.juv <- 0.3 # Juvenile annual survival
phi.ad <- 0.65 # Adult annual survival
p <- rep(0.5, n.occasions-1) # Recapture
phi.j <- c(phi.juv, rep(phi.ad, n.occasions-2))
phi.a <- rep(phi.ad, n.occasions-1)


#Define matrices with survival and recapture probabilities

PHI.J <- matrix(0, ncol = n.occasions-1, nrow = sum(marked.j))

for (i in 1:length(marked.j)){  # length(marked.j) correspond au nombre d'occasions de recaptures 
  
  PHI.J[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),  i:(n.occasions-1)] <-
    matrix(rep(phi.j[1:(n.occasions-i)],marked.j[i]), ncol = n.occasions-i, byrow = TRUE)
}

P.J <- matrix(rep(p, sum(marked.j)), ncol = n.occasions-1,
              nrow = sum(marked.j), byrow = TRUE)

PHI.A <- matrix(rep(phi.a, sum(marked.a)), ncol = n.occasions-1,
                nrow = sum(marked.a), byrow = TRUE)
P.A <- matrix(rep(p, sum(marked.a)), ncol = n.occasions-1,
              nrow = sum(marked.a), byrow = TRUE)


# Apply simulation function
CH.J <- simul.cjs(PHI.J, P.J, marked.j)
CH.A <- simul.cjs(PHI.A, P.A, marked.a)



# Create vector with occasion of marking   
get.first <- function(x) min(which(x!=0))
f.j <- apply(CH.J, 1, get.first)
f.a <- apply(CH.A, 1, get.first)


# Create matrices X indicating age classes
x.j <- matrix(NA, ncol = dim(CH.J)[2]-1, nrow = dim(CH.J)[1])
x.a <- matrix(NA, ncol = dim(CH.A)[2]-1, nrow = dim(CH.A)[1])
for (i in 1:dim(CH.J)[1]){
  for (t in f.j[i]:(dim(CH.J)[2]-1)){
    x.j[i,t] <- 2
    x.j[i,f.j[i]] <- 1 # on ne pourrait pas placer cette ligne à l'exrt. de la boucle t ?!!
  } #t
} #i
for (i in 1:dim(CH.A)[1]){
  for (t in f.a[i]:(dim(CH.A)[2]-1)){
    x.a[i,t] <- 2
  } #t
} #i

# We combine the two data sets into a common set.
CH <- rbind(CH.J, CH.A)
f <- c(f.j, f.a)
x <- rbind(x.j, x.a)



# m array -----------------------------------------------------------------

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

marray_res = marray(CH)


# m-array time-dependant model -----------------------------------------------------------



library(nimble)

CJSCode = nimbleCode(code =
                       {
                         # Priors and constraints
                         for (t in 1:(n.occasions-1)){
                           phi[t] ~ dunif(0, 1) # Priors for survival
                           p[t] ~ dunif(0, 1) # Priors for recapture
                         }
                         # Define the multinomial likelihood
                         for (t in 1:(n.occasions-1)){
                           marr[t,1:n.occasions] ~ dmulti(pr[t,1:10], r[t])
                         }
        
                         # Define the cell probabilities of the m-array
                         # Main diagonal
                         for (t in 1:(n.occasions-1)){
                           q[t] <- 1-p[t] # Probability of non-recapture
                           pr[t,t] <- phi[t]*p[t]
                           # Above main diagonal
                           for (j in (t+1):(n.occasions-1)){
                             pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
                           } #j
                           # Below main diagonal
                           for (j in 1:(t-1)){
                             pr[t,j] <- 0
                           } #j
                         } #t
                         # Last column: probability of non-recapture
                         for (t in 1:(n.occasions-1)){
                           pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
                         } #t
                      
                       }#,
                     
              #    dimensions = list(pr = c((n.occasions-1), (n.occasions-1)),
              #                      marr = c((n.occasions-1), (n.occasions)))
)



# Create the m-array from the capture-histories
marr <- marray_res

CJSData <- list(marr = marr)

r = c()
# Calculate the number of birds released each year
for (t in 1:(n.occasions-1)){
  r[t] <- sum(marr[t,1:10])
}

CJSConsts <- list(n.occasions = dim(marr)[2], r = r )

# Initial values
inits <- function(){list(phi = runif(dim(marr)[2]-1, 0, 1),
                         p = runif(dim(marr)[2]-1, 0, 1))}
inits = inits()

# Parameters monitored
parameters <- c("phi", "p")

# MCMC settings
ni <- 10000
nt <- 3
nb <- 5000
nc <- 3

# -------------------------------------------------------------------------
# Version compact

mcmc.out <- nimbleMCMC(code = CJSCode, constants = CJSConsts,
                       data = CJSData, inits = inits,
                       monitors = parameters,
                       niter = ni, nchains = nc, nburnin = nb, thin = nt,
                       summary = TRUE, WAIC = TRUE)

mcmc.out$summary


hist(mcmc.out$samples$chain1[,1], nclass = 30, col = "gray", main = "",
     xlab = "Juvenile survival", ylab = "Frequency")
abline(v = 0.5, col = "red", lwd = 2)







