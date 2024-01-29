
set.seed(1248)

# State code
# 1 : Nestling born in colony 1
# 2 : Nestling born in colony 2
# 3 : Nestling born in colony 3
# 4 : Nestling born in colony 4
# 5 : Nestling born in colony 5
# 6 : Breeders in colony 1
# 7 : Breeders in colony 2
# 8 : Breeders in colony 3
# 9 : Breeders in colony 4
# 10 : Pre-Breeders in colony 5
# 11 : Pre-Breeders in colony 1
# 12 : Pre-Breeders in colony 2
# 13 : Pre-Breeders in colony 3
# 14 : Pre-Breeders in colony 4
# 15 : Pre-Breeders in colony 5
# 16 : DEAD

# Constants in simulation

n.state = 16
n.colony = 5 # Number of colonies
nyears = 40 # Number of years
nmarked = 50 # Number of marked individuals each year and site

# Apparent survival probability
# - age class dependence
phi1 = 0.6 # first year
phiA = 0.85 # subadult and adult

# Natal dispersion from one colony to another
# - colony dependence
eta_ = matrix(runif(n = n.colony*n.colony , min = 0, max = 1), 
             nrow = n.colony, ncol = n.colony) 
eta = eta_ / rowSums(eta_)

# Breeding dispersion from one colony to another
# - colony dependence
nu_ = matrix(runif(n = n.colony*n.colony , min = 0, max = 1), 
              nrow = n.colony, ncol = n.colony) 
nu = nu_ / rowSums(nu_)

# Recruitment probability
kappaLR = 0.7 # LaRonze
kappaSAT = 0.8 # Satellite colonies
kappa = c(kappaLR, kappaSAT)

# Productivity # not used in CMR
#rhoLR = 3 # LaRonze
#rhoSAT = 2 # Satellite colonies

# recapture probability
#- time and colony dependence
p = matrix(runif(n = n.colony*n.colony , min = 0.2, max = 0.5),
             nrow = nyears, ncol = n.colony) 

# Determine occasion when an individual first captured and marked
f = rep(rep(1:(nyears-1), each=nmarked), n.colony)
nind = length(f) # Total number of marked individuals


# Construct the transition probability matrix (TPM),
# Includes the dead state as state number 6
# Departure state in rows (time t), arrival state in columns (t+1)


# Devenir des nestlings : ligne poussin de TPM
TMP_nestling = function(birth_colony, n.colony, phi1, eta){
  
  next_state = matrix(0, nrow = 1, ncol = n.state)
  for (next_col in 1:n.colony){
    next_state[1,(2*n.colony+next_col)] = phi1 * eta[birth_colony,next_col] # nestling become prebreeders 
  }
  next_state[1,n.state] = 1 - phi1 # become dead
  return(next_state)
}

# Devenir des breeders : ligne breeders de TPM
TMP_breeder = function(former_colony, n.colony, phiA, nu){
  
  next_state = matrix(0, nrow = 1, ncol = n.state)
  for (next_col in 1:n.colony){
    next_state[1,(1*n.colony+next_col)] = phiA * nu[former_colony,next_col] # breeders stay breeders
  }
  next_state[1,n.state] = 1 - phiA # become dead
  return(next_state)
}

# Devenir des prebreeders : ligne prebreeders de TPM
TMP_prebreeder = function(settlement_colony, n.colony, phiA, kappa){
  
  next_state = matrix(0, nrow = 1, ncol = n.state)
  
  # prebreeders established in LR colony
  if (settlement_colony == 1){
  next_state[1,(1*n.colony+1)] = phiA * kappa[1] # prebreeders become breeders
  next_state[1,(2*n.colony+1)] = phiA * (1-kappa[1]) # prebreeders stay prebreeders
  }
  # prebreeders established in satellite colonies
  else {
    next_state[1,(1*n.colony+settlement_colony)] = phiA * kappa[2] # prebreeders become breeders
    next_state[1,(2*n.colony+settlement_colony)] = phiA * (1-kappa[2]) # prebreeders stay prebreeders
  }
  next_state[1,n.state] = 1 - phiA # become dead
  return(next_state)
}


TPM = matrix(NA , nrow= n.state, ncol=n.state)

# Building TPM
for (col in 1:n.colony){
  TPM[col,] = TMP_nestling(birth_colony = col, n.colony = n.colony, phi1 = phi1, eta = eta)
  TPM[(n.colony + col),] = TMP_breeder(former_colony = col, n.colony = n.colony, phiA = phiA, nu = nu)
  TPM[(2* n.colony + col),] = TMP_prebreeder(settlement_colony = col, n.colony = n.colony, phiA = phiA, kappa = kappa)
}

# Last line of TPM : Dead stay dead
TPM[n.state, ] = 0
TPM[n.state, n.state] = 1

# Construct the observation probability matrix (OPM), above called THETA
# True state is in rows, observed state is in columns

OPM = matrix(0 , nrow= n.state, ncol=n.state)
for (col in 1:n.colony){
  OPM[n.colony + col, n.colony + col] = diag(p[1,col],1)
}
OPM[,ncol(OPM)] = c(1-rowSums(OPM))

# State or ecological process
# Simulate true system state
z = array(NA, dim=c(nind, nyears)) # Empty alive/dead matrix
# Initial conditions: all individuals alive at f(i)
initial.state = rep(1:5, each = nind/5)

for (i in 1:nind){
  z[i,f[i]] = initial.state[i]
}

# Propagate alive/dead process forwards via transition rule (=TPM=OMEGA)
for (i in 1:nind){
  for (t in (f[i]+1):nyears){
    departure.state = z[i,t-1]
    arrival.state = which(rmultinom(1,1, TPM[departure.state,])==1)
    z[i,t] = arrival.state
  } #t
} #i


# Observation process: simulate observations using observation matrix OPM (=THETA)
y = array(16, dim=c(nind, nyears))
for (i in 1:nind){
  y[i,f[i]] = z[i,f[i]]
  for (t in (f[i]+1):nyears){
    true.state = z[i,t-1]
    observed.state = which(rmultinom(1,1, OPM[true.state,])==1)
    y[i,t] = observed.state
  } #t
} #i
