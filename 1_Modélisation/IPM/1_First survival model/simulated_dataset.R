
# Simulated Dataset  ------------------------------------------------------

# How to use
#res = simul_data()
#dataset = res$CH.A


simul_data = function(total.marked.j = 55000, n.occasions = (2023 - 1986)){
 

# Dataset properties : define parameter values

#total.marked.j <- 55000 # 55000   # Total number of marked juveniles

n.occasions <-  (2023 - 1986) # Number of capture occasions  # p32 thèse Perron
marked.j <- rep( trunc(total.marked.j / n.occasions-1), n.occasions-1) # Annual number of newly marked juveniles
marked.a <- rep(1, n.occasions-1) # Annual number of newly marked adults
phi.juv <- 0.243  # Juvenile annual survival # p42 thèse Perron
phi.ad <- 0.844 # Adult annual survival # p42 thèse Perron

prob.recap.mean = 0.2
#p <- rep(prob.recap.mean, n.occasions-1) # Recapture # p137 Perron

# Alternative pour la probabilité de recapture
# Paramètres de la loi beta
esperance <- prob.recap.mean
ecart_type <- 0.1

# Calcul des paramètres alpha et beta de la loi beta
alpha <- (esperance^2 * (1 - esperance)) / ecart_type^2 - esperance
beta <- alpha * (1 / esperance - 1)

# # Recapture variable
p <-rbeta(n.occasions-1, alpha, beta)

#mean(p)

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
    if (mark.occ[i]==n.occasions) next 
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

return(list(CH.J = CH.J, CH.A = CH.A))
}

