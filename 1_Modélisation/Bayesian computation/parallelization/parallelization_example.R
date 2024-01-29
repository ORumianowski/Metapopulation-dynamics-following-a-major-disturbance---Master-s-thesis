

library(parallel)


this_cluster <- makeCluster(4)

set.seed(10120)

# Simulate some data
myData <- rgamma(1000, shape = 0.4, rate = 0.8)

library(nimble)

myCode <- nimbleCode({
  a ~ dunif(0, 100)
  b ~ dnorm(0, 100)
  
  for (i in 1:length_y) {
    y[i] ~ dgamma(shape = a, rate = b)
  }
})

# Create a function with all the needed code
run_MCMC_allcode <- function(seed, data, code, useWAIC = TRUE) {
  library(nimble)
  
  myModel <- nimbleModel(code = code,
                         data = list(y = data),
                         constants = list(length_y = 1000),
                         inits = list(a = 0.5, b = 0.5))
  
  CmyModel <- compileNimble(myModel)
  
  ## One may also wish to add additional monitors
  
  myMCMC <- buildMCMC(CmyModel)
  CmyMCMC <- compileNimble(myMCMC)
  
  results <- runMCMC(CmyMCMC, niter = 10000, setSeed = seed)
  
  return(results)
}

chain_output <- parLapply(cl = this_cluster, X = 1:4, 
                          fun = run_MCMC_allcode, 
                          data = myData, code = myCode,
                          useWAIC = useWAIC)

# It's good practice to close the cluster when you're done with it.
stopCluster(this_cluster)

par(mfrow = c(2,2))
for (i in 1:4) {
  this_output <- chain_output[[i]]
  plot(this_output[,"b"], type = "l", ylab = 'b')
}
