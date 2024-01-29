
# Ce qu'il faut sauvegarder à chaque fois qu'on veut faire tourner un modèle --------------------------------
for (chain in 1:3){
  name <- paste0("NimbleInFile/Set", set, "Rep", rep, "Chain",chain,".RData")
  save(modelCode, nimConstants, nimData, nimInits, file = name)
}



# Dans le script qu'on lance pour faire tourner un modèle on met cette partie là -----------------------------
path<-"NimbleInFile/" # définir le chemin où il y a les chaines 

path3<-"NimbleOutFile/" # définir le chemin pour sauvegarder les chaines 

pausesecs <- sample(seq(5,10, by=0.1), 1) # pause aléatoire de quelques secondes pour éviter que la même chaîne soit lancé plusieur fois 
Sys.sleep(pausesecs)


## ----- PICK A FILE AT RANDOM -----
file.list <- list.files(path)
set <- sample(file.list, 1)

## ----- LOAD THE FILE -----
load(paste(path, set, sep="")) # charger la chaine 
file.remove(paste(path, set, sep="")) # la supprimer du dossier d'entrée pour éviter quelle soit reprise

## ----- CREATE NIMBLE OBJETS -----
#nimbleOptions(clearNimbleFunctionsAfterCompiling = TRUE)
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      inits = nimInits,
                      data = nimData,
                      check = FALSE,
                      calculate = FALSE)

cmodel <- compileNimble(model)

calc <- cmodel$calculate() 
print(calc)


if(calc =="-Inf"){
  #file.remove(paste(path2, set, sep=""))                                    
  stop("There is an INF")
}


MCMCconf <- configureMCMC(model = model, 
                          monitors = c("p0","sigma","psi","gamma", "N","alpha","sd","gamma","phi"),
                          thin = 1,
                          monitors2 = c("s","z"),
                          thin2 = 500)

MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE,showCompilerOutput = TRUE)

### ====    4.2. RUN THE MCMC ====
MCMCRuntime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC,
                                                      nburnin = 1000,
                                                      niter = 50000,
                                                      nchains = 1,
                                                      samplesAsCodaMCMC = TRUE))

myNimbleOutput[["samples2"]] <- round(myNimbleOutput[["samples2"]],6) # make object smaller to store

outname <- paste(path3, "NimbleOutFOR", set, sep="")
save(myNimbleOutput, MCMCRuntime, file = outname)

# Combiner toutes les chaines a posteriori -------------------------------------
setwd("/Users/maeliskervellec/Documents/Projets these/4 - OPSCR/NonEuclideanOPSCR/Simulations/sims4/NimbleOutFile")
Output <- Output_samples2 <- Output1_1 <- Output1_2 <- Output2_1 <- Output2_2 <- list()
OutputNA <- tibble(name = NA )
for(set in 1:4){
  for(rep in 1:5){
    for(chain in 1:3){
      name <- paste0("NimbleOutFORSet",set,"Rep",rep,"Chain",chain,".RData")
      if(any(grepl(name,list.files()))){
        load(name)
        Output2_1[[paste0("Chain",chain)]] <- myNimbleOutput$samples
        Output2_2[[paste0("Chain",chain)]] <- myNimbleOutput$samples2
        # print(name)
        # print(MCMCRuntime)
      }
      else {
        OutputNA <- OutputNA %>% rbind(name) # garder le nom des chaines qui n'ont pas tournée 
      }
      
    }
    Output1_1[[paste0("Rep",rep)]] <- Output2_1
    Output1_2[[paste0("Rep",rep)]] <- Output2_2
  }
  Output[[paste0("Set",set)]] <- Output1_1
  Output_samples2[[paste0("Set",set)]] <- Output1_2
}


# visualiser les chaines 
MCMCvis::MCMCtrace(object = Output[["Set4"]][[Rep]],
                   pdf = FALSE, # no export to PDF
                   ind = TRUE,
                   Rhat = TRUE, # add Rhat
                   n.eff = TRUE,
                   params = c("alpha","sd"),
                   iter = 50000) # only plot the 9000 last iteration 


MCMCvis::MCMCsummary(Output[["Set1"]][[1]]) 


# Commande sbatch pour lancer plusieurs fois le même jobs => $ sbatch -a 1-60%10 Sims.sbatch (lance les 60 jobs en utilisant 10 coeurs)