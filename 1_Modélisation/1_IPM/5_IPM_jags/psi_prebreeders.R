
for (prev.state in prebr.states){
  # Prebreeders may become breeders
  for (next_state in breed.states){
    if (state2col[next_state] == state2col[prev.state]){
      psi[prev.state,next_state] <- phi[2] * kappa[state2kappa[prev.state]]
    }
    else {
      psi[prev.state,next_state] <- 0
    }
  }
  # prebreeders may stay prebreeders
  for (next_state in prebr.states){
    if (state2col[next_state] == state2col[prev.state]){
      psi[prev.state,next_state] <- phi[2] * (1 - kappa[state2kappa[prev.state]])
    }
    else {
      psi[prev.state,next_state] <- 0
    }
  }
  # prebreeders cannot become nestling
  for (next_state in nest.states){
    psi[prev.state,next_state] <- 0
  }
}

