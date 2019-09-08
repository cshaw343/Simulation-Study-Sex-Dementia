compute_Cij <- function(obs){
  Cij <- matrix(nrow = nrow(obs), ncol = length(visit_times))
  for(i in 1:ncol(Cij)){
    centered_age <- visit_times[i]
    Cij[, i] <- obs[, "a0"] + obs[, "a1"]*centered_age + 
      obs[, "a2"]*centered_age^2
  }
  
  return(Cij)
}
