compute_Ci <- function(obs){
  Ci <- matrix(nrow = nrow(obs), ncol = length(visit_times))
  for(i in 1:ncol(Ci)){
    centered_age <- visit_times[i]
    Ci[, i] <- obs[, "a0"] + obs[, "a1"]*centered_age + 
      obs[, "a2"]*centered_age^2
  }
  
  return(Ci)
}
