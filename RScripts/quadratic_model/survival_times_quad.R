#---- Model for Survival Time ----
survival <- function(obs_matrix){
  #Calculate survival times for each interval
  Sij <- matrix(ncol = ncol(obs_matrix), nrow = (length(visit_times) - 1))
  for(i in 1:ncol(obs_matrix)){
    survtimes <- matrix(NA, nrow = (length(visit_times) - 1), ncol = 1)
    for(j in 1:length(survtimes)){
      r_name <- variable_names$rij_varnames[j]
      
      survtime = -log(obs_matrix[r_name, i])/
        (lambda[j]*exp(g1[j]*obs_matrix["female", i] + 
                         g2*obs_matrix["U", i] + 
                         g3*obs_matrix["female", i]*obs_matrix["U", i]))
      
      survtimes[j] <- as.numeric(survtime)
      if(survtime < 5){
        break
      }
    }
    Sij[, i] <- survtimes
  }
  Sij[Sij > 5] <- 5
  return(Sij)
}

