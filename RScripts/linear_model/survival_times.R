#---- Model for Survival Time ----
survival <- function(obs_matrix){
  #Calculate survival times for each interval
  Sij <- matrix(ncol = ncol(obs_matrix), nrow = (length(visit_times) - 1))
  for(i in 1:ncol(obs_matrix)){
    survtimes <- matrix(NA, nrow = (length(visit_times) - 1), ncol = 1)
    for(j in 1:length(survtimes)){
      r_name <- variable_names$rij_varnames[j]
      agec_name <- variable_names$agec_varnames[j]
      cij_slope_name <- variable_names$cij_slopeij_varnames[j]
      C_name <- variable_names$Cij_varnames[j]
      dem_name <- variable_names$dem_varnames[j]
      survtime = -log(obs_matrix[r_name, i])/
        (lambda[j]*exp(g1[j]*obs_matrix["female", i] + 
                         g2*obs_matrix["U", i] + 
                         g3*obs_matrix["female", i]*obs_matrix["U", i] 
                         # + g4*obs_matrix[cij_slope_name, i] + 
                         #g5*obs_matrix[C_name, i] + 
                         #g6[j]*obs_matrix[dem_name, i]
                         ))
      
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

