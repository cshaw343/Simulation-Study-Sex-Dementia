#---- Model for Survival Time ----
survival <- function(obs_matrix){
  #Calculate survival times for each interval
  Sij <- matrix(nrow = nrow(obs_matrix), ncol = (length(visit_times) - 1))
  for(i in 1:nrow(obs_matrix)){
    survtimes <- matrix(NA, ncol = (length(visit_times) - 1), nrow = 1)
    for(j in 1:length(survtimes)){
      r_name <- variable_names$rij_varnames[j]
      agec_name <- variable_names$agec_varnames[j]
      cij_slope_name <- variable_names$cij_slopeij_varnames[j]
      C_name <- variable_names$Cij_varnames[j]
      survtime = -log(obs_matrix[i, r_name])/
        (lambda[j]*exp(g1[j]*obs_matrix[i, "sex"] + 
                         g2*obs_matrix[i, agec_name] + g3*obs_matrix[i, "U"] + 
                         g4*obs_matrix[i, "sex"]*obs_matrix[i, agec_name] + 
                         g5*obs_matrix[i, cij_slope_name] + 
                         g6*obs_matrix[i, C_name]))
      
      survtimes[j] <- as.numeric(survtime)
      if(survtime < 5){
        break
      }
    }
    Sij[i, ] <- survtimes
  }
  Sij[Sij > 5] <- 4.999999
  return(list("Sij" = Sij, "survtimes" = rowSums(Sij, na.rm = TRUE)))
}

