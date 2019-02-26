#---- Model for Survival Time ----
survival <- function(obs_matrix, lambda){
  #Calculate survival times for each interval
  Sij <- matrix(nrow = nrow(obs_matrix), ncol = (length(visit_times) - 1))
  for(i in 1:nrow(obs_matrix)){
    survtimes <- matrix(NA, ncol = (length(visit_times) - 1), nrow = 1)
    for(j in 1:length(survtimes)){
      r <- variable_names$rij_varnames[j]
      agec <- variable_names$agec_varnames[j]
      cij_slope <- variable_names$cij_slopeij_varnames[j]
      C <- variable_names$Cij_varnames[j]
      survtime = -log(obs_matrix[i, r])/
        (lambda[j]*exp(g1[j]*obs_matrix[i, "sex"] + g2*obs_matrix[i, agec] + 
                         g3*obs_matrix[i, "U"] + 
                         g4*obs_matrix[i, "sex"]*obs_matrix[i, agec] + 
                         g5*obs_matrix[i, cij_slope] + g6*obs_matrix[i, C]))
      survtimes[j] <- as.numeric(survtime)
      if(survtime < 5){
        break
      }
    }
    Sij[i, ] <- survtimes
  }
  return("Sij" = Sij)
}

