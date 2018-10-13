#---- Model for Survival Time ----
survival <- function(obs, lambda){
  #Calculate survival times for each interval
  Sij <- vector(length = (length(visit_times) - 1))
  for(j in 1:length(Sij)){
    r <- variable_names$rij_varnames[j]
    agec <- variable_names$agec_varnames[j]
    cij_slope <- variable_names$cij_slopeij_varnames[j]
    C <- variable_names$Cij_varnames[j]
    Sij[j] = -log(obs[r])/
      (lambda[j]*exp(g1[j]*obs["sex"] + g2*obs[agec] + g3*obs["U"] + 
                       g4*obs["sex"]*obs[agec] + g5*obs[cij_slope] + 
                       g6*obs[C]))
  }
  return("Sij" = Sij)
}

