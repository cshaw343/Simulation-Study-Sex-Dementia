#---- Model for Survival Time ----
survival <- function(obs, lambda){
  #Calculate survival times for each interval
  Sij <- vector(length = (length(visit_times) - 1))
  for(j in 1:length(Sij)){
    test_num = j - 1
    US <- paste("U", test_num, test_num + 1, sep = "")
    agec <- paste("age", test_num, "_c50", sep = "")
    slope <- paste("slope", test_num, test_num + 1, sep = "")
    C <- paste("Ci", test_num, sep = "")
    Sij[j] = -log(obs[US])/
      (lambda[j]*exp(g1[j]*obs["sex"] + g2*obs[agec] + g3*obs["U"] + 
                       g4*obs["sex"]*obs[agec] + g5*obs[slope] + 
                       g6*obs[C]))
  }
  return("Sij" = Sij)
}