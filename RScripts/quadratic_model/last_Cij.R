last_Cij <- function(obs){
  last_Cij <- vector(length = ncol(obs))
  
  for(i in 1:length(last_Cij)){
    if(obs["study_death", i] == 1){
      survtime <- obs["survtime", i]
      last_Cij[i] <- 
        obs["a0", i] + obs["a1", i]*survtime + obs["a2", i]*survtime^2
    } 
  }
  
  return(last_Cij)
}
