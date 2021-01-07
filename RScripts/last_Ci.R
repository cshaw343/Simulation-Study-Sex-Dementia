last_Ci <- function(obs){
  last_Ci <- vector(length = ncol(obs))
  
  for(i in 1:length(last_Ci)){
    if(obs["study_death", i] == 1){
      survtime <- obs["survtime", i]
      last_Ci[i] <- 
        obs["a0", i] + obs["a1", i]*survtime + obs["a2", i]*survtime^2
    } else{
      last_Ci[i] <- NA
    }
  }
  
  return(last_Ci)
}
