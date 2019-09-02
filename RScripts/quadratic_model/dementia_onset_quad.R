#---- Function for computing dementia onset ----
dem_onset <- function(obs_matrix, dem_cut){
  
  timetodem <- vector(length = ncol(obs_matrix))
  for(i in 1:ncol(obs_matrix)){
    if(is.na(obs_matrix["dem_wave", i])){
      timetodem[i] = as.double(obs_matrix["survtime", i])
    } else {
      wave = obs_matrix["dem_wave", i]
      
      #Function we are trying to find the root of
      quad_function <- function(x){
        abs(obs["a0", i] + obs["a1", i]*x + obs["a2", i]*x^2 - dem_cut)
      }
      
      obs["timetodem", i] <- 
        optimize(quad_function, lower = 5*(wave - 1), upper = 5*(wave))$minimum
    }
  }
  return(timetodem)
}

