#---- Function for computing dementia onset ----
dem_onset <- function(obs_matrix, dem_cut){
  #Function we are trying to find the root of
  event_est <- function(t, M, B, dem_cut){
    cij = B + M*t
    return(abs(cij - dem_cut))
  }
    
  timetodem <- vector(length = ncol(obs_matrix))
  for(i in 1:ncol(obs_matrix)){
    if(is.na(obs_matrix["dem_wave", i])){
      timetodem[i] = as.double(obs_matrix["survtime", i])
    } else {
      wave = as.double(obs_matrix["dem_wave", i])

      #Uniform assignment
      timetodem[i] = (wave - 1)*int_time + runif(n = 1, min = 0, max = 5)
      
    }
  }
  return(timetodem)
}
