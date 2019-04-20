#---- Function for computing dementia onset ----
dem_onset <- function(obs_matrix, dem_cuts){
  #Function we are trying to find the root of
  event_est <- function(t, M, B, dem_cut){
    cij = B + M*t
    return(abs(cij - dem_cut))
  }
    
  timetodem <- vector(length = nrow(obs_matrix))
  for(i in 1:nrow(obs_matrix)){
    if(is.na(obs_matrix[i, "dem_wave"])){
      timetodem[i] = as.double(obs_matrix[i, "survtime"])
    } else {
      wave = as.double(obs_matrix[i, "dem_wave"])
      slope <- paste("cij_slope", wave - 1, "-", wave, sep = "")
      int <- paste("Ci", wave - 1, sep = "")
      m = as.double(obs_matrix[i, slope])
      b = as.double(obs_matrix[i, int])
      interp_dem_cuts = seq(dem_cuts[wave], dem_cuts[wave + 1], length = 5)
      
      for(j in 1:length(interp_dem_cuts)){
        test_time = optimize(event_est, M = m, B = b, 
                             dem_cut = interp_dem_cuts[j], 
                             interval = c(0, 5))$minimum
        
        if(test_time > (j - 1) & test_time <= j){
          timetodem[i] = (wave - 1)*int_time + test_time
          break
        }
      }
    }
  }
  return(timetodem)
}
