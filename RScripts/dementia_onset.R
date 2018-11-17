#---- Function for computing dementia onset ----
dem_onset <- function(obs, dem_cuts){
  #Function we are trying to find the root of
  event_est <- function(t, M, B, dem_cut){
    cij = B + M*t
    return(abs(cij - dem_cut))
  }
  timetodem <- vector(length = nrow(obs))
  for(i in 1:nrow(obs)){
    data <- obs[i, ]
    if(is.na(data["dem_wave"])){
      timetodem[i] = as.double(data["survtime"])
    } else if(data["dem_wave"] == 0){
      timetodem[i] = 0
    } else {
      wave = as.double(data["dem_wave"])
      slope <- paste("cij_slope", wave - 1, "-", wave, sep = "")
      int <- paste("Ci", wave - 1, sep = "")
      m = as.double(data[slope])
      b = as.double(data[int])
      timetodem[i] = (wave - 1)*int_time + 
        optimize(event_est, M = m, 
                 B = b, dem_cut = dem_cuts[(wave + 1)], 
                 interval = c(0, 5))$minimum
    }
  }
  return(timetodem)
}
