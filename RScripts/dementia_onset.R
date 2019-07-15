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

      #Regular linear interpolation
       slope <- paste("cij_slope", wave - 1, "-", wave, sep = "")
       int <- paste("Ci", wave - 1, sep = "")
       m = as.double(obs_matrix[slope, i])
       b = as.double(obs_matrix[int, i])
      
       time = optimize(event_est, M = m, B = b, dem_cut = dem_cut,
                       interval = c(0, 5))$minimum
       timetodem[i] = (wave - 1)*int_time + time
      
      
      # #Linear interpolation of dementia cutoffs
      # interp_dem_cuts = seq(dem_cuts[wave], dem_cuts[wave + 1], length = 6)
      # for(j in 1:(length(interp_dem_cuts) - 1)){
      #   test_time = optimize(event_est, M = m, B = b,
      #                        dem_cut = interp_dem_cuts[j + 1],
      #                        interval = c(0, 5))$minimum
      #   
      #   if(test_time > (j - 1) & test_time <= j){
      #     timetodem[i] = (wave - 1)*int_time + test_time
      #     break
      #   }
      # }
      # 
      # #Probabilistic draw of dementia cutoffs
      # for(j in 1:length(interp_dem_cuts)){
      #   if(rbernoulli(1, p = 0.75)){
      #     dem_cut = interp_dem_cuts[j + 1]
      #   } else {
      #     dem_cut = interp_dem_cuts[j]
      #   }
      #   test_time = optimize(event_est, M = m, B = b,
      #                        dem_cut = dem_cut,
      #                        interval = c(0, 5))$minimum
      #   if(test_time > (j - 1) & test_time <= j){
      #     timetodem[i] = (wave - 1)*int_time + test_time
      #     break
      #   }
      # }
      
      # #Probabilistic draw of bins and timetodem
      # timetodem[i] = (wave - 1)*int_time + 
      #   #Randomly draw the bin
      #   #sample(seq(0, 4, by = 1), size = 1, 
      #          #prob = c(0.10, 0.10, 0.20, 0.30, 0.30)) + 
      #   #Randomly draw the time within the bin  
      #   5*rbeta(1, shape1 = 5, shape2 = 1)
      #   #runif(n = 1, min = 0, max = 1)
      
      # #Uniform assignment
      # timetodem[i] = (wave - 1)*int_time + runif(n = 1, min = 0, max = 5)
      
    }
  }
  return(timetodem)
}
