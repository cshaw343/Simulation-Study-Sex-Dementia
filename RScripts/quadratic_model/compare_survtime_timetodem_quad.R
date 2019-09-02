compare_survtime_timetodem <- function(obs_matrix){
  indices = obs_matrix["survtime", ] < obs_matrix["timetodem", ]
  if(sum(indices) > 0){
    obs_matrix["dem", indices] <- 0
    obs_matrix["dem_wave", indices] <- NA
    obs_matrix[variable_names$dem_varnames, indices] <- 
      obs_matrix[variable_names$dem_varnames, indices]*0
  }
  return(obs_matrix)
}
