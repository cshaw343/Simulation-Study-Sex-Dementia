compare_survtime_timetodem <- function(obs_matrix){
 obs_matrix[(obs_matrix$survtime < obs_matrix$timetodem), "dem"] <- 0
 obs_matrix[(obs_matrix$survtime < obs_matrix$timetodem), "dem_wave"] <- NA
 obs_matrix[(obs_matrix$survtime < obs_matrix$timetodem), 
            variable_names$dem_varnames] <- 
   obs_matrix[(obs_matrix$survtime < obs_matrix$timetodem), 
              variable_names$dem_varnames]*0
 
 return(obs_matrix)
}
