survival_censor <- function(obs){
  max_Cij = length(variable_names$Cij_varnames)
  max_dem = length(variable_names$dem_varnames)
  
  for(i in 1:ncol(obs)){
    if(obs["study_death", i] == 1){
      last_wave <- floor(obs["survtime", i]/5)
      obs[variable_names$Cij_varnames[(last_wave + 2):max_Cij], i] <- NA
      obs[variable_names$dem_varnames[(last_wave + 2):max_dem], i] <- NA
    }
  }
  
  return(obs)
}
