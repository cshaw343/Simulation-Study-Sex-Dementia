survival_censor <- function(obs){
  max_Ci = length(variable_names$Ci_varnames)
  
  for(i in 1:ncol(obs)){
    if(obs["study_death", i] == 1){
      last_wave <- floor(obs["survtime", i]/5)
      obs[variable_names$Ci_varnames[(last_wave + 2):max_Ci], i] <- NA
    }
  }
  
  return(obs)
}
