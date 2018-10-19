#*******************************************************************************
# This is the function that performs the search for the beta values in the 
# logistic regression of dementia diagnosis on Cij, Fij, and sex
#
# This function is attempting to find betas that will create dementia incidence
# rates that match the Supplemental Table 2 in Inequalities in dementia 
# incidence between six racial and ethnic groups over 14 years 
# (Mayeda et al, 2016)
#*******************************************************************************

#---- Source Files ----
source("RScripts/dementia_incidence2000-2013.R")

#---- Function that we are trying to optimize ----
dem_diag <- function(BETAS, 
                     sex_i, Cij, Fij, 
                     risk_match){
  
  model <- BETAS[1] + BETAS[2]*Cij + BETAS[3]*Fij + BETAS[4]*sex_i + 
    BETAS[5]*Cij*Fij
  pi_hat <- exp(model)/(1 + exp(model))
  dem_diag <- rbernoulli(1, p = pi_hat)*1
  risk_percent <- mean(dem_diag)*100
  
  return(abs(risk_percent - risk_match))
}

#---- Function that looks for the betas ----
find_betas <- function(obs, dem_risk){
  beta_results <- matrix(nrow = length(dem_risk), ncol = 6)
  start = 4 #represent starting age where we have data (age 70)
  data <- obs %>% 
    dplyr::select("sex", 
                  variable_names$Cij_varnames[(start + 1):nrow(variable_names)], 
                  variable_names$Fij_varnames[(start + 1):
                                                nrow(variable_names)]) %>% 
    mutate("dem_case" = rep(0, nrow(obs)))
  
  for(i in 1:nrow(beta_results)){
    #target_measurement <- variable_names$Cij_varnames[start + 1] %>% noquote()
    data %<>% filter(dem_case == 0)
    
    sex_i <- data$sex
    Cij <- data[, variable_names$Cij_varnames[start + 1]] 
    Fij <- data[, variable_names$Fij_varnames[start + 1]] 

    # if(i == 1)
    {
      many_opts <- replicate(10, 
                             {optim(par = rep(-13, 5), 
                                    fn = dem_diag, 
                                    sex_i = sex_i, Cij = Cij, Fij = Fij, 
                                    risk_match = dem_risk[[i]], 
                                    upper = rep(0, 5), 
                                    lower = rep(-15, 5))})
      
      mean_diff <- many_opts["value", ] %>% unlist() %>% mean()
      mean_betas <- matrix(unlist(many_opts["par", ]), ncol = 5, byrow = TRUE) %>%
        colMeans()
    # } else {
    #   many_opts <- replicate(10, 
    #                          {optim(par = rep(-4.5, 5), 
    #                                 fn = dem_diag, 
    #                                 sex_i = sex_i_nomiss, Cij = Cij_nomiss, 
    #                                 Fij = Fij_nomiss, 
    #                                 risk_match = dem_risk[[i]], 
    #                                 upper = rep(-4, 5), 
    #                                 lower = rep(-5, 5))})
    }
    
    beta_results[i, ] <- c(mean_betas, mean_diff)
    
    mean_model <- mean_betas[1] + mean_betas[2]*Cij + mean_betas[3]*Fij + 
      mean_betas[4]*sex_i + mean_betas[5]*Cij*Fij
    
    pi_hats <- exp(mean_model)/(1 + exp(mean_model))
    data %<>% mutate("dem_case" = rbernoulli(1, p = pi_hats)*1)
    
    start = start + 1
  }
  return(beta_results)
}

fingers_crossed <- find_betas(results_mat, dem_risk = dem_risk)
