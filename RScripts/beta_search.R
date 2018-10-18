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
  risk <- sum(dem_diag)/length(dem_diag)
  
  return(abs(risk - risk_match))
}

#---- Function that looks for the betas ----
find_betas <- function(obs, dem_risk){
  beta_results <- matrix(nrow = length(dem_risk), ncol = 6)
  start = 4 #represent starting age where we have data (age 70)
  for(i in 1:nrow(beta_results)){
    sex_i <- obs$sex
    Cij <- obs[, variable_names$Cij_varnames[start + 1]] 
    Fij <- obs[, variable_names$Fij_varnames[start + 1]] 
    
    sex_i_nomiss <- sex_i[-which(is.na(Cij))]
    Cij_nomiss <- Cij[-which(is.na(Cij))]
    Fij_nomiss <- Fij[-which(is.na(Cij))]
    
    many_opts <- replicate(10, 
                           {optim(par = rep(-0.05, 5), 
                                  fn = dem_diag, 
                                  sex_i = sex_i_nomiss, Cij = Cij_nomiss, 
                                  Fij = Fij_nomiss, 
                                  risk_match = dem_risk[[i]], 
                                  upper = rep(0, 5), 
                                  lower = rep(-1, 5))})
    
    mean_diff <- many_opts["value", ] %>% unlist() %>% mean()
    mean_betas <- matrix(unlist(many_opts["par", ]), ncol = 5, byrow = TRUE) %>%
      colMeans()
    
    beta_results[i, ] <- c(mean_betas, mean_diff)
    
    start = start + 1
  }
}