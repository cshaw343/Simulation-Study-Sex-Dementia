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

#---- Function that looks for the betas ----
find_betas <- function(obs, dem_risk){
  beta_mat <- matrix(nrow = length(dem_risk), ncol = 5)
  #Doing this for the first case-- will use loop later
  sex_i <- obs$sex
  Cij <- obs$Ci4 
  Fij <- obs$Fi4
  
  #Don't consider the cases who have died by age 65
  sex_i_nomiss <- sex_i[-which(is.na(Cij))]
  Cij_nomiss <- Cij[-which(is.na(Cij))]
  Fij_nomiss <- Fij[-which(is.na(Cij))]
  #Test with values for B first
  B_0 <- 0
  B_1 <- -0.5
  B_2 <- -0.5
  B_3 <- -0.5
  B_4 <- -0.25
  
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
  
  optim(par = rep(-0.05, 5), 
        fn = dem_diag, 
        sex_i = sex_i_nomiss, Cij = Cij_nomiss, Fij = Fij_nomiss, 
        risk_match = dem_risk[[1]],
        upper = rep(0, 5), 
        lower = rep(-0.5, 5))
  
}