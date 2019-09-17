gen_random_effects <- function(obs, cij_var0, cij_var1, cij_cov, cij_var3){
  #Covariance matrices for random slope and intercept terms
  cij_slope_int_cov <- lapply(1:(num_tests + 1), 
                              function(x) matrix(NA, nrow = 2, ncol = 2))
  for(i in 1:(num_tests + 1)){
    cij_slope_int_cov[[i]] <- matrix(c(cij_var0, cij_cov, 
                                       cij_cov, cij_var1[i]), 
                                     nrow = 2, byrow = TRUE)
  }
  
  #Generate random terms for each individual
  for(i in 1:(num_tests + 1)){
    noise <- mvrnorm(n = num_obs, mu = rep(0, 2), 
                     Sigma = cij_slope_int_cov[[i]]) 
    if(i == 1){
      obs[, c("z0_i", paste0("z1_", (i - 1), "i"))] <- noise
    } else{
      obs[, paste0("z1_", (i - 1), "i")] <- noise[, 2]
    }
  }
  
  #Generating noise term (unexplained variance in Cij) for each visit
  #Creating AR(1) correlation matrix
  num_visits = num_tests + 1
  powers <- abs(outer(1:(num_visits), 1:(num_visits), "-")) #Exponents for autoregressive terms
  corr <- sqrt(cij_var3)*(cij_r1^powers)                    #Correlation matrix
  S <- diag(rep(sqrt(cij_var3)), nrow(corr))                #Diagonal matrix of SDs
  cij_cov_mat <- S%*%corr%*%S                               #Covariance matrix
  
  #Generating noise terms
  obs[, variable_names$eps_varnames] <- 
    mvrnorm(n = num_obs, mu = rep(0, num_visits), Sigma = cij_cov_mat)
  
  return(obs)
}

