#---- Package loading, options, seed ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "MASS")

#Using standard notation (as opposed to scientific), rounded to three 
#decimal places
options(scipen = 999)
options(digits = 3)

set.seed(10789)

#---- Specify the parameter file ----
source("sex-demensia_sim_parA.R")

#---- Generating variable names for each assessment timepoint ----
#Age labels at each assessment timepoint
age_varnames <- c("id", "age0", vector(length = num_tests))
for(i in 1:num_tests){
  age_varnames[i + 2] = paste("age", i, sep = "")
}

#Centered age labels at each assessment timepoint
agec_varnames <- c("id", "age0_c50", vector(length = num_tests))
for(i in 1:num_tests){
  agec_varnames[i + 2] = paste(age_varnames[i + 2], "_c50", sep = "")
}

#Alpha (for autoregressive noise) labels at each assessment timepoint
eps_varnames <- c("id", "eps0", vector(length = num_tests))
for(i in 1:num_tests){
  eps_varnames[i + 2] = paste("eps", i, sep = "")
}

#---- Generating assessment timepoint data ----
visit_times <- seq(from = 0, to = int_time*num_tests, by = int_time)

#---- Model for Cognitive Function ----
Cij <- function(obs, t){
  knots = c(0, 20, 35)
  if(t >= knots[1] & t < knots[2]){
    return(b00 + z0i + boi*sex + b02*age)
  }
}

#---- Generate Covariance Matrix for random slope and intercept terms ----
slope_int_cov <- matrix(c(var0, cov, cov, var1), nrow = 2, byrow = TRUE)

#---- The simulation function ----
sex_dem_sim <- function(){
  #---- Generating IDs, sex, U ----
  obs <- tibble("id" = seq(from = 1, to = num_obs, by = 1), 
                "sex" = rbinom(num_obs, size = 1, prob = psex), 
                "U" = rnorm(num_obs, mean = 0, sd = 1))
  
  #---- Generating age data ----
  #Creating ages at each timepoint j
  ages <- as_tibble(matrix(NA, nrow = num_obs, ncol = length(age_varnames)))
  for(i in 1:length(age_varnames)){
    if(i == 1){
      ages[, i] = seq(from = 1, to = num_obs, by = 1)
    } else if(i == 2){
      ages[, i] = age0
    } else ages[, i] = ages[, (i-1)] + int_time
  }
  colnames(ages) <- age_varnames
  
  #---- Generating centered age data ----
  #Creating centered ages at each timepoint j
  c_ages <- as_tibble(ages - mean(age0)) %>% 
    mutate("id" = seq(from = 1, to = num_obs, by = 1))
  colnames(c_ages) <- agec_varnames
  
  #---- Generating "true" cognitive function Cij ----
  #Cij = b00 + z0i + bo1*sexi + b02*age0i + b03*Ui + (b10 + z1i + b11*sexi + 
  #b12*age0i + b13*Ui)*timej + epsilonij
  
    #---- Generating random terms for slope and intercept ----
    #Generate random terms for each individual
    slope_int_noise <- as_tibble(mvrnorm(n = num_obs, mu = rep(0, 2), 
                               Sigma = slope_int_cov))
    colnames(slope_int_noise) <- c("z0i", "z1i")
  
    #---- Generating noise term (unexplained variance in Cij) for each visit ----
    sd_eps <- sqrt(var3)
    eps <- replicate(num_tests, rnorm(n = num_obs, mean = 0, sd = sd_eps))
    colnames(eps) <- eps_varnames
  
    #---- Creating full matrix of data for each individual ----
    obs <- cbind(obs, cbind(slope_int_noise, eps))
    
    #---- Calculating Cij for each individual ----
    cog_func <- as_tibble(matrix(NA, nrow = num_obs, 
                                 ncol = length(visit_times)))
    for(i in 1:length(visit_times)){
      t = visit_times[i]
      cog_func[, i] = apply(noise, 1, Cij)
    }
    
  
  
  
}

