#***************************************************************
# Performs a search for an appropriate dementia cutpoint 
# Reference data found in dem_calcs(Gross, et.al).R
#***************************************************************

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "magrittr")

#Suppress warnings
options(warn = -1)

#---- Specify source files ----
source("sex-dementia_sim_parA.R")
source("dementia_incidence2000-2013.R")
source("sex-dementia_sim_script.R")

#---- Dementia search function ----
find_demcut <- function(dem_table){
  #Function we are trying to optimize
  dem_rates <- function(parameters, J, dem_val){
    CslopeB <- parameters[1]
    CslopeC <- parameters[2]
    dem_cut <- parameters[3]
  #---- Model for Cognitive Function ----
  opt_cog_func <- function(CslopeB, CslopeC, obs){
    knots = c(0, 20, 35)
    Cij <- vector(length = length(visit_times))
    for(j in 1:length(Cij)){
      t = visit_times[j]
      test_num = j - 1
      eps <- paste("eps", test_num, sep = "")
      Cij[j] = b00 + obs[, "z0i"] + b01*obs[, "sex"] + b02*obs[, "age0_c50"] + 
        b03*obs[, "U"] + obs[, eps] + 
        (b10a - CslopeB)*knots[2]*(t >= knots[2]) + 
        (CslopeB - CslopeC)*knots[3]*(t >= knots[3]) + 
        (obs[, "z1i"] + b11*obs[, "sex"] + b12*obs[, "age0_c50"] + 
          b13*obs[, "U"] + b10a*(t >= knots[1] & t< knots[2]) + 
          CslopeB*(t >= knots[2] & t< knots[3]) +
          CslopeC*(t >= knots[3]))*t
    }
    slopes <- matrix(NA, nrow = num_obs, ncol= (length(visit_times) - 1))
    #Cij is stored as a list in this function environment so use list indexing
    for(j in 1:ncol(slopes)){
      b <- Cij[[j + 1]] 
      a <- Cij[[j]]
      slopes[, j] = (b-a)/int_time
    }
    return(list("Cij" = Cij, "slopes" = slopes))
  }

  #---- Generate Covariance Matrix for random slope and intercept terms ----
  slope_int_cov <- matrix(c(var0, cov, cov, var1), nrow = 2, byrow = TRUE)
  
  #---- The simulation function ----
  opt_sex_dem_sim <- function(dem_cut){
    #---- Generating IDs, sex, U ----
    obs <- tibble("id" = seq(from = 1, to = num_obs, by = 1), 
                  "sex" = rbinom(num_obs, size = 1, prob = psex), 
                  "U" = rnorm(num_obs, mean = 0, sd = 1))
    
    #---- Generating age data ----
    #Creating ages at each timepoint j
    ages <- as_tibble(matrix(NA, nrow = num_obs, ncol = length(age_varnames)))
    for(j in 1:length(age_varnames)){
      if(j == 1){
        ages[, j] = seq(from = 1, to = num_obs, by = 1) #Creates column of ids
      } else if(j == 2){
        ages[, j] = age0 #Creates column of baseline ages
      } else ages[, j] = ages[, (j-1)] + int_time #Creates ages at following timepoints
    }
    colnames(ages) <- age_varnames
    
    #---- Generating centered age data ----
    #Creating centered ages at each timepoint j
    c_ages <- as_tibble(ages - mean(age0)) %>% 
      mutate("id" = seq(from = 1, to = num_obs, by = 1)) #Creates column of ids
    colnames(c_ages) <- agec_varnames
    
    #---- Generating "true" cognitive function Cij ----
    #Cij = b00 + z0i + bo1*sexi + b02*age0i + b03*Ui + (b10 + z1i + b11*sexi + 
    #b12*age0i + b13*Ui)*timej + epsilonij
    
    #---- Generating random terms for slope and intercept ----
    #Generate random terms for each individual
    slope_int_noise <- as_tibble(mvrnorm(n = num_obs, mu = rep(0, 2), 
                                         Sigma = slope_int_cov)) %>% 
      cbind("id" = seq(from = 1, to = num_obs, by = 1), .) #Creates column of ids
    colnames(slope_int_noise) <- c("id", "z0i", "z1i")
    
    #---- Generating noise term (unexplained variance in Cij) for each visit ----
    sd_eps <- sqrt(var3)
    eps <- as_tibble(replicate(num_tests + 1, 
                               rnorm(n = num_obs, mean = 0, sd = sd_eps))) %>% 
      cbind("id" = seq(from = 1, to = num_obs, by = 1), .) #Creates column of ids
    colnames(eps) <- eps_varnames
    
    #---- Creating complete matrix of observation data ----
    obs <- left_join(obs, ages, by = "id") %>% left_join(c_ages, by = "id") %>%
      left_join(slope_int_noise, by = "id") %>% left_join(eps, by = "id")
    
    #---- Calculating Cij for each individual ----
    #Store Cij values
    Cij <- as.data.frame(opt_cog_func(CslopeB, CslopeC, obs)$Cij) %>% 
      cbind("id" = seq(from = 1, to = num_obs, by = 1), .) #Creating column of ids
    colnames(Cij) <- Cij_varnames
    
    #Store slope values per interval per individual
    slopeij <- as.data.frame(cog_func(obs)$slopes) %>% 
      cbind("id" = seq(from = 1, to = num_obs, by = 1), .) #Creating column of ids
    colnames(slopeij) <- slopeij_varnames
    
    #---- Generate survival time for each person ----
    #Individual hazard functions
    #h(tij|x) = lambda*exp(g1*sexi + g2*ageij + g3*Ui + g4*sexi + 
    #g5*slopeij + g6Cij)
    #See Additional notes in README file
    
    #---- Generating uniform random variables per interval for Sij ----
    USij <- as_tibble(replicate(num_tests, 
                                runif(num_obs, min = 0, max = 1))) %>%
      cbind("id" = seq(from = 1, to = num_obs, by = 1), .) #Creating column of ids
    colnames(USij) <- USij_varnames
    
    #---- Merging Cij, slopeij, and USij with observation data ----
    #Used as input for survival function
    obs <- left_join(obs, Cij, by = "id") %>% left_join(slopeij, by = "id") %>% 
      left_join(USij, by = "id")
    
    #---- Calculating Sij for each individual ----
    #Store Sij values
    Sij <- as.data.frame(survival(obs, lambda))
    colnames(Sij) <- Sij_varnames
    
    #---- Calculating death data for each individual ----
    #Compute death indicator for each interval
    deathij <- (Sij < int_time)*1 
    for(i in 1:nrow(deathij)){
      death <- min(which(deathij[i, ] == 1))
      if(is.finite(death)){
        deathij[i, death:ncol(deathij)] = 1 #Changes death indicators to 1 after death
      }
    }
    colnames(deathij) <- deathij_varnames
    
    #Compute study death indicators
    study_death <- (rowSums(deathij) > 0)*1
    
    #Compute overall survival times
    survtime <- vector(length = num_obs)
    survtime[which(study_death == 0)] = num_tests*int_time
    for(i in 1:length(survtime)){
      if(survtime[i] == 0){
        death_int <- min(which(deathij[i, ] == 1)) 
        survtime[i] = int_time*(death_int - 1) + Sij[i, death_int]
      } 
    }
    
    #Computing age at death
    age_death <- age0 + survtime
    
    #---- Censor Cij based on death data ----
    for(i in 1:num_obs){
      death_int <- (min(which(deathij[i, ] == 1)) - 1)
      if(is.finite(death_int)){
        Cs <- vector(length = num_tests)
        for(j in death_int:num_tests){
          Cs[j] <- paste("Ci", j, sep = "")
        }
        Cs <- Cs[Cs != "FALSE"]
        obs[i, dput(Cs)] <- NA
      }
    }
    #---- Create a competing risk outcome ----
    demij <- obs %>% dplyr::select(dput(Cij_varnames[-1])) %>% 
      mutate_all(funs((. < dem_cut)*1))
    colnames(demij) <- dem_varnames
    
    #---- Compute person years ----
    Sij[Sij > 5] <- 5
    dem_cases <- colSums(demij, na.rm = TRUE)[-1]
    person_years <- colSums(Sij, na.rm = TRUE)
    cases_py1000 <- 1000*dem_cases/person_years
  
    #---- Values to return ----
    return(cases_py1000)
  }
  py_1000_vec <- opt_sex_dem_sim(dem_cut)
  target <- py_1000_vec[J]
  
  return(abs(target - dem_val))
  }
  
  #Creating vector search to return
  dem_cut_vals <- vector("list", length = num_tests)
  for(j in 1:length(dem_cut_vals)){
    testnum = j - 1
    age = age0[1] + testnum*int_time
    index <- (which(dem_table$LowAge <= age & dem_table$HighAge > age)) - 1
    if(length(index) > 0){
      if(index > 0){
        val_match <- as.double(dem_table[index, "Rate"])
        dem_cut_vals[[j]] <- optim(par = c(-0.15, -0.4, -1.5), 
                                   fn = dem_rates, J = j, dem_val = val_match, 
                                   upper = c(0, 0, 0))
      }
    }
  }
  return(dem_cut_vals)
}

dem_cut_vals <- find_demcut(dem_rates_whites)
