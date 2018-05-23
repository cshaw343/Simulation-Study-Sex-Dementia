#---- Package loading, options, seed ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "MASS", "reshape")

options(scipen = 999) #Standard Notation
options(digits = 3)   #Round to 3 decimal places
options(warn = -1)    #Suppress warnings

set.seed(10789)

#---- Generating variable names for each assessment timepoint ----
age_varnames <- c("id", "age0", vector(length = num_tests)) #Age labels
agec_varnames <- c("id", "age0_c50", vector(length = num_tests)) #Centered age labels
eps_varnames <- c("id", "eps0", vector(length = num_tests)) #Epsilon labels
dem_varnames <- c("dem0", vector(length = num_tests)) #Dementia labels
Cij_varnames <- c("id", "Ci0", vector(length = num_tests)) #Cij labels
mean_Cij_varnames <- c("sex", "meanCi0", vector(length = num_tests)) #Mean Cij labels
slopeij_varnames <- c("id", vector(length = num_tests)) #interval slope labels
USij_varnames <- c("id", vector(length = num_tests)) #Uniform survival noise labels
deathij_varnames <- vector(length = num_tests) #Death indicator labels
Sij_varnames <- vector(length = num_tests) #Survival time labels


for(j in 1:num_tests){
  age_varnames[j + 2] = paste("age", j, sep = "")
  agec_varnames[j + 2] = paste(age_varnames[j + 2], "_c50", sep = "")
  eps_varnames[j + 2] = paste("eps", j, sep = "")
  dem_varnames[j + 1] = paste("dem", j, sep = "")
  Cij_varnames[j + 2] = paste("Ci", j, sep = "")
  mean_Cij_varnames[j + 2] = paste("meanCi", j, sep = "")
  slopeij_varnames[j + 1] = paste("slope",j-1, j, sep = "")
  USij_varnames[j + 1] = paste("U",j-1, j, sep = "")
  deathij_varnames[j] = paste("death",j-1, j, sep = "")
  Sij_varnames[j] = paste("survtime",j-1, j, sep = "")
}

#---- Generating assessment timepoint data ----
visit_times <- seq(from = 0, to = int_time*num_tests, by = int_time)

#---- Model for Cognitive Function ----
cog_func <- function(obs){
  knots = c(0, 20, 35)
  Cij <- vector(length = length(visit_times))
  for(j in 1:length(Cij)){
    t = visit_times[j]
    test_num = j - 1
    eps <- paste("eps", test_num, sep = "")
    Cij[j] = b00 + obs[, "z0i"] + b01*obs[, "sex"] + b02*obs[, "age0_c50"] + 
      b03*obs[, "U"] + obs[, eps] + 
      (b10a - b10b)*knots[2]*(t >= knots[2]) + 
      (b10b - b10c)*knots[3]*(t >= knots[3]) + 
      (obs[, "z1i"] + b11*obs[, "sex"] + b12*obs[, "age0_c50"] + 
         b13*obs[, "U"] + b10a*(t >= knots[1] & t< knots[2]) + 
                          b10b*(t >= knots[2] & t< knots[3]) +
                          b10c*(t >= knots[3]))*t
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

#---- Model for Survival Time ----
survival <- function(obs, lambda){
  #Calculate survival times for each interval
  Sij <- vector(length = (length(visit_times) - 1))
  for(j in 1:length(Sij)){
    test_num = j - 1
    US <- paste("U", test_num, test_num + 1, sep = "")
    agec <- paste("age", test_num, "_c50", sep = "")
    slope <- paste("slope", test_num, test_num + 1, sep = "")
    C <- paste("Ci", test_num, sep = "")
    Sij[j] = -log(obs[US])/
      (lambda[j]*exp(g1[j]*obs["sex"] + g2*obs[agec] + g3*obs["U"] + 
                    g4*obs["sex"]*obs[agec] + g5*obs[slope] + 
                    g6*obs[C]))
  }
  return("Sij" = Sij)
}

#---- Function for computing dementia onset ----
dem_onset <- function(obs){
  #Function we are trying to find the root of
  event_est <- function(t, M, B){
    cij = B + M*t
    return(abs(cij - dem_cut))
  }
  timetodem <- vector(length = nrow(obs))
  for(i in 1:nrow(obs)){
    data <- obs[i, ]
    if(is.na(data["dem_wave"])){
      timetodem[i] = as.double(data["survtime"])
    } else if(data["dem_wave"] == 0){
      timetodem[i] = 0
    } else {
      wave = as.double(data["dem_wave"])
      slope <- paste("slope", wave - 1, wave, sep = "")
      int <- paste("Ci", wave - 1, sep = "")
      m = as.double(data[slope])
      b = as.double(data[int])
      timetodem[i] = (wave - 1)*int_time + 
        optimize(event_est, M = m, 
                 B = b, interval = c(0, 5))$minimum
    }
  }
  return(timetodem)
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
    Cij <- as.data.frame(cog_func(obs)$Cij) %>% 
      cbind("id" = seq(from = 1, to = num_obs, by = 1), .) #Creating column of ids
    colnames(Cij) <- Cij_varnames
    
    #Store slope values per interval per individual
    slopeij <- as.data.frame(cog_func(obs)$slopes) %>% 
      cbind("id" = seq(from = 1, to = num_obs, by = 1), .) #Creating column of ids
    colnames(slopeij) <- slopeij_varnames
    
    #---- Calculating mean Cij by sex ----
    mean_Cij <- Cij %>% mutate("sex" = obs$sex) %>% 
      mutate_at("sex", as.factor) %>% group_by(sex) %>% 
      dplyr::select(-id) %>% summarise_all(mean)
    colnames(mean_Cij) <- mean_Cij_varnames
    
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
    dem_wave <- vector(length = num_obs)  #Wave at which dementia was diagnosed
    for(i in 1:nrow(demij)){
      dem_time <- min(which(demij[i, ] == 1))
        if(is.finite(dem_time)){
          demij[i, dem_time:ncol(demij)] = 1  #Changes dementia indicators to 1 after initial diagnosis
          dem_wave[i] = dem_time - 1          #Fills in wave of dementia diagnosis
        } else{
          dem_wave[i] = NA
        }
    }
    
    #Dementia diagnosis indicator
    dem <- (1 - is.na(dem_wave))
      
    #Compute time to dementia
    timetodem <- obs %>% dplyr::select(Cij_varnames, slopeij_varnames) %>% 
      cbind(., dem_wave, survtime) %>% mutate("timetodem" = dem_onset(.)) %>%
      dplyr::select("timetodem")
    
    #Age at dementia diagnosis
    ageatdem <- age0 + timetodem
    names(ageatdem) <- c("ageatdem")
    
    #Dementia status at death
    dem_death <- 
      as_tibble(cbind(dem, timetodem, survtime, study_death)) %>% 
      mutate("dem_death" = 
               case_when(dem == 1 & timetodem <= survtime ~ 1, 
                         study_death == 1 & 
                           (dem == 0 | (dem == 1 & timetodem > survtime)) 
                         ~ 2)) %>% 
      mutate_at("dem_death", funs(replace(., is.na(.), 0))) %>% 
      dplyr::select("dem_death")
    
    timetodem_death <- as_tibble(cbind(timetodem, survtime, dem)) %>% 
      mutate("timetodem_death" = 
               if_else(dem == 1, pmin(timetodem, survtime), survtime)) %>%
      dplyr::select("timetodem_death")
    
    ageatdem_death <- age0 + timetodem_death %>% 
      mutate("ageatdem_death" = timetodem_death) %>% 
      dplyr::select("ageatdem_death")
    
    dem_alive <- as_tibble((dem_death == 1)*1) %>% 
      mutate("dem_alive" = dem_death) %>% dplyr::select("dem_alive")
    
#---- Combine all variables ----
obs <- cbind(obs, Sij, deathij, study_death, demij, dem_wave, dem, 
             survtime, age_death, timetodem, ageatdem, dem_death, 
             timetodem_death, ageatdem_death, dem_alive) %>% mutate() 
    
#---- Edit for actual simulation ----
#Comment out for simulation checks
#obs <- obs %>% filter(dem_wave != 0)

#Desired simulation return values    
#return(list("mean_Cij" = mean_Cij))
    
#Alternative return function for code checking
return(list("obs" = obs, "mean_Cij" = mean_Cij))
}



