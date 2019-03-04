#*******************************************************************************
# This is the function that performs the search for appropriate slope values in 
# the model for Cij as well as age-specific dementia cut-off points
#
# This function is attempting to find slopes that will create dementia incidence
# rates that match the EURODEM from Table 2 of The EURODEM Studies 
#(Andersen 1999)
#
#*******************************************************************************

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "magrittr", "MASS", "optimParallel")

options(scipen = 999)

#---- Source Files ----
source("RScripts/dementia_incidence_EURODEM_pooled.R")
source("RScripts/sex-dementia_sim_parA.R")
source("RScripts/variable_names.R")
source("RScripts/cognitive_function_model.R")
source("RScripts/survival_times.R")
source("RScripts/dementia_onset.R")
source("RScripts/misc_custom_functions.R")

#---- Data generation ----
#Can generate data all the way up to generating cognitive function
generate_base_data <- function(n){
  #---- Create a blank dataset ----
  obs <- matrix(NA, nrow = n, ncol = length(column_names)) %>% 
    as.data.frame() %>% set_colnames(column_names)
  
  #---- Generating IDs, sex, U ----
  obs$id <- seq(from = 1, to = n, by = 1)
  obs$sex <- rbinom(n, size = 1, prob = psex)
  obs$U <- rnorm(n, mean = 0, sd = 1)
  
  #---- Generating age data ----
  #Creating ages at each timepoint j
  for(j in 1:length(variable_names$age_varnames)){
    if(j == 1){
      obs[, variable_names$age_varnames[j]] = age0 #Creates column of baseline ages
    } else 
      obs[, variable_names$age_varnames[j]] = 
        obs[, variable_names$age_varnames[j - 1]] + int_time #Creates ages at following timepoints
  }
  
  #---- Generating centered age data ----
  #Creating baseline-mean-centered ages at each timepoint j
  obs[, variable_names$agec_varnames] <- 
    obs[, variable_names$age_varnames] - mean(age0)
  
  #---- Generating "true" cognitive function Cij ----
  #Refer to Manuscript/manuscript_equations.pdf for equation
  
  #Generating random terms for slope and intercept
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
    noise <- mvrnorm(n = n, mu = rep(0, 2), 
                     Sigma = cij_slope_int_cov[[i]]) 
    obs[, c(paste0("z0_", (i - 1), "i"), paste0("z1_", (i - 1), "i"))] <- noise
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
    mvrnorm(n = n, mu = rep(0, num_visits), Sigma = cij_cov_mat)
  
  return(obs)
}

#---- The function we want to optimize ----
#This function should optimize the incidence rate up to a certain point
dem_irate_1000py <- function(NEWSLOPE_NEWDEMCUT, 
                             old_slopes = NA, old_demcuts = NA, 
                             age, pub_inc, obs){
  half_point <- length(NEWSLOPE_NEWDEMCUT)/2
  if(is.na(old_slopes)){
    slopes <- c(NEWSLOPE_NEWDEMCUT[1:half_point], 
                rep(0, (num_tests - half_point)))
    demcuts <- c(NEWSLOPE_NEWDEMCUT[(half_point + 1)], #for baseline Cij
                 NEWSLOPE_NEWDEMCUT[(half_point + 1):length(NEWSLOPE_NEWDEMCUT)], 
                 rep(0, (num_tests - half_point)))
  } else {
    slopes <- c(old_slopes, NEWSLOPE_NEWDEMCUT[1], 
                rep(0, (num_tests - (1 + length(old_slopes)))))
    demcuts <- c(0, #for baseline Cij
                 old_demcuts, NEWSLOPE_NEWDEMCUT[2], 
                 rep(0, (num_tests - (1 + length(old_demcuts)))))
  }
  
  #Calculating Cij for each individual
  #Store Cij values and slope values for each assessment
  compute_Cij <- cog_func(cij_knots, cij_slopes, obs)
  obs[, variable_names$Cij_varnames] <- compute_Cij$Cij
  obs[, na.omit(variable_names$cij_slopeij_varnames)] <- compute_Cij$slopes
  
  #---- Generate survival time for each person ----
  #Refer to Manuscript/manuscript_equations.pdf for equation
  
  #---- Generating uniform random variables per interval for Sij ----
  obs[, na.omit(variable_names$rij_varnames)]<- 
    replicate(num_tests, runif(small_batch_n, min = 0, max = 1))
  
  #---- Calculating Sij for each individual ----
  #Store Sij values and survival time
  survival_data <- survival(obs)
  obs[, na.omit(variable_names$Sij_varnames)] <- survival_data$Sij
  obs[, "survtime"] <- survival_data$survtimes
  
  #---- Calculating death data for each individual ----
  #Indicator of 1 means the individual died in that interval
  #NAs mean the individual died in a prior interval
  obs[, "death0"] <- 0
  obs[, na.omit(variable_names$deathij_varnames)] <- 
    (obs[, na.omit(variable_names$Sij_varnames)] < int_time)*1 
  
  obs[, "study_death"] <- 
    rowSums(obs[, na.omit(variable_names$deathij_varnames)], na.rm = TRUE) #Study death indicator
  
  obs[, "age_death"] <- obs[, "age0"] + obs[, "survtime"]
  
  # #---- Standardize Cij values ----
  # std_Cij <- obs %>% dplyr::select(variable_names$Cij_varnames) %>%
  #   map_df(~(. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE)) %>%
  #   set_colnames(variable_names$std_Cij_varnames)
  # 
  # obs %<>% bind_cols(., std_Cij)
  
  #---- Create a competing risk outcome ----
  dem_cuts_mat <- matrix(dem_cuts, nrow = nrow(obs), ncol = length(dem_cuts), 
                         byrow = TRUE)
  
  obs[, variable_names$dem_varnames] <- 
    (obs[, variable_names$Cij_varnames] < dem_cuts_mat)*1
  
  obs %<>% filter(dem0 == 0)
  
  #---- Survival censoring matrix ----
  censor <- (obs[, na.omit(variable_names$Sij_varnames)] == 5)*1
  censor[censor == 0] <- NA
  censor %<>% cbind(1, .)
  
  #---- Censor Cij and dem data ----
  obs[, variable_names$Cij_varnames] <- 
    obs[, variable_names$Cij_varnames]*censor
  
  obs[, variable_names$dem_varnames[-1]] <- 
    obs[, variable_names$dem_varnames]*censor[, 1:(ncol(censor) - 1)]
  
  #---- Dementia indicators ----
  for(i in 1:nrow(obs)){
    dem_int <- min(which(obs[i, variable_names$dem_varnames] == 1))
    if(is.finite(dem_int)){
      obs[i, "dem_wave"] <- (dem_int - 1)
      first_censor <- min(which(is.na(obs[i, variable_names$dem_varnames])))
      if(dem_int < 10 & is.finite(first_censor)){
        obs[i, variable_names$dem_varnames[dem_int:(first_censor - 1)]] <- 1 #Changes dementia indicator to 1 after dementia diagnosis
      }
      if(dem_int < 10 & !is.finite(first_censor)){
        last_1 <- length(variable_names$dem_varnames)
        obs[i, variable_names$dem_varnames[dem_int:last_1]] <- 1 #Changes dementia indicator to 1 after dementia diagnosis
      }
    } else {
      obs[i, "dem_wave"] = NA
    }
  }
  
  #---- Dementia calcs ----
  obs[, "dem"] <- (1 - is.na(obs[, "dem_wave"])) #Dementia diagnosis indicator
  obs[, "timetodem"] <- dem_onset(obs, dem_cuts) #Time to dementia diagnosis
  obs <- compare_survtime_timetodem(obs)
  obs[, "ageatdem"] <- obs[, "age0"] + obs[, "timetodem"] #Age at dementia diagnosis
  
  #Dementia status at death
  for(i in 1:nrow(obs)){
    if(obs[i, "dem"] == 1 & obs[i, "timetodem"] <= obs[i, "survtime"]){
      obs[i, "dem_death"] <- 1
    } else if(obs[i, "study_death"] == 1 & 
              (obs[i, "dem"] == 0 | (obs[i, "dem"] == 1 & 
                                     obs[i, "timetodem"] > obs[i, "survtime"]))){
      obs[i, "dem_death"] <- 2
    } else {
      obs[i, "dem_death"] <- 0
    }
  }
  
  #Time to dem_death
  for(i in 1:nrow(obs)){
    if(obs[i, "dem"] == 0){
      obs[i, "timetodem_death"] <- obs[i, "survtime"]
    } else {
      obs[i, "timetodem_death"] <- min(obs[i, "timetodem"], obs[i, "survtime"])
    }
  }
  
  obs[, "ageatdem_death"] <- obs[, "age0"] + obs[, "timetodem_death"]
  obs[obs[, "dem_death"] == 1, "dem_alive"] <- 1
  obs[is.na(obs[, "dem_alive"]), "dem_alive"] <- 0
  
  #---- Compute person years ----
  contributed <- (obs$timetodem_death)%%5
  slot <- (age - 50)/5
  dem_last_wave <- paste("dem", (slot - 1), sep = "")
  dem_this_wave <- paste("dem", slot, sep = "")
  death_last_wave <- paste("death", (slot - 2), "-", (slot - 1), sep = "")
  death_this_wave <- paste("death", (slot - 1), "-", (slot), sep = "")
  
  PY_data <- obs %>% dplyr::select(death_last_wave, death_this_wave, 
                                   dem_last_wave, dem_this_wave) %>% 
    cbind(., contributed) %>% 
    filter(!! as.name(death_last_wave) == 0 & 
             !! as.name(dem_last_wave) == 0) %>% 
    mutate("person_years" = 
             case_when(!! as.name(death_this_wave) == 1 | 
                         !! as.name(dem_this_wave) == 1 ~ contributed, 
           TRUE ~ 5))
  
  cases_py1000 = 
    1000*sum(PY_data[, dem_this_wave], na.rm = TRUE)/sum(PY_data$person_years)

  return(abs(cases_py1000 - pub_inc))
}

#---- Running the optimization ----
dem_inc_table <- EURODEM_inc_rates

best_slopes_cuts <- matrix(nrow = num_tests, ncol = 4) %>% as.data.frame()
colnames(best_slopes_cuts) <- c("age", "slope", "dem_cut", #"dem_inc_rate", 
                                "diff")
best_slopes_cuts[, "age"] <- seq(55, 100, by = 5)

#Setting up cluster for parallel optimization
#If using a Mac/Linux system, it's highly recommended to use the type = "FORK"
#option instead and comment out the clusterEvalQ() lines

cluster <- makeCluster(detectCores() - 2, type = "PSOCK")
clusterEvalQ(cl = cluster, {
  if (!require("pacman")) 
    install.packages("pacman", repos='http://cran.us.r-project.org')
  
  p_load("tidyverse", "magrittr", "MASS")
  
  source("RScripts/sex-dementia_sim_parA.R") 
  source("RScripts/variable_names.R") 
  source("RScripts/cognitive_function_model.R") 
  source("RScripts/survival_times.R") 
  source("RScripts/dementia_onset.R")
  source("Rscripts/misc_custom_functions.R")
  
}) 

clusterExport(cl = cluster, 
              varlist = c("dem_inc_table", "best_slopes_cuts"), 
              envir = environment())

#Finding the parameters based on the young cohort
slopes_cuts = c(seq(0, 0.05, by = 0.015)*-1, seq(2.75, 4.25, by = 0.5)*-1)
first_search <- 
  replicate(10, 
            optimParallel(par = slopes_cuts, fn = dem_irate_1000py, 
                          age = dem_inc_table[[1, "Visit_Age"]], 
                          pub_inc = 
                            dem_inc_table[[1, "Total_All_Dementia_1000PY"]], 
                          obs = generate_base_data(n = 10000),
                          upper = c(rep(0, 4), rep(-2.5, 4)), 
                          lower = c(rep(-0.05, 4), rep(-5.5, 4)), 
                          parallel = list(cl = cluster)))
avg_first_pars <- as_tibble(do.call(rbind, first_search["par", ])) %>%
  colMeans()
avg_first_diffs <- as_tibble(do.call(rbind, first_search["value", ])) %>%
  colMeans()
first_slopes <- avg_first_pars[1:4]
first_dem_cuts <- avg_first_pars[5:8] %>% force_dec()
best_slopes_cuts[1:(length(avg_first_pars)/2), "slope"] <- 
  first_slopes
best_slopes_cuts[1:(length(avg_first_pars)/2), "dem_cut"] <- 
  first_dem_cuts
best_slopes_cuts[length(avg_first_pars)/2, "diff"] <- avg_first_diffs
write_csv(best_slopes_cuts[1:(length(avg_first_pars)/2), ], 
          paste0("Data/best_slopes_cuts_", gsub("-", "", Sys.Date()), ".csv"))

#Finding parameters based on the older cohorts
#Trying to repeat the optimization many times and take an average (like lambda search)

#For Visit Age 75:
last_slot <- max(which(!is.na(best_slopes_cuts[, "slope"])))
this_slot <- last_slot + 1
index <- this_slot - 3
slopes_cuts = c(0, (best_slopes_cuts[[last_slot, "dem_cut"]] - 0.5))
search_75 <- 
  replicate(5, 
            optimParallel(par = slopes_cuts, fn = dem_irate_1000py, 
                          age = dem_inc_table[[index, "Visit_Age"]], 
                          pub_inc = 
                            dem_inc_table[[index, "Total_All_Dementia_1000PY"]], 
                          obs = generate_base_data(n = 20000), 
                          old_slopes = best_slopes_cuts[1:max(which(!is.na(
                                       best_slopes_cuts[, "slope"]))), "slope"], 
                          old_demcuts = best_slopes_cuts[1:max(which(!is.na(
                            best_slopes_cuts[, "dem_cut"]))), "dem_cut"], 
                          upper = c(0, best_slopes_cuts[[last_slot, "dem_cut"]]), 
                          lower = c(-0.15, best_slopes_cuts[[last_slot, "dem_cut"]] - 1), 
                          parallel = list(cl = cluster)))

avg_pars_75 <- as_tibble(do.call(rbind, search_75["par", ])) %>%
  colMeans()
avg_diffs_75 <- as_tibble(do.call(rbind, search_75["value", ])) %>%
  colMeans()
best_slopes_cuts[this_slot, 2:4] <- c(avg_pars_75, avg_diffs_75)
write_csv(best_slopes_cuts[this_slot, ], 
          paste0("Data/best_slopes_cuts_", gsub("-", "", Sys.Date()), ".csv"), 
          append = TRUE)

#For Visit Age 80:
last_slot <- max(which(!is.na(best_slopes_cuts[, "slope"])))
this_slot <- last_slot + 1
index <- this_slot - 3
slopes_cuts = c(0, (best_slopes_cuts[[last_slot, "dem_cut"]] - 0.5))
search_80 <- 
  replicate(5, 
            optimParallel(par = slopes_cuts, fn = dem_irate_1000py, 
                          age = dem_inc_table[[index, "Visit_Age"]], 
                          pub_inc = 
                            dem_inc_table[[index, "Total_All_Dementia_1000PY"]], 
                          obs = generate_base_data(n = 20000), 
                          old_slopes = best_slopes_cuts[1:max(which(!is.na(
                            best_slopes_cuts[, "slope"]))), "slope"], 
                          old_demcuts = best_slopes_cuts[1:max(which(!is.na(
                            best_slopes_cuts[, "dem_cut"]))), "dem_cut"], 
                          upper = c(0, (best_slopes_cuts[[last_slot, "dem_cut"]])), 
                          lower = c(-0.15, best_slopes_cuts[[last_slot, "dem_cut"]] - 1), 
                          parallel = list(cl = cluster)))

avg_pars_80 <- as_tibble(do.call(rbind, search_80["par", ])) %>%
  colMeans()
avg_diffs_80 <- as_tibble(do.call(rbind, search_80["value", ])) %>%
  colMeans()
best_slopes_cuts[this_slot, 2:4] <- c(avg_pars_80, avg_diffs_80)
write_csv(best_slopes_cuts[this_slot, ], 
          paste0("Data/best_slopes_cuts_", gsub("-", "", Sys.Date()), ".csv"), 
          append = TRUE)

#For Visit Age 85:
last_slot <- max(which(!is.na(best_slopes_cuts[, "dem_cut"])))
this_slot <- last_slot + 1
index <- this_slot - 3
slopes_cuts = c(0, (best_slopes_cuts[[last_slot, "dem_cut"]] - 0.5))
search_85 <- 
  replicate(5, 
            optimParallel(par = slopes_cuts, fn = dem_irate_1000py, 
                          age = dem_inc_table[[index, "Visit_Age"]], 
                          pub_inc = 
                            dem_inc_table[[index, "Total_All_Dementia_1000PY"]], 
                          obs = generate_base_data(n = 20000), 
                          old_slopes = best_slopes_cuts[1:max(which(!is.na(
                            best_slopes_cuts[, "slope"]))), "slope"], 
                          old_demcuts = best_slopes_cuts[1:max(which(!is.na(
                            best_slopes_cuts[, "dem_cut"]))), "dem_cut"], 
                          upper = c(0, best_slopes_cuts[[last_slot, "dem_cut"]]), 
                          lower = c(-0.15, (best_slopes_cuts[[last_slot, "dem_cut"]] - 1)), 
                          parallel = list(cl = cluster)))

avg_pars_85 <- as_tibble(do.call(rbind, search_85["par", ])) %>%
  colMeans()
avg_diffs_85 <- as_tibble(do.call(rbind, search_85["value", ])) %>%
  colMeans()
best_slopes_cuts[this_slot, 2:4] <- c(avg_pars_85, avg_diffs_85)
write_csv(best_slopes_cuts[this_slot, ], 
          paste0("Data/best_slopes_cuts_", gsub("-", "", Sys.Date()), ".csv"), 
          append = TRUE)

#For Visit Age 90:
last_slot <- max(which(!is.na(best_slopes_cuts[, "dem_cut"])))
this_slot <- last_slot + 1
index <- this_slot - 3
slopes_cuts = c(0, (best_slopes_cuts[[last_slot, "dem_cut"]] - 0.5))
search_90 <- 
  replicate(1, 
            optimParallel(par = slopes_cuts, fn = dem_irate_1000py, 
                          age = dem_inc_table[[index, "Visit_Age"]], 
                          pub_inc = 
                            dem_inc_table[[index, "Total_All_Dementia_1000PY"]], 
                          obs = generate_base_data(n = 100000), 
                          old_slopes = best_slopes_cuts[1:max(which(!is.na(
                            best_slopes_cuts[, "slope"]))), "slope"], 
                          old_demcuts = best_slopes_cuts[1:max(which(!is.na(
                            best_slopes_cuts[, "dem_cut"]))), "dem_cut"], 
                          upper = c(0, (best_slopes_cuts[[last_slot, "dem_cut"]])), 
                          lower = c(-0.4,(best_slopes_cuts[[last_slot, "dem_cut"]] - 1)), 
                          parallel = list(cl = cluster)))

avg_pars_90 <- as_tibble(do.call(rbind, search_90["par", ])) %>%
  colMeans()
avg_diffs_90 <- as_tibble(do.call(rbind, search_90["value", ])) %>%
  colMeans()
best_slopes_cuts[this_slot, 2:4] <- c(avg_pars_90, avg_diffs_90)
write_csv(best_slopes_cuts[this_slot, ], 
          paste0("Data/best_slopes_cuts_", gsub("-", "", Sys.Date()), ".csv"), 
          append = TRUE)

#For Visit Age 95:
last_slot <- max(which(!is.na(best_slopes_cuts[, "dem_cut"])))
this_slot <- last_slot + 1
index <- this_slot - 3
slopes_cuts = c(0, (best_slopes_cuts[[last_slot, "dem_cut"]] - 0.5))
search_95 <- 
  replicate(1, 
            optimParallel(par = slopes_cuts, fn = dem_irate_1000py, 
                          age = dem_inc_table[[index, "Visit_Age"]], 
                          pub_inc = 
                            dem_inc_table[[index, "Total_All_Dementia_1000PY"]], 
                          obs = generate_base_data(n = 100000), 
                          old_slopes = best_slopes_cuts[1:max(which(!is.na(
                            best_slopes_cuts[, "slope"]))), "slope"], 
                          old_demcuts = best_slopes_cuts[1:max(which(!is.na(
                            best_slopes_cuts[, "dem_cut"]))), "dem_cut"], 
                          upper = c(0, (best_slopes_cuts[[last_slot, "dem_cut"]])), 
                          lower = c(-0.4,(best_slopes_cuts[[last_slot, "dem_cut"]] - 1)), 
                          parallel = list(cl = cluster)))

avg_pars_95 <- as_tibble(do.call(rbind, search_95["par", ])) %>%
  colMeans()
avg_diffs_95 <- as_tibble(do.call(rbind, search_95["value", ])) %>%
  colMeans()
best_slopes_cuts[this_slot, ] <- c(age, avg_pars_95, avg_diffs_95)
write_csv(best_slopes_cuts[this_slot, ], 
          paste0("Data/best_slopes_cuts_", gsub("-", "", Sys.Date()), ".csv"), 
          append = TRUE)


#---- testing code ----
obs <- generate_base_data(n = 1000)
NEWSLOPE <- rep(0, 4)
NEWDEMCUT <- rep(-4.5, 4)
NEWSLOPE_NEWDEMCUT <- c(NEWSLOPE, NEWDEMCUT)
old_slopes <- NA
old_demcuts <- NA
age <- 70
pub_inc <- 3.15

test <- dem_irate_1000py (NEWSLOPE_NEWDEMCUT = NEWSLOPE_NEWDEMCUT, 
                          age = age, pub_inc = pub_inc, obs = young_cohort)

slopes_cuts = c(rep(0, 4), rep(-5, 4))

test_optim <- optim(par = slopes_cuts, fn = dem_irate_1000py, 
                    age = 70, pub_inc = 3.45, obs = cohort,
                    upper = slopes_cuts, 
                    lower = c(rep(-0.15, 4), rep(-6, 4)))

