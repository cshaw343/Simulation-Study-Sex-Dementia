#*******************************************************************************
# This is the function that performs the search for appropriate slope values in 
# the model for Cij as well as age-specific dementia cut-off points
#
# This function is attempting to find slopes that will create dementia incidence
# rates that match the Supplemental Table 2 in Inequalities in dementia 
# incidence between six racial and ethnic groups over 14 years 
# (Mayeda et al, 2016)
#
# We're starting with this US data, but may move to a different country's data
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

#---- Data generation ----
#Can generate data all the way up to generating cognitive function
generate_base_data <- function(n){
  #---- Generating IDs, sex, U ----
  obs <- tibble("id" = seq(from = 1, to = n, by = 1),
                "sex" = rbinom(n, size = 1, prob = psex), 
                "U" = rnorm(n, mean = 0, sd = 1), 
                "age0" = rep(50, n))
  
  #---- Generating age data ----
  #Creating ages at each timepoint j
  ages <- as_tibble(matrix(NA, nrow = n, 
                           ncol = length(nrow(variable_names))))
  for(j in 1:nrow(variable_names)){
    if(j == 1){
      ages[, j] = obs$age0 #Creates column of baseline ages
    } else ages[, j] = ages[, (j-1)] + int_time #Creates ages at following timepoints
  }
  colnames(ages) <- variable_names$age_varnames
  obs %<>% bind_cols(., ages)
  
  #---- Generating centered age data ----
  #Creating centered ages at each timepoint j
  c_ages <- as_tibble(ages - mean(age0)) %>% 
    set_colnames(., variable_names$agec_varnames)
  obs %<>% bind_cols(., c_ages)
  
  #---- Generating "true" cognitive function Cij ----
  #Refer to Manuscript/manuscript_equations.pdf for equation
  
  #Generating random terms for slope and intercept
  #Covariance matrix for random slope and intercept terms
  cij_slope_int_cov <- matrix(c(cij_var0, cij_cov, cij_cov, cij_var1), 
                              nrow = 2, byrow = TRUE)
  
  #Generate random terms for each individual
  cij_slope_int_noise <- as_tibble(mvrnorm(n = n, mu = rep(0, 2), 
                                           Sigma = cij_slope_int_cov)) %>% 
    set_colnames(., c("z0i", "z1i"))
  obs %<>% bind_cols(., cij_slope_int_noise)
  
  #Generating noise term (unexplained variance in Cij) for each visit
  #Creating AR(1) correlation matrix
  num_visits = num_tests + 1
  powers <- abs(outer(1:(num_visits), 1:(num_visits), "-")) #Exponents for autoregressive terms
  corr <- sqrt(cij_var3)*(cij_r1^powers)                    #Correlation matrix
  S <- diag(rep(sqrt(cij_var3)), nrow(corr))                #Diagonal matrix of SDs
  cij_cov_mat <- S%*%corr%*%S                               #Covariance matrix
  
  #Generating noise terms
  eps <- as_tibble(mvrnorm(n = n, 
                           mu = rep(0, num_visits), Sigma = cij_cov_mat)) %>%
    set_colnames(., variable_names$eps_varnames)
  obs %<>% bind_cols(., eps)
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
  compute_Cij <- cog_func(cij_knots, slopes, obs)
  Cij <- as.data.frame(compute_Cij$Cij) %>% 
    set_colnames(., variable_names$Cij_varnames)
  cij_slopeij <- as.data.frame(compute_Cij$slopes) %>% 
    #remove the last variable name because there are only 10 intervals
    set_colnames(., head(variable_names$cij_slopeij_varnames, -1)) 
  obs %<>% bind_cols(., Cij, cij_slopeij)
  
  #---- Generate survival time for each person ----
  #Refer to Manuscript/manuscript_equations.pdf for equation
  
  #---- Generating uniform random variables per interval for Sij ----
  rij <- as_tibble(replicate(num_tests, 
                             runif(n = nrow(obs), min = 0, max = 1))) %>%
    #remove the last variable name because there are only 10 intervals
    set_colnames(head(variable_names$rij_varnames, -1))
  obs %<>% bind_cols(., rij)
  
  #---- Calculating Sij for each individual ----
  #Store Sij values
  Sij <- as.data.frame(survival(obs, lambda)) %>% 
    set_colnames(head(variable_names$Sij_varnames, -1))
  obs %<>% bind_cols(., Sij)
  
  #---- Calculating death data for each individual ----
  #Compute death indicator for each interval
  #Change death indicators to 1 after death
  deathij <- as.tibble((Sij < int_time)*1) %>% 
    set_colnames(head(variable_names$deathij_varnames, -1))
  
  for(i in 1:nrow(deathij)){
    death <- min(which(deathij[i, ] == 1))
    if(is.finite(death)){
      deathij[i, death:ncol(deathij)] = 1 
    }
  }
  obs %<>% bind_cols(., deathij) %>% 
    mutate("study_death" = (rowSums(deathij) > 0)*1) #Study death indicator
  
  #Compute overall survival times
  survtime <- vector(length = nrow(obs))
  survtime[which(obs$study_death == 0)] = num_tests*int_time
  for(i in 1:length(survtime)){
    if(survtime[i] == 0){
      death_int <- min(which(deathij[i, ] == 1)) 
      survtime[i] = int_time*(death_int - 1) + Sij[i, death_int]
    } 
  }
  obs %<>% mutate("survtime" = survtime, 
                  "age_death" = age0 + survtime) #Age at death
  
  #---- Censor Cij based on death data ----
  for(i in 1:nrow(obs)){
    if(obs[i, "study_death"] == 1){
      death_int <- (min(which(deathij[i, ] == 1)))
      Cs <- c(variable_names$Cij_varnames[(death_int + 1):nrow(variable_names)])
      obs[i, Cs] <- NA
    }
  }
  
  #---- Standardize Cij values ----
  std_Cij <- obs %>% dplyr::select(variable_names$Cij_varnames) %>%
    map_df(~(. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE)) %>%
    set_colnames(variable_names$std_Cij_varnames)

  obs %<>% bind_cols(., std_Cij)
  
  #---- Create a competing risk outcome ----
  dem_cuts_mat <- matrix(demcuts, nrow = nrow(obs), ncol = length(demcuts), 
                         byrow = TRUE)
  
  demij <- obs %>% dplyr::select(variable_names$std_Cij_varnames) %>% 
    set_colnames(variable_names$dem_varnames)
  demij <- (demij < dem_cuts_mat)*1

  dem_wave <- vector(length = nrow(obs))  #Wave at which dementia was diagnosed
  for(i in 1:nrow(demij)){
    dem_time <- min(which(demij[i, ] == 1))
    if(is.finite(dem_time)){
      demij[i, dem_time:ncol(demij)] = 1  #Changes dementia indicators to 1 after initial diagnosis
      dem_wave[i] = dem_time - 1          #Fills in wave of dementia diagnosis
    } else{
      dem_wave[i] = NA
    }
  }

  obs %<>% cbind(., demij) %>%
    mutate("dem_wave" = dem_wave) %>%
    mutate("dem" = (1 - is.na(dem_wave)), #Dementia diagnosis indicator
           "timetodem" = dem_onset(., demcuts),    #Time to dementia diagnosis
           "ageatdem" = age0 + timetodem, #Age at dementia diagnosis
           "dem_death" =                  #Dementia status at death
             case_when(dem == 1 & timetodem <= survtime ~ 1,
                       study_death == 1 &
                         (dem == 0 | (dem == 1 & timetodem > survtime)) ~
                         2)) %>%
    mutate_at("dem_death", funs(replace(., is.na(.), 0))) %>%
    mutate("timetodem_death" = if_else(dem == 1, pmin(timetodem, survtime),
                                       survtime),
           "ageatdem_death" = age0 + timetodem_death,
           "dem_alive" = case_when(dem_death == 1 ~ 1,
                                   TRUE ~ 0))

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
young_cohort <- generate_base_data(n = 5000)
old_cohort <- generate_base_data(n = 50000)
very_old_cohort <- generate_base_data(n = 100000)

best_slopes_cuts <- matrix(nrow = num_tests, ncol = 4) %>% as.data.frame()
colnames(best_slopes_cuts) <- c("age", "slope", "dem_cut", #"dem_inc_rate", 
                                "diff")
best_slopes_cuts[, "age"] <- seq(55, 100, by = 5)

#Setting up cluster for parallel optimization
#If using a Mac/Linux system, it's highly recommended to use the type = "FORK"
#option instead and comment out the clusterEvalQ() lines

cluster <- makeCluster(detectCores() - 5, type = "PSOCK")
clusterEvalQ(cl = cluster, {
  source("RScripts/sex-dementia_sim_parA.R") 
  source("RScripts/variable_names.R") 
  source("RScripts/cognitive_function_model.R") 
  source("RScripts/survival_times.R") 
  source("RScripts/dementia_onset.R")
  
  if (!require("pacman")) 
    install.packages("pacman", repos='http://cran.us.r-project.org')
  
  p_load("tidyverse", "magrittr", "MASS")
}) 

clusterExport(cl = cluster, 
              varlist = c("dem_inc_table", "young_cohort", "old_cohort", 
                          "very_old_cohort", "best_slopes_cuts"), 
              envir = environment())

#Finding the parameters based on the young cohort
slopes_cuts = c(rep(0, 4), rep(-2.25, 4))
opt <- optimParallel(par = slopes_cuts, fn = dem_irate_1000py, 
                     age = dem_inc_table[[1, "Visit_Age"]], 
                     pub_inc = dem_inc_table[[1, "Total_All_Dementia_1000PY"]], 
                     obs = young_cohort,
                     upper = c(rep(0, 4), rep(-1.5, 4)), 
                     lower = c(rep(-0.15, 4), rep(-3, 4)), 
                     parallel = list(cl = cluster))
parameters <- opt$par
best_slopes_cuts[1:(length(parameters)/2), "slope"] <- 
  parameters[1:(length(parameters)/2)]
best_slopes_cuts[1:(length(parameters)/2), "dem_cut"] <- 
  parameters[(length(parameters)/2 + 1):length(parameters)]
#best_slopes_cuts[length(parameters)/2, "dem_inc_rate"] <- 
  #dem_inc_table[[1, "Total_All_Dementia_1000PY"]] + opt$value
best_slopes_cuts[length(parameters)/2, "diff"] <- opt$value
write_csv(best_slopes_cuts[1:(length(parameters)/2), ], 
          "Results/slopes_dem-cut_search.csv")

#Finding parameters based on the older cohorts
for(i in 2:nrow(dem_inc_table)){
  if(i <= 4){
    last_slot <- max(which(!is.na(best_slopes_cuts[, "dem_cut"])))
    this_slot <- last_slot + 1
    slopes_cuts = c(-0.05, (best_slopes_cuts[[last_slot, "dem_cut"]]))
    opt <- optimParallel(par = slopes_cuts, fn = dem_irate_1000py, 
                         age = dem_inc_table[[i, "Visit_Age"]], 
                         pub_inc = 
                           dem_inc_table[[i, "Total_All_Dementia_1000PY"]], 
                         obs = old_cohort, 
                         old_slopes = best_slopes_cuts[1:max(which(!is.na(
                           best_slopes_cuts[, "slope"]))), "slope"], 
                         old_demcuts = best_slopes_cuts[1:max(which(!is.na(
                           best_slopes_cuts[, "dem_cut"]))), "dem_cut"],
                         upper = c(0, (slopes_cuts[2] + 0.5)), 
                         lower = c(-0.15, slopes_cuts[2]), 
                         parallel = list(cl = cluster))
    parameters <- opt$par
    best_slopes_cuts[this_slot, "slope"] <- 
      parameters[1]
    best_slopes_cuts[this_slot, "dem_cut"] <- 
      parameters[2]
    #best_slopes_cuts[this_slot, "dem_inc_rate"] <- 
      #dem_inc_table[[i, "Total_All_Dementia_1000PY"]] + opt$value
    best_slopes_cuts[this_slot, "diff"] <- opt$value
    write_csv(best_slopes_cuts[this_slot, ], 
              "Results/slopes_dem-cut_search.csv", append = TRUE)
  } else {
    last_slot <- max(which(!is.na(best_slopes_cuts[, "dem_cut"])))
    slopes_cuts = c(-0.15, best_slopes_cuts[[last_slot, "dem_cut"]])
    opt <- optimParallel(par = slopes_cuts, fn = dem_irate_1000py, 
                         age = dem_inc_table[[i, "Visit_Age"]], 
                         pub_inc = 
                           dem_inc_table[[i, "Total_All_Dementia_1000PY"]], 
                         obs = very_old_cohort, 
                         old_slopes = best_slopes_cuts[1:max(which(!is.na(
                           best_slopes_cuts[, "slope"]))), "slope"], 
                         old_demcuts = best_slopes_cuts[1:max(which(!is.na(
                           best_slopes_cuts[, "dem_cut"]))), "dem_cut"],
                         upper = c(0, (slopes_cuts[2] + 1)),
                         lower = c(-0.40, (slopes_cuts[2])), 
                         parallel = list(cl = cluster))
    parameters <- opt$par
    best_slopes_cuts[this_slot, "slope"] <- 
      parameters[1]
    best_slopes_cuts[this_slot, "dem_cut"] <- 
      parameters[2]
    #best_slopes_cuts[this_slot, "dem_inc_rate"] <- 
      #dem_inc_table[[i, "Total_All_Dementia_1000PY"]] + opt$value
    best_slopes_cuts[this_slot, "diff"] <- opt$value
    write_csv(best_slopes_cuts[this_slot, ], 
              "Results/slopes_dem-cut_search.csv", append = TRUE)
  }
}

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

