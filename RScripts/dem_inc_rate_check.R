#---- Package loading, options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "here", "magrittr", "MASS")

options(scipen = 999) #Standard Notation
options(digits = 6)   #Round to 6 decimal places
options(warn = -1)    #Suppress warnings

#---- Source files ----
source(here(
  "RScripts", "quadratic_model",
  "sex-dementia_sim_parA_onedemcut_nodemkill_maleAD_quad.R"))
source(here("RScripts", "quadratic_model", "variable_names_quad.R"))
source(here("RScripts", "quadratic_model", "sex-dementia_sim_data_gen_quad.R"))
source(here("RScripts", "dementia_incidence_ACT.R"))
source(here("RScripts", "US_life_table_calcs.R"))

#---- Plug in newly optimized data ----
# #linear splines model
# cij_slopes <- opt_cij_slopes
# cij_var1 <- opt_cij_var1

#Quadratic model
b10 <- opt_linear_term 
b20 <- opt_quadratic_term
cij_var1 <- opt_linear_var
cij_var2 <- opt_quadratic_var
dem_cut <- opt_dem_cut

#---- Generate the data ----
num_obs = 500000
obs <- data_gen(num_obs)

#---- Data by sex ----
male_data <- obs %>% filter(female == 0)
female_data <- obs %>% filter(female == 1)

#---- Compute incidence rates ----
# all_sim_inc_rates <- matrix(ncol = 9, nrow = 1)
# colnames(all_sim_inc_rates) <- variable_names$interval_ages[1:num_tests]
# rownames(all_sim_inc_rates) <- c("")
# 
# for(slot in 1:num_tests){
#   if(slot == 1){
#     dem_last_wave <- paste0("dem", (slot - 1))
#     dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
#     death_last_wave <- paste0("death", (slot - 1))
#     death_this_wave <- paste0("death", (slot - 1), "-", slot)
#     contributed <- paste0("contributed", (slot - 1), "-", slot)
#   } else {
#     dem_last_wave <- paste0("dem", (slot - 2), "-", (slot - 1))
#     dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
#     death_last_wave <- paste0("death", (slot - 2), "-", (slot - 1))
#     death_this_wave <- paste0("death", (slot - 1), "-", slot)
#     contributed <- paste0("contributed", (slot - 1), "-", slot)
#   }
#   PY_data <- obs %>% 
#     dplyr::select(death_last_wave, death_this_wave, 
#                   dem_last_wave, dem_this_wave, contributed) %>% 
#     filter(!! as.name(death_last_wave) == 0 & 
#              !! as.name(dem_last_wave) == 0) 
#   
#   all_sim_inc_rates[1, slot] = round(1000*(sum(PY_data[, dem_this_wave], 
#                                                 na.rm = TRUE)/
#                                               sum(PY_data[, contributed])), 3)
# }

male_sim_inc_rates <- matrix(ncol = 9, nrow = 1)
colnames(male_sim_inc_rates) <- variable_names$interval_ages[1:num_tests]
rownames(male_sim_inc_rates) <- c("")

for(slot in 1:num_tests){
  if(slot == 1){
    dem_last_wave <- paste0("dem", (slot - 1))
    dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
    death_last_wave <- paste0("death", (slot - 1))
    death_this_wave <- paste0("death", (slot - 1), "-", slot)
    contributed <- paste0("contributed", (slot - 1), "-", slot)
  } else {
    dem_last_wave <- paste0("dem", (slot - 2), "-", (slot - 1))
    dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
    death_last_wave <- paste0("death", (slot - 2), "-", (slot - 1))
    death_this_wave <- paste0("death", (slot - 1), "-", slot)
    contributed <- paste0("contributed", (slot - 1), "-", slot)
  }
  PY_data <- male_data %>% 
    dplyr::select(death_last_wave, death_this_wave, 
                  dem_last_wave, dem_this_wave, contributed) %>% 
    filter(!! as.name(death_last_wave) == 0 & 
             !! as.name(dem_last_wave) == 0) 
  
  male_sim_inc_rates[1, slot] = round(1000*(sum(PY_data[, dem_this_wave], 
                                       na.rm = TRUE)/
                                     sum(PY_data[, contributed])), 3)
}

female_sim_inc_rates <- matrix(ncol = 9, nrow = 1)
colnames(female_sim_inc_rates) <- variable_names$interval_ages[1:num_tests]
rownames(female_sim_inc_rates) <- c("")

for(slot in 1:num_tests){
  if(slot == 1){
    dem_last_wave <- paste0("dem", (slot - 1))
    dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
    death_last_wave <- paste0("death", (slot - 1))
    death_this_wave <- paste0("death", (slot - 1), "-", slot)
    contributed <- paste0("contributed", (slot - 1), "-", slot)
  } else {
    dem_last_wave <- paste0("dem", (slot - 2), "-", (slot - 1))
    dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
    death_last_wave <- paste0("death", (slot - 2), "-", (slot - 1))
    death_this_wave <- paste0("death", (slot - 1), "-", slot)
    contributed <- paste0("contributed", (slot - 1), "-", slot)
  }
  PY_data <- female_data %>% 
    dplyr::select(death_last_wave, death_this_wave, 
                  dem_last_wave, dem_this_wave, contributed) %>% 
    filter(!! as.name(death_last_wave) == 0 & 
             !! as.name(dem_last_wave) == 0) 
  
  female_sim_inc_rates[1, slot] = round(1000*(sum(PY_data[, dem_this_wave], 
                                                  na.rm = TRUE)/
                                              sum(PY_data[, contributed])), 3)
}

#---- Compute survival data ----

#---- Cohort size ----
num_obs <- nrow(obs)
num_males <- nrow(male_data)
num_females <- nrow(female_data)

#---- Survival probabilities ----
p_alive <- obs %>% 
  dplyr::select(variable_names$deathij_varnames[1:num_tests]) %>% 
  map_dbl(~sum(. == 0, na.rm = TRUE))/num_obs

#---- Survival probabilities by sex ----
p_alive_males <- male_data %>%
  dplyr::select(variable_names$deathij_varnames[1:num_tests]) %>% 
  map_dbl(~sum(. == 0, na.rm = TRUE))/num_males
 
p_alive_females <- female_data %>%
  dplyr::select(variable_names$deathij_varnames[1:num_tests]) %>% 
  map_dbl(~sum(. == 0, na.rm = TRUE))/num_females

#---- Checking values ----
#head(EURODEM_inc_rates$Total_All_Dementia_1000PY, -1)
head(ACT_inc_rates$Male_AD_1000PY)
male_sim_inc_rates

male_life_US$cum_surv_cond50[-1]
p_alive_males

female_life_US$cum_surv_cond50[-1]
p_alive_females

# #Linear splines model
# #Calculate observed slopes
# Cij_check <- obs %>% group_by(female) %>% 
#   dplyr::select(c("female", variable_names$Cij_varnames)) %>% 
#   summarise_all(mean, na.rm = TRUE)
# 
# slopes_check <- matrix(nrow = 2, ncol = (length(visit_times)))
# colnames(slopes_check) = c("female", variable_names$interval_ages[1:9])
# slopes_check[, "female"] = Cij_check$female
# 
# for(i in 2:ncol(slopes_check)){
#   slopes_check[, i] <- as.matrix((Cij_check[, i + 1] - Cij_check[, i])/int_time)
# }

# slopes_check 
# colMeans(slopes_check)
