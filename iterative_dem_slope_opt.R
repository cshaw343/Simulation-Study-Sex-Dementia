#---- Package loading, options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("here")

#---- Source files ----
source(here("RScripts", "sex-dementia_sim_parA.R"))
source(here("RScripts", "dementia_incidence_EURODEM_pooled.R"))
source(here("RScripts", "euro_life_tables.R"))
source(here("RScripts", "variable_names.R"))
source(here("RScripts", "cognitive_function_model.R"))
source(here("RScripts", "survival_times.R"))
source(here("RScripts", "dementia_onset.R"))
source(here("RScripts", "compare_survtime_timetodem.R"))

#---- Pre-allocation ----
opt_cij_slopes <- rep(0, 10)    #Start with 0 slopes
opt_cij_var1 <- rep(0.001, 10)  #Start with tiny variances in slopes
opt_base_haz <- lambda          #Start with original hazards


#---- Values to match ----
#The first values are inputs based on climbing to desired inc rate at age 70
dem_inc <- c(0.25, 0.5, 1, head(EURODEM_inc_rates$Total_AD_1000PY, -1))
survival_data <- female_life_netherlands$CP

#---- Objective Function ----
dem_inc_rate_match <- function(PARAMETERS, #cij_slope[j], cij_var1, lambda[j] 
                               small_batch_n, timepoint, dem_inc_rate, 
                               opt_cij_slopes, opt_base_haz){
  
  return(abs(sim_dem_inc - dem_inc_rate))
}
