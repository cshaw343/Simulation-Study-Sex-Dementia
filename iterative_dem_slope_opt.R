#---- Package loading, options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("here")

#---- Source files ----
source(here("RScripts", "sex-dementia_sim_parA.R"))
source(here("RScripts", "dementia_incidence_EURODEM_pooled.R"))

#---- Pre-allocation ----
opt_cij_slopes <- rep(0, 10)  #Start with 0 slopes
opt_base_haz <- lambda        #Start with original hazards

#---- Values to match ----
#The first values are inputs based on climbing to desired inc rate at age 70
dem_inc <- c(0.25, 0.5, 1, head(EURODEM_inc_rates$Total_AD_1000PY, -1))

#---- Objective Function ----

