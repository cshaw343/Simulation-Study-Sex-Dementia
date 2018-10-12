#---- Package loading, options, seed ----
if (!require("pacman")) 
install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse")

#---- Generating variable names for dataset ----
variable_names <- tibble("timepoints" = seq(from = 0, to = num_tests, by = 1),
                         #have to go to timepoint 11 just to make enough rows
                         "timepoints_nobase" = seq(from = 1, 
                                                   to = (num_tests + 1), 
                                                   by = 1),
                         "c50" = rep("c50", num_tests + 1),
                         "age" = rep("age", num_tests + 1), 
                         "eps" = rep("eps", num_tests + 1), 
                         "dem" = rep("dem", num_tests + 1), 
                         "Ci" = rep("Ci", num_tests + 1), 
                         "mean" = rep("mean", num_tests + 1), 
                         "cij_slope" = rep("cij_slope", num_tests + 1),
                         "fij_slope" = rep("fij_slope", num_tests + 1),
                         "rij" = rep("rij", num_tests + 1), 
                         "death" = rep("death", num_tests + 1), 
                         "survtime" = rep("survtime", num_tests + 1)) %>% 
  #Interval timepoints
  unite("interval_times", c(timepoints, timepoints_nobase), sep = "-", 
        remove = FALSE) %>%
  #Age labels
  unite("age_varnames", c(age, timepoints), sep = "", remove = FALSE) %>% 
  #Centered age labels
  unite("agec_varnames", c(age_varnames, c50), sep = "_", remove = FALSE) %>% 
  #Epsilon labels
  unite("eps_varnames", c(eps, timepoints), sep = "", remove = FALSE) %>% 
  #Dementia labels
  unite("dem_varnames", c(dem, timepoints), sep = "", remove = FALSE) %>% 
  #Cij labels
  unite("Cij_varnames", c(Ci, timepoints), sep = "", remove = FALSE) %>%
  #Mean Cij labels
  unite("mean_Cij_varnames", c(mean, Cij_varnames), sep = "", 
        remove = FALSE) %>%
  #Interval slope labels
  unite("cij_slopeij_varnames", c(cij_slope, interval_times), sep = "", 
        remove = FALSE) %>%
  unite("fij_slopeij_varnames", c(fij_slope, interval_times), sep = "", 
        remove = FALSE) %>%
  #Uniform survival noise labels
  unite("rij_varnames", c(rij, interval_times), sep = "", remove = FALSE) %>%
  #Death indicator labels
  unite("deathij_varnames", c(death, interval_times), sep = "", 
        remove = FALSE) %>%
  #Survival time labels
  unite("Sij_varnames", c(survtime, interval_times), sep = "", remove = FALSE)
