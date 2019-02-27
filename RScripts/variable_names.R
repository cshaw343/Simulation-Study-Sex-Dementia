#---- Package loading, options, seed ----
if (!require("pacman")) 
install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse")

#---- Generating variable names for dataset ----
variable_names <- tibble("exo_var" = c("id", "sex", "U", 
                                       rep(NA, (num_tests + 1) - 3)),
                         "timepoints" = seq(from = 0, to = num_tests, by = 1),
                         #have to go to timepoint 11 just to make enough rows
                         "timepoints_nobase" = seq(from = 1, 
                                                   to = (num_tests + 1), 
                                                   by = 1),
                         "start_ages" = seq(50, 95, by = 5), 
                         "end_ages" = seq(55, 100, by = 5),
                         "c50" = rep("c50", num_tests + 1),
                         "age" = rep("age", num_tests + 1), 
                         "eps" = rep("eps", num_tests + 1),
                         "delta" = rep("delta", num_tests + 1),
                         "dem" = rep("dem", num_tests + 1), 
                         "Ci" = rep("Ci", num_tests + 1), 
                         "std_Ci" = rep("std_Ci", num_tests + 1),
                         "cij_slope" = rep("cij_slope", num_tests + 1),
                         "rij" = rep("rij", num_tests + 1), 
                         "death" = rep("death", num_tests + 1), 
                         "survtime" = rep("survtime", num_tests + 1), 
                         "contributed" = rep("contributed", num_tests + 1)) %>% 
  #Interval timepoints
  unite("interval_times", c(timepoints, timepoints_nobase), sep = "-", 
        remove = FALSE) %>%
  #Interval timepoints + baseline
  mutate("interval_times_base" = c(0, head(interval_times, -1))) %>%
  #Interval ages
  unite("interval_ages", c(start_ages, end_ages), sep = "-", 
        remove = FALSE) %>%
  #Age labels
  unite("age_varnames", c(age, timepoints), sep = "", remove = FALSE) %>% 
  #Centered age labels
  unite("agec_varnames", c(age_varnames, c50), sep = "_", remove = FALSE) %>% 
  #Random noise labels
  unite("eps_varnames", c(eps, timepoints), sep = "", remove = FALSE) %>% 
  #Dementia labels
  unite("dem_varnames", c(dem, interval_times_base), sep = "", 
        remove = FALSE) %>% 
  #Cij labels
  unite("Cij_varnames", c(Ci, timepoints), sep = "", remove = FALSE) %>%
  #Standardized Cij labels
  unite("std_Cij_varnames", c(std_Ci, timepoints), sep = "", remove = FALSE) %>%
  #Interval slope labels
  unite("cij_slopeij_varnames", c(cij_slope, interval_times), sep = "", 
        remove = FALSE) %>%
  #Uniform survival noise labels
  unite("rij_varnames", c(rij, interval_times), sep = "", remove = FALSE) %>%
  #Death indicator labels
  unite("deathij_varnames", c(death, interval_times), sep = "", 
        remove = FALSE) %>%
  #Survival time labels
  unite("Sij_varnames", c(survtime, interval_times), sep = "", 
        remove = FALSE) %>%
  #Contributed time labels
  unite("contributed_varnames", c(contributed, interval_times), sep = "", 
        remove = FALSE) %>% 
  dplyr::select("exo_var", "cij_slopeij_varnames", "rij_varnames", 
                "deathij_varnames", "Sij_varnames", "contributed_varnames", 
                "age_varnames", "agec_varnames", "eps_varnames", "Cij_varnames")

#NAs for those intervals that don't exist in the data set
variable_names[nrow(variable_names), 
               c("cij_slopeij_varnames", "rij_varnames", "deathij_varnames", 
                 "Sij_varnames", "contributed_varnames")] <- NA

column_names <- c(na.omit(variable_names$exo_var), 
                  na.omit(variable_names$age_varnames), 
                  na.omit(variable_names$agec_varnames))

