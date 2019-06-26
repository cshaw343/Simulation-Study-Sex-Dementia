#---- Package loading, options, seed ----
if (!require("pacman")) 
install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse")

#---- Generating variable names for dataset (5 year) ----
variable_names <- tibble("exo_var" = c("id", "sex", "female", "U", 
                                       rep(NA, (num_tests + 1) - 4)),
                         "timepoints" = seq(from = 0, to = num_tests, by = 1),
                         #have to go to timepoint 1 + num_tests just to make enough rows (because of baseline)
                         "timepoints_nobase" = seq(from = 1, 
                                                   to = (num_tests + 1), 
                                                   by = 1),
                         "start_ages" = seq(50, 95, by = 5), 
                         "end_ages" = seq(55, 100, by = 5),
                         "c50" = rep("c50", num_tests + 1),
                         "age" = rep("age", num_tests + 1), 
                         "eps" = rep("eps", num_tests + 1),
                         "z0" = rep("z0", num_tests + 1),
                         "z1" = rep("z1", num_tests + 1),
                         "i" = rep("i", num_tests + 1),
                         "delta" = rep("delta", num_tests + 1),
                         "dem" = rep("dem", num_tests + 1), 
                         "Ci" = rep("Ci", num_tests + 1), 
                         "std_Ci" = rep("std_Ci", num_tests + 1),
                         "cij_slope" = rep("cij_slope", num_tests + 1),
                         "rij" = rep("rij", num_tests + 1), 
                         "death" = rep("death", num_tests + 1), 
                         "survtime" = rep("survtime", num_tests + 1), 
                         "contributed" = rep("contributed", num_tests + 1), 
                         "p_alive" = rep("p_alive", num_tests + 1), 
                         "females" = rep("females", num_tests + 1), 
                         "males" = rep("males", num_tests + 1), 
                         "mortality_logHR" = 
                           rep("mortality_logHR", num_tests + 1), 
                         "at_risk" = rep("at_risk", num_tests + 1), 
                         "dem_inc_rate" = 
                           rep("dem_inc_rate", num_tests + 1), 
                         "inc_cases" = rep("inc_cases", num_tests + 1), 
                         "PY" = rep("PY", num_tests + 1), 
                         "logIRR" = rep("logIRR", num_tests + 1), 
                         "dementia_logHR" = 
                           rep("dementia_logHR", num_tests + 1), 
                         "SE" = rep("SE", num_tests + 1),
                         "CI_95_Upper" = rep("CI_95_Upper", num_tests + 1),
                         "CI_95_Lower" = rep("CI_95_Lower", num_tests + 1),
                         "CI_95_Coverage" = 
                           rep("CI_95_Coverage", num_tests + 1),
                         "dem_cases" = rep("dem_cases", num_tests + 1), 
                         "prop_dem" = rep("prop_dem", num_tests + 1), 
                         "mean_U" = rep("mean_U", num_tests + 1)) %>% 
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
  #Slope-int noise labels
  unite("int_noise", c(z0, timepoints), sep = "_", remove = FALSE) %>%
  unite("int_noise_names", c(int_noise, i), sep = "", remove = FALSE) %>%
  unite("slope_noise", c(z1, timepoints), sep = "_", remove = FALSE) %>%
  unite("slope_noise_names", c(slope_noise, i), sep = "", remove = FALSE) %>%
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
  #For return values
  unite("p_alive_varnames", c(p_alive, end_ages), sep = "_", remove = FALSE) %>%
  unite("p_alive_females_varnames", c(p_alive, females, end_ages), sep = "_", 
        remove = FALSE) %>%
  unite("p_alive_males_varnames", c(p_alive, males, end_ages), sep = "_", 
        remove = FALSE) %>% 
  unite("mortality_logHR_varnames", c(mortality_logHR, interval_ages), sep = "_", 
        remove = FALSE) %>% 
  unite("at_risk_females_varnames", c(at_risk, females, start_ages), sep = "_", 
        remove = FALSE) %>% 
  unite("at_risk_males_varnames", c(at_risk, males, start_ages), sep = "_", 
        remove = FALSE) %>% 
  unite("dem_inc_rate_varnames", c(dem_inc_rate, interval_ages), sep = "_", 
        remove = FALSE) %>% 
  unite("dem_inc_rate_females_varnames", 
        c(dem_inc_rate, females, interval_ages), sep = "_", remove = FALSE) %>% 
  unite("dem_inc_rate_males_varnames", 
        c(dem_inc_rate, males, interval_ages), sep = "_", remove = FALSE) %>% 
  unite("inc_cases_females_varnames", 
        c(inc_cases, females, interval_ages), sep = "_", remove = FALSE) %>% 
  unite("inc_cases_males_varnames", 
        c(inc_cases, males, interval_ages), sep = "_", remove = FALSE) %>% 
  unite("PY_females_varnames", 
        c(PY, females, interval_ages), sep = "_", remove = FALSE) %>% 
  unite("PY_males_varnames", 
        c(PY, males, interval_ages), sep = "_", remove = FALSE) %>% 
  unite("logIRR_varnames", 
        c(logIRR, interval_ages), sep = "_", remove = FALSE) %>% 
  unite("logIRR_SE_varnames", 
        c(logIRR, SE, interval_ages), sep = "_", remove = FALSE) %>% 
  unite("logIRR_95CI_Upper_varnames", 
        c(logIRR, CI_95_Upper, interval_ages), sep = "_", remove = FALSE) %>% 
  unite("logIRR_95CI_Lower_varnames", 
        c(logIRR, CI_95_Lower, interval_ages), sep = "_", remove = FALSE) %>% 
  unite("logIRR_95CI_Coverage_varnames", 
        c(logIRR, CI_95_Coverage, interval_ages), sep = "_", remove = FALSE) %>% 
  unite("dementia_logHR_varnames", c(dementia_logHR, interval_ages), sep = "_", 
        remove = FALSE) %>% 
  unite("dementia_logHR_SE_varnames", c(dementia_logHR, SE, interval_ages), 
        sep = "_", remove = FALSE) %>%
  unite("dementia_logHR_95CI_Upper_varnames", 
        c(dementia_logHR, CI_95_Upper, interval_ages), 
        sep = "_", remove = FALSE) %>%
  unite("dementia_logHR_95CI_Lower_varnames", 
        c(dementia_logHR, CI_95_Lower, interval_ages), 
        sep = "_", remove = FALSE) %>%
  unite("dementia_logHR_95CI_Coverage_varnames", 
        c(dementia_logHR, CI_95_Coverage, interval_ages), 
        sep = "_", remove = FALSE) %>%
  unite("dem_cases_females_varnames", 
        c(dem_cases, females, interval_ages), sep = "_", remove = FALSE) %>%
  unite("dem_cases_males_varnames", 
        c(dem_cases, males, interval_ages), sep = "_", remove = FALSE) %>%
  unite("prop_dem_females_varnames", 
        c(prop_dem, females, end_ages), sep = "_", remove = FALSE) %>%
  unite("prop_dem_males_varnames", 
        c(prop_dem, males, end_ages), sep = "_", remove = FALSE) %>%
  unite("mean_U_at_risk_females_varnames", 
        c(mean_U, at_risk, females, start_ages), sep = "_", remove = FALSE) %>%
  unite("mean_U_at_risk_males_varnames", 
        c(mean_U, at_risk, males, start_ages), sep = "_", remove = FALSE) %>%
  dplyr::select("exo_var", "int_noise_names", "slope_noise_names", 
                "cij_slopeij_varnames", "rij_varnames", 
                "deathij_varnames", "Sij_varnames", "contributed_varnames", 
                "age_varnames", "agec_varnames", "eps_varnames", "Cij_varnames", 
                "dem_varnames", "interval_ages", "p_alive_varnames", 
                "p_alive_females_varnames", "p_alive_males_varnames", 
                "mortality_logHR_varnames", "at_risk_females_varnames", 
                "at_risk_males_varnames", "dem_inc_rate_varnames", 
                "dem_inc_rate_females_varnames", "dem_inc_rate_males_varnames", 
                "inc_cases_females_varnames", "inc_cases_males_varnames", 
                "PY_females_varnames", "PY_males_varnames", "logIRR_varnames", 
                "logIRR_SE_varnames", "logIRR_95CI_Upper_varnames", 
                "logIRR_95CI_Lower_varnames", "logIRR_95CI_Coverage_varnames",
                "dementia_logHR_varnames", "dementia_logHR_SE_varnames",
                "dementia_logHR_95CI_Upper_varnames", 
                "dementia_logHR_95CI_Lower_varnames", 
                "dementia_logHR_95CI_Coverage_varnames",
                "dem_cases_females_varnames", 
                "dem_cases_males_varnames", "prop_dem_females_varnames", 
                "prop_dem_males_varnames", "mean_U_at_risk_females_varnames", 
                "mean_U_at_risk_males_varnames")

#---- Generating variable names for dataset (1 year) ----
variable_names_1year <- tibble("dem" = rep("dem", num_tests*5),
                               "females" = rep("females", num_tests*5), 
                               "males" = rep("males", num_tests*5),
                               "inc_cases" = rep("inc_cases", num_tests*5), 
                               "PY" = rep("PY", num_tests*5), 
                               "contributed" = rep("contributed", num_tests*5), 
                               "logIRR" = rep("logIRR", num_tests*5),
                               "SE" = rep("SE", num_tests*5),
                               "CI_95_Upper" = rep("CI_95_Upper", num_tests*5),
                               "CI_95_Lower" = rep("CI_95_Lower", num_tests*5),
                               "CI_95_Coverage" = 
                                 rep("CI_95_Coverage", num_tests*5),
                               "start_ages" = seq(50, 94, by = 1), 
                               "end_ages" = seq(51, 95, by = 1)) %>%
  unite("interval_ages", c(start_ages, end_ages), sep = "-", 
        remove = FALSE) %>%
  unite("contributed_varnames", c(contributed, interval_ages), sep = "", 
        remove = FALSE) %>%
  unite("dem_varnames", c(dem, interval_ages), sep = "", remove = FALSE) %>%
  unite("inc_cases_females_varnames", 
        c(inc_cases, females, interval_ages), sep = "_", remove = FALSE) %>% 
  unite("inc_cases_males_varnames", 
        c(inc_cases, males, interval_ages), sep = "_", remove = FALSE) %>% 
  unite("PY_females_varnames", 
        c(PY, females, interval_ages), sep = "_", remove = FALSE) %>% 
  unite("PY_males_varnames", 
        c(PY, males, interval_ages), sep = "_", remove = FALSE) %>% 
  unite("logIRR_varnames", 
        c(logIRR, interval_ages), sep = "_", remove = FALSE) %>% 
  unite("logIRR_SE_varnames", 
        c(logIRR, SE, interval_ages), sep = "_", remove = FALSE) %>% 
  unite("logIRR_95CI_Upper_varnames", 
        c(logIRR, CI_95_Upper, interval_ages), sep = "_", remove = FALSE) %>% 
  unite("logIRR_95CI_Lower_varnames", 
        c(logIRR, CI_95_Lower, interval_ages), sep = "_", remove = FALSE) %>% 
  unite("logIRR_95CI_Coverage_varnames", 
        c(logIRR, CI_95_Coverage, interval_ages), sep = "_", remove = FALSE) %>%
  dplyr::select("contributed_varnames", "dem_varnames", 
                "inc_cases_females_varnames", "inc_cases_males_varnames", 
                "PY_females_varnames", "PY_males_varnames", "logIRR_varnames", 
                "logIRR_SE_varnames", "logIRR_95CI_Upper_varnames", 
                "logIRR_95CI_Lower_varnames", "logIRR_95CI_Coverage_varnames")

#---- NAs for those intervals that don't exist in the data sets ----
variable_names[nrow(variable_names), 
               c("cij_slopeij_varnames", "rij_varnames", "deathij_varnames", 
                 "Sij_varnames", "contributed_varnames", "interval_ages", 
                 "p_alive_varnames", "p_alive_females_varnames", 
                 "p_alive_males_varnames", "mortality_logHR_varnames", 
                 "at_risk_females_varnames", "at_risk_males_varnames", 
                 "dem_inc_rate_varnames", "dem_inc_rate_females_varnames", 
                 "dem_inc_rate_males_varnames", "inc_cases_females_varnames", 
                 "inc_cases_males_varnames", "PY_females_varnames", 
                 "PY_males_varnames", "logIRR_varnames", "logIRR_SE_varnames", 
                 "logIRR_95CI_Upper_varnames", 
                 "logIRR_95CI_Lower_varnames", "logIRR_95CI_Coverage_varnames",
                 "dementia_logHR_varnames", "dementia_logHR_SE_varnames", 
                 "dementia_logHR_95CI_Upper_varnames", 
                 "dementia_logHR_95CI_Lower_varnames", 
                 "dementia_logHR_95CI_Coverage_varnames",
                 "dem_cases_females_varnames", 
                 "dem_cases_males_varnames", "prop_dem_females_varnames", 
                 "prop_dem_males_varnames", "mean_U_at_risk_females_varnames", 
                 "mean_U_at_risk_males_varnames")] <- NA

#---- Define column names for simulation dataset ----
column_names <- c(na.omit(variable_names$exo_var), variable_names$age_varnames, 
                  variable_names$agec_varnames, variable_names$int_noise_names, 
                  variable_names$slope_noise_names, variable_names$eps_varnames, 
                  variable_names$Cij_varnames, 
                  variable_names$cij_slopeij_varnames[1:num_tests], 
                  variable_names$rij_varnames[1:num_tests], 
                  variable_names$Sij_varnames[1:num_tests], 
                  "death0", variable_names$deathij_varnames[1:num_tests], 
                  "study_death", "survtime", "age_death", 
                  variable_names$dem_varnames, "dem_wave", "dem", "timetodem", 
                  "ageatdem", "dem_death", "timetodem_death", "ageatdem_death", 
                  "dem_alive", variable_names$contributed_varnames[1:num_tests], 
                  variable_names_1year$contributed_varnames, 
                  variable_names_1year$dem_varnames)


