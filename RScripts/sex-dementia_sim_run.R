#---- Package Loading, Options, and Seed ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("future.apply", "parallel", "here")

set.seed(20190311)

#---- Source Files ----
source(here("RScripts", "sex-dementia_sim_parA.R"))          #The parameter file
source(here("RScripts", "sex-dementia_sim_data_gen.R"))      #The data generation script
source(here("RScripts", "sex-dementia_sim_data_analysis.R")) #The data analysis script
source(here("RScripts", "misc_custom_functions.R"))          #Other functions needed

#---- Generating one cohort ----
data_gen() %>% saveRDS(here("Data", "test_sim_results_A_20190311"))

#---- Running the simulation in parallel----
start_time <- Sys.time()
runs = 100
plan(multiprocess, workers = (detectCores() - 2))
sim_results <- future_replicate(runs, sex_dem_sim())
Sys.time() - start_time

#Means of results
mean_num_obs <- mean(unlist(sim_results["num_obs", ]))
mean_num_females <- mean(unlist(sim_results["num_females", ]))
mean_num_males <- mean(unlist(sim_results["num_males", ]))

mean_at_risk_by_sex <- sim_results["at_risk_by_sex", ] %>% 
  Reduce(`+`, .)/runs
mean_inc_cases_by_sex <- sim_results["inc_cases_by_sex", ] %>% 
  Reduce(`+`, .)/runs
mean_dem_cases_by_sex <- sim_results["dem_cases_by_sex", ] %>% 
  Reduce(`+`, .)/runs
mean_PY_by_sex <- sim_results["PY_by_sex", ] %>% 
  Reduce(`+`, .)/runs
mean_dem_prob_by_sex <- sim_results["dem_prob_by_sex", ] %>% 
  Reduce(`+`, .)/runs

mean_sim_rates <- sim_results["sim_rates", ] %>% Reduce(`+`, .)/runs

#Female:Male
mean_IRRs <- exp(sim_results["logIRRs", ] %>% Reduce(`+`, .)/runs) 

Sys.time() - start_time








