#---- Package Loading, Options, and Seed ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("future.apply", "parallel", "here", "magrittr")

set.seed(20190311)

#---- Source Files ----
source(here("RScripts", "sex-dementia_sim_parA.R"))          #The parameter file
source(here("RScripts", "sex-dementia_sim_data_gen.R"))      #The data generation script
source(here("RScripts", "sex-dementia_sim_data_analysis.R")) #The data analysis script
source(here("RScripts", "misc_custom_functions.R"))          #Other functions needed

#---- Generating one cohort ----
data_gen() %>% saveRDS(here("Data", "test_sim_results_A_20190312"))

#---- Running the simulation in parallel----
start_time <- Sys.time()
plan(multiprocess, workers = (detectCores() - 2))
sim_results <- future_replicate(runs, sex_dem_sim())

#Mean results
mean_results_mat <- data.frame(matrix(NA, nrow = num_tests, 
                                      ncol = nrow(sim_results))) %>%
  set_colnames(c("num_obs", "num_females", "num_males", "p_alive", 
                 "p_alive_females", "p_alive_males", "mortality_HRs(F:M)",
                 "at_risk_females", "at_risk_males", "inc_cases_females", 
                 "inc_cases_males", "dem_cases_females", "dem_cases_males", 
                 "PY_females", "PY_males", "dem_prob_females", "dem_prob_males", 
                 "sim_rates","sim_rates_females", "sim_rates_males", 
                 "IRRs(F:M)", "dem_HRs(F:M)"))

mean_results_mat[1, "num_obs"] <- mean(unlist(sim_results["num_obs", ]))
mean_results_mat[1, "num_females"] <- mean(unlist(sim_results["num_females", ]))
mean_results_mat[1, "num_males"] <- mean(unlist(sim_results["num_males", ]))

mean_results_mat[, "p_alive"] <- sim_results["p_alive", ] %>% 
  Reduce(`+`, .)/runs
mean_results_mat[, "p_alive_females"] <- sim_results["p_alive_females", ] %>% 
  Reduce(`+`, .)/runs
mean_results_mat[, "p_alive_males"] <- sim_results["p_alive_males", ] %>% 
  Reduce(`+`, .)/runs
mean_results_mat[, "mortality_HRs(F:M)"] <- 
  exp(sim_results["simulated_mortality_logHRs", ] %>% Reduce(`+`, .)/runs)

mean_results_mat[, "at_risk_females"] <- 
  sim_results["at_risk_females", ] %>% Reduce(`+`, .)/runs
mean_results_mat[, "at_risk_males"] <- 
  sim_results["at_risk_males", ] %>% 
  Reduce(`+`, .)/runs

mean_results_mat[, "inc_cases_females"] <- 
  sim_results["inc_cases_females", ] %>% Reduce(`+`, .)/runs
mean_results_mat[, "inc_cases_males"] <- 
  sim_results["inc_cases_males", ] %>% 
  Reduce(`+`, .)/runs

mean_results_mat[, "dem_cases_females"] <- 
  sim_results["dem_cases_female", ] %>% Reduce(`+`, .)/runs
mean_results_mat[, "dem_cases_males"] <- 
  sim_results["dem_cases_male", ] %>% 
  Reduce(`+`, .)/runs

mean_results_mat[, "PY_females"] <- 
  sim_results["PY_females", ] %>% Reduce(`+`, .)/runs
mean_results_mat[, "PY_males"] <- 
  sim_results["PY_males", ] %>% Reduce(`+`, .)/runs

mean_results_mat[, "dem_prob_females"] <- 
  sim_results["dem_prob_females", ] %>% Reduce(`+`, .)/runs
mean_results_mat[, "dem_prob_males"] <- 
  sim_results["dem_prob_males", ] %>% 
  Reduce(`+`, .)/runs

mean_results_mat[, "sim_rates"] <- 
  sim_results["sim_rates", ] %>% Reduce(`+`, .)/runs
mean_results_mat[, "sim_rates_females"] <- 
  sim_results["sim_rates_females", ] %>% Reduce(`+`, .)/runs
mean_results_mat[, "sim_rates_males"] <- 
  sim_results["sim_rates_males", ] %>% Reduce(`+`, .)/runs

mean_results_mat[, "IRRs(F:M)"] <- 
  exp(sim_results["logIRRs", ] %>% Reduce(`+`, .)/runs) 

mean_results_mat[, "dem_HRs(F:M)"] <- 
  exp(sim_results["simulated_dementia_logHRs", ] %>% Reduce(`+`, .)/runs)

write_csv(mean_results_mat,
            here("Results", "Scenario_A_no_bias", 
                 paste0("sim_results_", gsub("-", "", Sys.Date()), ".csv")))

Sys.time() - start_time











