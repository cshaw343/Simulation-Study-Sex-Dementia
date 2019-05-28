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
data_gen() %>% saveRDS(here("Data", "dataset_A_onedemcut_20190528"))

#---- Running the simulation in parallel----
#Function to run simulation in batches
batch_runs <- function(x){
  plan(multiprocess, 
       workers = (floor(0.5*detectCores())), 
       gc = TRUE)
  sim_results <- future_replicate(5, sex_dem_sim()) %>% t() %>% 
    as.data.frame()
  
  return(sim_results)

  future:::ClusterRegistry("stop")
}

#---- Test Code ----
if(runs%%5 != 0){
  stop("Number of runs must be a multiple of 5.")
}
 
gc()
Start <- Sys.time()

sim_results <- replicate(runs/5, batch_runs())

Sys.time() - Start

#---- Output column names ----
output_column_names <- rownames(sim_results)

#---- Create results matrix ----
results_matrix <- 
  as.data.frame(matrix(NA, ncol = length(output_column_names), nrow = runs)) %>% 
  set_colnames(output_column_names)

for(i in output_column_names){
  results_matrix[, i] <- unlist(sim_results[i, ])
}

write_csv(results_matrix, here("Results", "Scenario_A_no_bias", 
                               "one_demcut_100_20190528.csv"))
  






