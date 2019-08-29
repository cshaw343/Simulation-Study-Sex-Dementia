#---- Package Loading, Options, and Seed ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("parallel", "here", "magrittr")

set.seed(20190624)

#---- Source Files ----
source(
  here("RScripts", 
       "sex-dementia_sim_parA_onedemcut_uniform_timetodem_nodemkill_UonInt.R"))  #The parameter file
source(here("RScripts", "sex-dementia_sim_data_gen_quad.R"))      #The data generation script
source(here("RScripts", "sex-dementia_sim_data_analysis.R")) #The data analysis script
source(here("RScripts", "misc_custom_functions.R"))          #Other functions needed

#---- Generating one cohort ----
data_gen(500000) %>%
  saveRDS(here(
    "Data", 
    "dataset_C_onedemcut_uniform_timetodem_nodemkill_UonInt_500000_20190824"))


#---- Running the simulation in parallel----
# #Function to run simulation in batches
# batch_runs <- function(){
#   plan(multiprocess, 
#        workers = (floor(0.5*detectCores())), 
#        gc = TRUE)
#   sim_results <- future_replicate(2, sex_dem_sim(num_obs)) %>% t() %>% 
#     as.data.frame()
#   
#   return(sim_results)
# 
#   future:::ClusterRegistry("stop")
# }

# #---- Test Code ----
# if(runs%%5 != 0){
#   stop("Number of runs must be a multiple of 5.")
# }
#  
gc()
Start <- Sys.time()

#sim_results <- replicate(runs/5, batch_runs())

cl <- makeCluster(0.5*detectCores())  

# get library support needed to run the code
clusterEvalQ(cl, {
  if (!require("pacman")) 
    install.packages("pacman", repos='http://cran.us.r-project.org')
  
  p_load("here", "magrittr")
  
  source(
    here("RScripts", 
         "sex-dementia_sim_parA_onedemcut_uniform_timetodem_nodemkill_UonInt.R"))  #The parameter file
  source(here("RScripts", "sex-dementia_sim_data_gen_quad.R"))      #The data generation script
  source(here("RScripts", "sex-dementia_sim_data_analysis.R")) #The data analysis script
  source(here("RScripts", "misc_custom_functions.R"))  
})

sim_results <- parSapply(cl, 1:runs, function(i) {sex_dem_sim(1000)}) %>% 
  t() %>% as.data.frame() %>% 
  write_csv(
    here("Results", "Scenario_A", 
         "one_demcut_uniform_timetodem_nodemkill_UonInt_100_quad_test.csv"))

#stop the cluster
stopCluster(cl)

Sys.time() - Start

# #---- Output column names ----
# output_column_names <- rownames(sim_results)
# 
# #---- Create results matrix ----
# results_matrix <- 
#   as.data.frame(matrix(NA, ncol = length(output_column_names), nrow = runs)) %>% 
#   set_colnames(output_column_names)
# 
# for(i in output_column_names){
#   results_matrix[, i] <- unlist(sim_results[i, ])
# }

# write_csv(results_matrix, 
#           here("Results", "Scenario_A_no_bias", 
#                "one_demcut_nodemkill_1000_20190814.csv"))

