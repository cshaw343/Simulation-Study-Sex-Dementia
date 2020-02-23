#---- Package Loading, Options, and Seed ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("future.apply", "here", "magrittr")

set.seed(20200113)

#---- Source Files ----
#Specify the parameter file
source(here("RScripts", "scenario_A_pars.R"))  #The parameter file
source(here("RScripts", "var_names.R"))
source(here("RScripts", "data_gen.R"))         #The data generation script
source(here("RScripts", "data_analysis.R"))    #The data analysis script

#---- Set filenames for output ----
#Paths are taken care of with the here package
#Specify string for filename
one_cohort_output <- "one_cohort_A"
simulation_output <- "sim_run_A.csv"

#---- Generating one cohort ----
#Use this to get one cohort for cohort specific checks and plots 
#(like the U plots)
#Set the desired number of people in the data_gen function
#Set the desired file name in the here function where it says "filename" 
data_gen(num_obs = 1000) %>%
  saveRDS(here("Data", one_cohort_output))

#---- Running the full simulation ----
#Simulation settings
runs = 2              #Number of simulation runs
num_obs <- 100000     #Size of each simulated cohort

gc()                  #Clear the environment of unecessary junk
start <- Sys.time()   #Start timing the code

#Run in parallel using half of the available cores
plan(multiprocess, workers = 0.5*availableCores()) 

sim_results <- future_replicate(runs, sex_dem_sim(num_obs)) %>% t() %>% 
  as.data.frame() 

write_csv(sim_results, here("Results", simulation_output))

#stop the multiprocess-- do this or else things will be SO slow after
plan(sequential)

#How long did it take
Sys.time() - start
