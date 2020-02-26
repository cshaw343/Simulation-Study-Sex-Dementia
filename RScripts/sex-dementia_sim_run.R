#---- Seed setting + package loading ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("parallel")

set.seed(20200113)

#---- Set filenames for output ----
#Paths are taken care of with the here package
#Specify string for filename
one_cohort_output <- "dataset_A_500000_20200223"
simulation_output <- "sim_run_B1_20200226.csv"

#Simulation settings
runs = 1000         #Number of simulation runs
num_obs <- 100000   #Size of each simulated cohort

#---- Create the cluster for parallel computing ----
#We run these simulations in parallel because each run is incredibly 
#computationally intensive

cl <- makeCluster(0.5*detectCores())  

#Export variables in this environment to cluster nodes
clusterExport(cl = cl, list("runs", "num_obs"))

#Evaluate code on cluster nodes
clusterEvalQ(cl, {
  if (!require("pacman")) 
    install.packages("pacman", repos='http://cran.us.r-project.org')
  
  p_load("here", "magrittr", "tidyverse")
  
  #Specify the parameter file
  source(here("RScripts", "scenario_B1_pars.R"))  #The parameter file
  source(here("RScripts", "var_names.R"))
  source(here("RScripts", "data_gen.R"))         #The data generation script
  source(here("RScripts", "data_analysis.R"))    #The data analysis script
})

#---- Generating one cohort ----
#Use this to get one cohort for cohort specific checks and plots 
#(like the U plots)
#Set the desired number of people in the data_gen function

parSapply(cl, 1:1, function(i) {data_gen(num_obs = 500000)}) %>% 
  saveRDS(here("Data", one_cohort_output))

#---- Running the full simulation ----
#Clear the environment of unecessary junk
gc()
#Start timing the code
start <- Sys.time()

sim_results <- parSapply(cl, 1:runs, function(i) {sex_dem_sim(num_obs)}) %>% 
  t() %>% as.data.frame() 

#stop the cluster-- do this or else things will be SO slow after
stopCluster(cl)

#How long did it take
Sys.time() - start
