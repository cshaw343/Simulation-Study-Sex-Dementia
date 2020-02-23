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

#---- Generating one cohort ----
#Use this to get one cohort for cohort specific checks and plots 
#(like the U plots)
#Set the desired number of people in the data_gen function
#Set the desired file name in the here function where it says "filename" 
data_gen(num_obs = 1000) %>%
  saveRDS(here("Data", "filename"))

#---- Running the full simulation ----
#Simulation settings
runs = 2         #Number of simulation runs
num_obs <- 1000   #Size of each simulated cohort

gc()                #Clear the environment of unecessary junk
Start <- Sys.time()

start <- Sys.time()
plan(multiprocess, workers = 0.5*availableCores())

cl <- makeCluster(0.5*detectCores())  

# get library support needed to run the code
clusterEvalQ(cl, {
  if (!require("pacman")) 
    install.packages("pacman", repos='http://cran.us.r-project.org')
  
  p_load("here", "magrittr")
  
  source(
    here("RScripts", "quadratic_model", 
         "sex-dementia_sim_parC_onedemcut_nodemkill_male_AllDem_quad.R"))  #The parameter file
  source(here("RScripts", "quadratic_model", "variable_names_quad.R"))
  source(here("RScripts", "quadratic_model", 
              "sex-dementia_sim_data_gen_quad.R"))      #The data generation script
  source(here("RScripts", "sex-dementia_sim_data_analysis.R")) #The data analysis script
})

sim_results <- parSapply(cl, 1:runs, function(i) {sex_dem_sim(num_obs)}) %>% 
  t() %>% as.data.frame() 

write_csv(sim_results, 
          here("Results", "quadratic_model", "Scenario_C", 
               paste0("sim_C_male_AllDem_", runs, "_20200212.csv")))

#stop the cluster
stopCluster(cl)

Sys.time() - Start
