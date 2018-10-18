#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("parallel")

#---- Source Files ----
source("RScripts/sex-dementia_sim_parA.R")    #The parameter file
source("RScripts/variable_names.R")           #Creates all the variable names
source("RScripts/sex-dementia_sim_script.R")  #The simulation script
source("RScripts/misc_custom_functions.R")    #Other functions needed

#---- Running the simulation in parallel----
use_cores <- detectCores() - 1                 #Use one less than number of available cores
cl <- makeCluster(use_cores, type = "FORK")    #Make a cluster from your cores 
                                               #"FORKING" environment only works on Mac
runs = 100
#Storing the results of the simulation
start <- proc.time()  #This times the code
sim_results <- parSapply(cl, 1:runs, sex_dem_sim)
#stop the cluster
stopCluster(cl)
proc.time() - start

#---- Running simulation on one core ----
runs = 100
sim_results <- replicate(runs, sex_dem_sim())

#---- Converting results to usable format ----
results_mat <- matrix(nrow = runs*num_obs, ncol = nrow(sim_results))
for(r in 1:nrow(sim_results)){
  results_mat[, r] <- unlist(sim_results[r, ])
}

results_mat %<>% as.data.frame() %>% set_colnames(dimnames(sim_results)[[1]])

  #---- Look at survival data ----
  test_life_table <- 
  sim_results[c("sex", head(variable_names$deathij_varnames,-1)), ] 

  results_mat <- matrix(nrow = runs*num_obs, ncol = nrow(test_life_table))
  for(r in 1:nrow(test_life_table)){
    results_mat[, r] <- unlist(test_life_table[r, ])
  }
  colnames(results_mat) <- c("sex", head(variable_names$deathij_varnames, -1))
  results_mat %<>% as.data.frame()
  
  #Want to look at unconditional and conditional probabilities of survival
    #Unconditional Survival
    death_intervals <- colnames(results_mat)[-1]
    all_death_counts <- results_mat %>% dplyr::select(-one_of("sex")) %>% 
      colSums() 
    all_survival_counts <- num_obs*runs - death_counts
    all_cp_survival <- cond_prob(all_survival_counts)
    all_cp_survival[1] <- all_survival_counts[1]/(num_obs*runs)
  
    #Female survival
    num_females <- results_mat %>% filter(sex == 0) %>% nrow()
    female_death_counts <- results_mat %>% filter(sex == 0) %>% 
      dplyr::select(-one_of("sex")) %>% colSums()
    female_survival_counts <- num_females - female_death_counts
    female_cp_survival <- cond_prob(female_survival_counts)
    female_cp_survival[1] <- female_survival_counts[1]/(num_females)
    
    #Male survival
    num_males <- results_mat %>% filter(sex == 1) %>% nrow()
    male_death_counts <- results_mat %>% filter(sex == 1) %>% 
      dplyr::select(-one_of("sex")) %>% colSums()
    male_survival_counts <- num_males - male_death_counts
    male_cp_survival <- cond_prob(male_survival_counts)
    male_cp_survival[1] <- male_survival_counts[1]/(num_males)
  






