#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("ggplot2", "tidyverse")

#Suppress warnings
options(warn = -1)

#---- Specify source file ----
par_file <- "RScripts/sex-dementia_sim_parA.R" #This syntax is used for file naming later
source(par_file)
source("RScripts/sex-dementia_sim_script.R")
source("RScripts/life_table2014.R")
source("RScripts/dementia_incidence2000-2013.R")
source("RScripts/sex-dementia_sim_run.R")

#---- Checking one simulated dataset----
#Storing the results of the simulation
sim_check <- sex_dem_sim()
obs_check <- as_tibble(sim_check$obs)

#---- Calculating mean Cij by sex ----
mean_Cij <- obs %>% mutate_at("sex", as.factor) %>% group_by(sex) %>% 
  dplyr::select(sex, Cij_varnames) %>% summarise_all(mean) %>%
  set_colnames(mean_Cij_varnames)

#Would I still need this?
#mean_Cij_check <- as_tibble(sim_check$mean_Cij)

#Check the distribution of Cij across timepoints
obs_check %>% dplyr::select(Cij_varnames) %>% gather() %>% 
  ggplot(aes(value)) + facet_wrap(~ key, scales = "free") +
  geom_histogram()

Cij_data <- obs_check %>% dplyr::select(Cij_varnames) 
Cij_summaries <- as_tibble(matrix(nrow = 3, ncol = length(Cij_varnames))) %>% 
  set_colnames(Cij_varnames) %>% set_rownames(c("Med", "Mean", "SD"))
Cij_summaries["Med", ] <- Cij_data %>% summarise_all(funs(median), na.rm = TRUE)
Cij_summaries["Mean", ] <- Cij_data %>% summarise_all(funs(mean), na.rm = TRUE)
Cij_summaries["SD", ] <- Cij_data %>% summarise_all(funs(sd), na.rm = TRUE)

cov_Cij <- cov(Cij_data, use = "complete.obs")
sd_mat <- diag(sqrt(diag(cov_Cij)))
sd_mat_inv <- solve(sd_mat)
corr_Cij <- sd_mat_inv%*%cov_Cij%*%sd_mat_inv

#Check means: proportion of males, U
means <- obs_check %>% summarise_at(c("sex", "U"), mean)

  #---- Creating plot data ----
  Cij_check <- obs %>% dplyr::select(sex, contains("Ci"))
  
  #Getting mean Cij by sex
  female_meanCij<- Cij_check %>% filter(sex == 0) %>% colMeans()
  male_mean_Cij <- Cij_check %>% filter(sex == 1) %>% colMeans()
  
  #Checking the mean slopes by sex
  slopes_check <- obs %>% dplyr::select(sex, contains("slope"))
  female_slopes_check <- slopes_check %>% filter(sex == 0) %>% colMeans()
  male_slopes_check <- slopes_check %>% filter(sex == 1) %>% colMeans()
  
  female_mean_plot <- tibble("t" = visit_times, "variable" = "female", 
                             "value" = female_meanCij[-1])
    
  male_mean_plot <- tibble("t" = visit_times, "variable" = "male", 
                           "value" = male_mean_Cij[-1])
    
    
  #Create all plot data (includes random sample of Cij)
  samp_Cij <- sample_n(mean_Cij_check, 10) %>% 
    dplyr::select(-sex) %>% t() %>%
    cbind(., "t" = visit_times) %>% as.data.frame() %>% 
    melt(., id.vars = "t") %>% rbind(., female_mean_plot) %>% 
    rbind(., male_mean_plot)

  #---- Plot a sample of Cij ----
  #Creating a plot with random sample in the background
  Cij_plot_samp <- ggplot(samp_Cij, aes(t, value)) + 
    geom_line(data = 
                subset(samp_Cij, variable != "female" & variable != "male"), 
              aes(group = variable), color = "gray") +
    geom_line(data = subset(samp_Cij, variable == "female"), 
              aes(color = variable), size = 1.25) + 
    geom_line(data = subset(samp_Cij, variable == "male"), 
              aes(color = variable), size = 1.25, alpha = 0.6) + 
    labs(y = "Cognitive Function", 
         x = "Visit Time", 
         color = "Mean Cognitive \n Function") + 
    theme_minimal()

  #Creating a plot without random sample in the background
  Cij_plot<- ggplot(samp_Cij, aes(t, value)) + 
    geom_line(data = subset(samp_Cij, variable == "female"), 
              aes(color = variable), size = 1.25) + 
    geom_line(data = subset(samp_Cij, variable == "male"), 
              aes(color = variable), size = 1.25, alpha = 0.6) + ylim(-15, 5) +
    labs(y = "Cognitive function", 
         x = "Visit Time", 
         color = "Mean Cognitive \n Function") + 
    theme_minimal()

  #Saving plot output
  lgd <- format(Sys.time(), "%Y_%m_%d_%H:%M:%S") #format the time/date for file creation
  ggsave(filename = paste("Plots/mean_Cij_samp_plot_parA_", lgd, ".jpeg", 
                          sep = ""), width = 10, height = 7, 
         plot = Cij_plot_samp)

  ggsave(filename = paste("Plots/mean_Cij_plot_parA_", lgd, ".jpeg", sep = ""), 
         width = 10, height = 7, 
         plot = Cij_plot)

#---- Comparing with life-table data ----
#Based on 2014 life table found in 
#National Vital Statistics Reports, Vol. 66, No. 4, August 14, 2017 (pg 48-49)

  #---- Look at survival data ----
  sim_data <- 
    results_mat[, c("sex", head(variable_names$deathij_varnames,-1))] 
  
  sim_data <- matrix(nrow = runs*num_obs, ncol = nrow(test_life_table))
  for(r in 1:nrow(test_life_table)){
    sim_data[, r] <- unlist(test_life_table[r, ])
  }
  colnames(sim_data) <- c("sex", head(variable_names$deathij_varnames, -1))
  sim_data %<>% as.data.frame()
  
  #Want to look at unconditional and conditional probabilities of survival
  #Unconditional Survival
  death_intervals <- colnames(sim_data)[-1]
  all_death_counts <- sim_data %>% dplyr::select(-one_of("sex")) %>% 
    colSums() 
  all_survival_counts <- num_obs*runs - death_counts
  all_cp_survival <- cond_prob(all_survival_counts)
  all_cp_survival[1] <- all_survival_counts[1]/(num_obs*runs)
  
  #Female survival
  num_females <- sim_data %>% filter(sex == 0) %>% nrow()
  female_death_counts <- sim_data %>% filter(sex == 0) %>% 
    dplyr::select(-one_of("sex")) %>% colSums()
  female_survival_counts <- num_females - female_death_counts
  female_cp_survival <- cond_prob(female_survival_counts)
  female_cp_survival[1] <- female_survival_counts[1]/(num_females)
  
  #Male survival
  num_males <- sim_data %>% filter(sex == 1) %>% nrow()
  male_death_counts <- sim_data %>% filter(sex == 1) %>% 
    dplyr::select(-one_of("sex")) %>% colSums()
  male_survival_counts <- num_males - male_death_counts
  male_cp_survival <- cond_prob(male_survival_counts)
  male_cp_survival[1] <- male_survival_counts[1]/(num_males)
  
#---- Comparing with dementia incidence data ----
#Make sure the appropriate return values are "turned on" in the simulation script
sample_sim <- replicate(35, sex_dem_sim())
dem_1000py <- sample_sim["obs", ] %>% do.call(rbind, .) %>% 
  dplyr::select(dput(dem_varnames)) %>% 
  map_dbl(.f = 
            ~1000*sum(., na.rm = TRUE)/(int_time*(length(.) - sum(is.na(.)))))



