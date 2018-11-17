#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("ggplot2", "tidyverse", "reshape", "magrittr")

options(warn = -1)    #Suppress warnings
options(scipen = 999) #Standard notation

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

#---- Cij Distributions ----
#Check the distribution of Cij across timepoints
dem_cut <- -4.75

Cij_data <- obs %>% dplyr::select(variable_names$Cij_varnames)
Cij_indicators <- matrix(nrow = nrow(Cij_data), ncol = ncol(Cij_data))
for(i in 1:ncol(Cij_data)){
  if(i == 1){
    indicators <- rep(0, nrow(Cij_indicators))
  } else {
    indicators <- (Cij_data[, i] <= dem_cut & Cij_data[, (i-1)] > dem_cut)*1
    indicators[is.na(indicators)] <- 0
  }
  Cij_indicators[, i] <- indicators
}
indicators <- unlist(as.data.frame(Cij_indicators))

hists <- Cij_data %>% gather() %>% mutate_at("key", as.factor) 
hists$indicators <- indicators
hists$key <- fct_relevel(hists$key, "Ci10", after = 10)


z <- data.frame(variable = levels(hists$key), 
                cutoff = rep(-4.75, length(levels(hists$key))))
  
Cij_dists <- hists %>% ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free") + 
  geom_vline(data = z, aes(xintercept = cutoff), color = "black", size = 1) + 
  geom_histogram(fill = "gray", alpha = 0.6) + 
  geom_histogram(data = hists %>% filter(indicators == 1), fill = "dodgerblue2", 
                 alpha = 0.6) + 
  xlim(-20, 10) + theme_minimal()

#Saving plot output
lgd <- format(Sys.time(), "%Y_%m_%d_%H:%M:%S") #format the time/date for file creation
ggsave(filename = paste("Plots/Cij_dists_dem_inc_parA_", lgd, ".jpeg", 
                        sep = ""), width = 10, height = 7, plot = Cij_dists)

#---- Mean Fij Plot ----
Fij_check <- obs %>% dplyr::select(sex, contains("Fi"))

#Getting mean Fij by sex
female_meanFij<- Fij_check %>% filter(sex == 0) %>% colMeans()
male_mean_Fij <- Fij_check %>% filter(sex == 1) %>% colMeans()

female_meanFij_plot <- tibble("t" = visit_times, "variable" = "female", 
                              "value" = head(female_meanFij[-1], 11))

male_meanFij_plot <- tibble("t" = visit_times, "variable" = "male", 
                            "value" = head(male_mean_Fij[-1], 11))

Fij_plot_data <- rbind(female_meanFij_plot, male_meanFij_plot)

Fij_plot <- ggplot(Fij_plot_data, aes(t, value)) + 
  geom_line(data = subset(Fij_plot_data, variable == "female"), 
            aes(color = variable), size = 1.25) +
  geom_line(data = subset(Fij_plot_data, variable == "male"), 
            aes(color = variable), size = 1.25, alpha = 0.6) + 
  labs(y = "Functional Ability", 
       x = "Visit Time", 
       color = "") + ylim(-20, 10) + 
  theme_minimal()

#Saving plot output
lgd <- format(Sys.time(), "%Y_%m_%d_%H:%M:%S") #format the time/date for file creation
ggsave(filename = paste("Plots/mean_Fij_parA_", lgd, ".jpeg", 
                        sep = ""), width = 10, height = 7, plot = Fij_plot)

#---- Comparing with life-table data ----
#Based on 2014 life table found in 
#National Vital Statistics Reports, Vol. 66, No. 4, August 14, 2017 (pg 48-49)

  #---- Look at survival data ----
  sim_data <- 
    results_mat[, c("sex", head(variable_names$deathij_varnames,-1))] 
  
  #Want to look at unconditional and conditional probabilities of survival
  #Unconditional Survival
  death_intervals <- colnames(sim_data)[-1]
  all_death_counts <- sim_data %>% dplyr::select(-one_of("sex")) %>% 
    colSums() 
  all_survival_counts <- nrow(sim_data) - all_death_counts
  all_cp_survival <- cond_prob(all_survival_counts)
  all_cp_survival[1] <- all_survival_counts[1]/(nrow(sim_data))
  
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
  
  #Visualize by plots
  plot_data <- tibble("age" = seq(55, 100, by = 5), 
                      "pub_all_survival" = life$CP[-c(1, 2)], 
                      "pub_female_survival" = female_life$CP[-c(1, 2)], 
                      "pub_male_survival" = male_life$CP[-c(1, 2)], 
                      "cohort_all_survival" = all_cp_survival, 
                      "cohort_female_survival" = female_cp_survival, 
                      "cohort_male_survival" = male_cp_survival) %>% 
    gather(key = "data_type", value = "cp", colnames(plot_data)[-1])
  
  all_survival_plot <- ggplot(plot_data, aes(age, cp)) + 
    geom_line(data = subset(plot_data, data_type == "cohort_all_survival"), 
              aes(color = "Simulated"), size = 1.25, color = "dodgerblue2", 
              alpha = 0.6) + 
    geom_line(data = subset(plot_data, data_type == "pub_all_survival"), 
              aes(color = "Published"), size = 1.25, 
              color = "gray", linetype = "longdash") + ylim(0, 1) + 
    labs(y = "Conditional Probability of Survival", x = "Age", 
         color = "") + ggtitle("Entire Cohort") + theme_minimal()
  
  survival_female_plot <- ggplot(plot_data, aes(age, cp)) + 
    geom_line(data = subset(plot_data, data_type == "cohort_female_survival"), 
              aes(color = "Simulated"), size = 1.25, color = "dodgerblue2", 
              alpha = 0.6) + 
    geom_line(data = subset(plot_data, data_type == "pub_female_survival"), 
              aes(color = "Published"), size = 1.25, 
              color = "gray", linetype = "longdash") + ylim(0, 1) + 
    labs(y = "Conditional Probability of Survival (Females)", x = "Age", 
         color = "") + ggtitle("Females") + theme_minimal()
  
  survival_male_plot <- ggplot(plot_data, aes(age, cp)) + 
    geom_line(data = subset(plot_data, data_type == "cohort_male_survival"), 
              aes(color = "Simulated"), size = 1.25, color = "dodgerblue2", 
              alpha = 0.6) + 
    geom_line(data = subset(plot_data, data_type == "pub_male_survival"), 
              aes(color = "Published"), size = 1.25, 
              color = "gray", linetype = "longdash") + ylim(0, 1) + 
    labs(y = "Conditional Probability of Survival (Males)", x = "Age", 
         color = "") + ggtitle("Males") + theme_minimal()
  
  #Saving plot output
  lgd <- format(Sys.time(), "%Y_%m_%d_%H:%M:%S") #format the time/date for file creation
  ggsave(filename = paste("Plots/all_survival_plot_parA_", lgd, ".jpeg", 
                          sep = ""), width = 10, height = 7, 
         plot = all_survival_plot)
  
  ggsave(filename = paste("Plots/female_survival_plot_parA_", lgd, ".jpeg", 
                          sep = ""), width = 10, height = 7, 
         plot = survival_female_plot)
  
  ggsave(filename = paste("Plots/male_survival_plot_parA_", lgd, ".jpeg", 
                          sep = ""), width = 10, height = 7, 
         plot = survival_male_plot)
  
  
#---- Comparing with dementia incidence data ----
#Make sure the appropriate return values are "turned on" in the simulation script
sample_sim <- replicate(35, sex_dem_sim())
dem_1000py <- sample_sim["obs", ] %>% do.call(rbind, .) %>% 
  dplyr::select(dput(dem_varnames)) %>% 
  map_dbl(.f = 
            ~1000*sum(., na.rm = TRUE)/(int_time*(length(.) - sum(is.na(.)))))



