#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("ggplot2", "tidyverse")

#Suppress warnings
options(warn = -1)

#---- Specify source file ----
par_file <- "sex-dementia_sim_parA.R" #This syntax is used for file naming later
source(par_file)
source("sex-dementia_sim_script.R")
source("life_table2014.R")
source("dementia_incidence2000-2013.R")

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
  #Defining mean Cij plot data for females
  female_meanCij <- mean_Cij_check %>% filter(sex == 1) %>% 
    dplyr::select(-sex) %>% t() %>% as.data.frame() %>% 
    mutate("female" = V1) %>% dplyr::select(-V1) %>% 
    cbind(., "t" = visit_times) %>% melt(., id.vars = "t")

  #Defining mean Cij plot data for males
  male_meanCij <- mean_Cij_check %>% filter(sex == 0) %>% 
    dplyr::select(-sex) %>% t() %>% as.data.frame() %>% mutate("male" = V1) %>% 
    dplyr::select(-V1) %>% cbind(., "t" = visit_times) %>% 
    melt(., id.vars = "t")

  #Checking the mean slopes by sex
  male_slopes <- data_frame("b0a" = (male_meanCij[5, "value"] - 
                                       male_meanCij[1, "value"])/
                              (male_meanCij[5, "t"] - male_meanCij[1, "t"]), 
                            "b0b" = (male_meanCij[8, "value"] - 
                                       male_meanCij[5, "value"])/
                              (male_meanCij[8, "t"] - male_meanCij[5, "t"]), 
                            "b0c" = (male_meanCij[11, "value"] - 
                                       male_meanCij[8, "value"])/
                              (male_meanCij[11, "t"] - male_meanCij[8, "t"]))
  female_slopes <- data_frame("b0a" = (female_meanCij[5, "value"] - 
                                       female_meanCij[1, "value"])/
                              (female_meanCij[5, "t"] - female_meanCij[1, "t"]), 
                            "b0b" = (female_meanCij[8, "value"] - 
                                       female_meanCij[5, "value"])/
                              (female_meanCij[8, "t"] - female_meanCij[5, "t"]), 
                            "b0c" = (female_meanCij[11, "value"] - 
                                       female_meanCij[8, "value"])/
                              (female_meanCij[11, "t"] - 
                                 female_meanCij[8, "t"]))
  
  #Combine all plot data into one dataframe (includes random sample of Cij)
  
  samp_Cij <- sample_n(obs_check, 10) %>% 
    dplyr::select(dput(Cij_varnames[-1])) %>% t() %>%
    cbind(., "t" = visit_times) %>% as.data.frame() %>% 
    melt(., id.vars = "t") %>% rbind(., female_meanCij) %>% 
    rbind(., male_meanCij)

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
    labs(y = "Cognitive function", 
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
  ggsave(filename = paste("mean_Cij_samp_plot", str_extract(par_file, ".\\."), 
                          "jpeg", sep = ""), width = 10, height = 7, 
         plot = Cij_plot_samp)

  ggsave(filename = paste("mean_Cij_plot", str_extract(par_file, ".\\."), 
                          "jpeg", sep = ""), width = 10, height = 7, 
         plot = Cij_plot)

#---- Comparing with life-table data ----
#Based on 2014 life table found in 
#National Vital Statistics Reports, Vol. 66, No. 4, August 14, 2017 (pg 48-49)
#Make sure the appropriate return values are "turned on" in the simulation script
sample_sim <- replicate(1, sex_dem_sim()) 
all_obs <- sample_sim["obs", ] %>% do.call(rbind, .)

#Conditional probability of survival at each timepoint by sex
all_alive <- all_obs[, dput(deathij_varnames)] %>% 
  map_dbl(.f = ~length(.) - sum(.)) %>% cond_prob() 
all_alive[1] <- 1 - sum(all_obs$death01)/nrow(all_obs)

female_alive <- all_obs[, c("sex", dput(deathij_varnames))] %>% 
  filter(sex == 0) %>% dplyr::select(-sex) %>%
  map_dbl(.f = ~length(.) - sum(.)) %>% cond_prob()
female_alive[1] <- filter(all_obs, sex == 0) %>% dplyr::select(death01) %>%
  map_dbl(.f = ~1 - sum(.)/length(.)) 
  
male_alive <- all_obs[, c("sex", dput(deathij_varnames))] %>% 
  filter(sex == 1) %>% dplyr::select(-sex) %>%
  map_dbl(.f = ~length(.) - sum(.)) %>% cond_prob()
male_alive[1] <- filter(all_obs, sex == 1) %>% dplyr::select(death01) %>%
  map_dbl(.f = ~1 - sum(.)/length(.)) 

#---- Comparing with dementia incidence data ----
#Make sure the appropriate return values are "turned on" in the simulation script
sample_sim <- replicate(35, sex_dem_sim())
dem_1000py <- sample_sim["obs", ] %>% do.call(rbind, .) %>% 
  dplyr::select(dput(dem_varnames)) %>% 
  map_dbl(.f = 
            ~1000*sum(., na.rm = TRUE)/(int_time*(length(.) - sum(is.na(.)))))



