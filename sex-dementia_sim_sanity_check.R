#---- Package Loading ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("ggplot2")

#---- Specify source file ----
par_file <- "sex-dementia_sim_parA.R" #This syntax is used for file naming later
source(par_file)
source("sex-dementia_sim_script.R")
source("life_table2014.R")

#---- Checking one simulated dataset----
#Storing the results of the simulation
sim_check <- sex_dem_sim()
obs_check <- as_tibble(sim_check$obs)
mean_Cij_check <- as_tibble(sim_check$mean_Cij)

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

death_check <- replicate(5, sex_dem_sim())

#Conditional probability of survival at each timepoint by sex
death_check_male <- obs[, c("sex", dput(deathij_varnames))] %>% 
  filter(sex == 1) %>% dplyr::select(-sex) %>% 
  map_dbl(.f = ~ length(.) - sum(.)) %>% cond_prob()
death_check_female <- obs[, c("sex", dput(deathij_varnames))] %>% 
  filter(sex == 0) %>% dplyr::select(-sex) %>% 
  map_dbl(.f = ~ length(.) - sum(.)) %>% cond_prob()





