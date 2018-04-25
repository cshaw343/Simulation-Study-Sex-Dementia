#---- Specify source file ----
source("sex-demensia_sim_script.R")

#---- The simulation function ----
sex_dem_sim_check <- function(){
  #---- Generating IDs, sex, U ----
  obs <- tibble("id" = seq(from = 1, to = num_obs, by = 1), 
                "sex" = rbinom(num_obs, size = 1, prob = psex), 
                "U" = rnorm(num_obs, mean = 0, sd = 1))
  
  #---- Generating age data ----
  #Creating ages at each timepoint j
  ages <- as_tibble(matrix(NA, nrow = num_obs, ncol = length(age_varnames)))
  for(j in 1:length(age_varnames)){
    if(j == 1){
      ages[, j] = seq(from = 1, to = num_obs, by = 1) #Creates column of ids
    } else if(j == 2){
      ages[, j] = age0 #Creates column of baseline ages
    } else ages[, j] = ages[, (j-1)] + int_time #Creates ages at following timepoints
  }
  colnames(ages) <- age_varnames
  
  #---- Generating centered age data ----
  #Creating centered ages at each timepoint j
  c_ages <- as_tibble(ages - mean(age0)) %>% 
    mutate("id" = seq(from = 1, to = num_obs, by = 1)) #Creates column of ids
  colnames(c_ages) <- agec_varnames
  
  #---- Generating "true" cognitive function Cij ----
  #Cij = b00 + z0i + bo1*sexi + b02*age0i + b03*Ui + (b10 + z1i + b11*sexi + 
  #b12*age0i + b13*Ui)*timej + epsilonij
  
  #---- Generating random terms for slope and intercept ----
  #Generate random terms for each individual
  slope_int_noise <- as_tibble(mvrnorm(n = num_obs, mu = rep(0, 2), 
                                       Sigma = slope_int_cov)) %>% 
    cbind("id" = seq(from = 1, to = num_obs, by = 1), .) #Creates column of ids
  colnames(slope_int_noise) <- c("id", "z0i", "z1i")
  
  #---- Generating noise term (unexplained variance in Cij) for each visit ----
  sd_eps <- sqrt(var3)
  eps <- as_tibble(replicate(num_tests + 1, 
                             rnorm(n = num_obs, mean = 0, sd = sd_eps))) %>% 
    cbind("id" = seq(from = 1, to = num_obs, by = 1), .) #Creates column of ids
  colnames(eps) <- eps_varnames
  
  #---- Creating complete matrix of observation data ----
  obs <- left_join(obs, ages, by = "id") %>% left_join(c_ages, by = "id") %>%
    left_join(slope_int_noise, by = "id") %>% left_join(eps, by = "id")
  
  #---- Calculating Cij for each individual ----
  #Store Cij values
  Cij <- as.data.frame(cog_func(obs)$Cij) %>% 
    cbind("id" = seq(from = 1, to = num_obs, by = 1), .) #Creating column of ids
  colnames(Cij) <- Cij_varnames
  
  #Store slope values per interval per individual
  slopeij <- as.data.frame(cog_func(obs)$slopes) %>% 
    cbind("id" = seq(from = 1, to = num_obs, by = 1), .) #Creating column of ids
  colnames(slopeij) <- slopeij_varnames
  
  #---- Calculating mean Cij by sex ----
  mean_Cij <- Cij %>% mutate("sex" = obs$sex) %>% 
    mutate_at("sex", as.factor) %>% group_by(sex) %>% 
    dplyr::select(-id) %>% summarise_all(mean)
  
  #---- Generate survival time for each person ----
  #Individual hazard functions
  #h(tij|x) = lambda*exp(g1*sexi + g2*ageij + g3*Ui + g4*sexi + 
  #g5*slopeij + g6Cij)
  #See Additional notes in README file
  
  #---- Generating uniform random variables per interval for Sij ----
  USij <- as_tibble(replicate(num_tests, 
                              runif(num_obs, min = 0, max = 1))) %>%
    cbind("id" = seq(from = 1, to = num_obs, by = 1), .) #Creating column of ids
  colnames(USij) <- USij_varnames
  
  #---- Merging Cij, slopeij, and USij with observation data ----
  #Used as input for survival function
  obs <- left_join(obs, Cij, by = "id") %>% left_join(slopeij, by = "id") %>% 
    left_join(USij, by = "id")
  
  #---- Calculating Sij for each individual ----
  #Store Sij values
  Sij <- as.data.frame(survival(obs)) 
  
  #---- Calculating death data for each individual ----
  #Compute death indicator for each interval
  deathij <- (Sij < int_time)*1 
  for(i in 1:nrow(deathij)){
    death <- min(which(deathij[i, ] == 1))
    if(is.finite(death)){
      deathij[i, death:ncol(deathij)] = 1 #Changes death indicators to 1 after death
    }
  }
  
  #Compute study death indicators
  study_death <- (rowSums(deathij) > 0)*1
  
  #Compute overall survival times
  survtime <- vector(length = num_obs)
  survtime[which(study_death == 0)] = num_tests*int_time
  for(i in 1:length(survtime)){
    if(survtime[i] == 0){
      death_int <- min(which(deathij[i, ] == 1)) 
      survtime[i] = int_time*(death_int - 1) + Sij[i, death_int]
    } 
  }
  
  #Computing age at death
  age_death <- age0 + survtime
  
  #Labelling datasets and appending IDs
  Sij <- cbind("id" = seq(from = 1, to = num_obs, by = 1), Sij) #Creating column of ids
  colnames(Sij) <- Sij_varnames
  
  deathij <- cbind("id" = seq(from = 1, to = num_obs, by = 1), deathij) #Creating column of ids
  colnames(deathij) <- deathij_varnames
  
  survtime <- cbind("id" = seq(from = 1, to = num_obs, by = 1), survtime) #Creating column of ids
  colnames(survtime) <- c("id", "survtime")
  
  age_death <- cbind("id" = seq(from = 1, to = num_obs, by = 1), age_death) #Creating column of ids
  colnames(age_death) <- c("id", "age_death")
  
  #---- Censor Cij based on death data ----
  for(i in 1:num_obs){
    death_int <- (min(which(deathij[i, ] == 1)) - 1) 
    if(is.finite(death_int)){
      Cs <- vector(length = num_tests)
      for(j in death_int:num_tests){
        Cs[j] <- paste("Ci", j, sep = "")
      }
      Cs <- Cs[Cs != "FALSE"]
      obs[i, dput(Cs)] <- NA
    }
  }
  return(list("obs" = obs, "mean_Cij" = mean_Cij))
}


#---- Checking the simulated data----
#Storing the results of the simulation
sim_check <- sex_dem_sim_check()
obs_check <- as_tibble(sim_check$obs)
mean_Cij_check <- as_tibble(sim_check$mean_Cij)

#Check means: proportion of males, U
means <- obs_check %>% summarise_at(c("sex", "U"), mean)

#---- Checking by plots ----
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




