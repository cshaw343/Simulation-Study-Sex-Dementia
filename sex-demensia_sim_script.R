#---- Package loading, options, seed ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "MASS", "reshape")

#Using standard notation (as opposed to scientific), rounded to three 
#decimal places
options(scipen = 999)
options(digits = 3)

set.seed(10789)

#---- Specify the parameter file ----
par_file <- "sex-demensia_sim_parA.R"
source(par_file)

#---- Generating variable names for each assessment timepoint ----
age_varnames <- c("id", "age0", vector(length = num_tests)) #Age labels
agec_varnames <- c("id", "age0_c50", vector(length = num_tests)) #Centered age labels
eps_varnames <- c("id", "eps0", vector(length = num_tests)) #Epsilon labels
Cij_varnames <- c("id", "Ci0", vector(length = num_tests)) #Cij labels
for(j in 1:num_tests){
  age_varnames[j + 2] = paste("age", j, sep = "")
  agec_varnames[j + 2] = paste(age_varnames[j + 2], "_c50", sep = "")
  eps_varnames[j + 2] = paste("eps", j, sep = "")
  Cij_varnames[j + 2] = paste("Ci", j, sep = "")
}

#---- Generating assessment timepoint data ----
visit_times <- seq(from = 0, to = int_time*num_tests, by = int_time)

#---- Model for Cognitive Function ----
cog_func_1 <- function(obs){
  knots = c(0, 20, 35)
  Cij <- vector(length = length(visit_times))
  for(j in 1:length(Cij)){
    t = visit_times[j]
    test_num = j - 1
    eps <- paste("eps", test_num, sep = "")
    if(t >= knots[1] & t < knots[2]){
      Cij[j] = b00 + obs[, "z0i"] + b01*obs[, "sex"] + b02*obs[, "age0_c50"] + 
        b03*obs[, "U"] + obs[, eps] + 
        (b10a + obs[, "z1i"] + b11*obs[, "sex"] + b12*obs[, "age0_c50"] + 
           b13*obs[, "U"])*t
    } else if (t >= knots[2] & t < knots[3]){
      Cij[j] = b00 + obs[, "z0i"] + b01*obs[, "sex"] + b02*obs[, "age0_c50"] + 
        b03*obs[, "U"] + (b10a - b10b)*knots[2] + obs[, eps] + 
        (b10b + obs[, "z1i"] + b11*obs[, "sex"] + b12*obs[, "age0_c50"] + 
           b13*obs[, "U"])*t
    } else {
      Cij[j] = b00 + obs[, "z0i"] + b01*obs[, "sex"] + b02*obs[, "age0_c50"] + 
        b03*obs[, "U"] + (b10a - b10b)*knots[2] + (b10b - b10c)*knots[3] + 
        obs[, eps] + 
        (b10c + obs[, "z1i"] + b11*obs[, "sex"] + b12*obs[, "age0_c50"] + 
           b13*obs[, "U"])*t
    }
  }
  return(Cij)
}

cog_func_2 <- function(obs){
  knots = c(0, 20, 35)
  Cij <- vector(length = length(visit_times))
  for(j in 1:length(Cij)){
    t = visit_times[j]
    test_num = j - 1
    eps <- paste("eps", test_num, sep = "")
    Cij[j] = b00 + obs[, "z0i"] + b01*obs[, "sex"] + b02*obs[, "age0_c50"] + 
      b03*obs[, "U"] + obs[, eps] + 
      (b10a - b10b)*knots[2]*(t >= knots[2]) + 
      (b10b - b10c)*knots[3]*(t >= knots[3]) + 
      (obs[, "z1i"] + b11*obs[, "sex"] + b12*obs[, "age0_c50"] + b13*obs[, "U"] + 
         b10a*(t >= knots[1] & t< knots[2]) + 
         b10b*(t >= knots[2] & t< knots[3]) +
         b10c*(t >= knots[3]))*t
  }
  return(Cij)
}

#---- Generate Covariance Matrix for random slope and intercept terms ----
slope_int_cov <- matrix(c(var0, cov, cov, var1), nrow = 2, byrow = TRUE)

#---- The simulation function ----
sex_dem_sim <- function(){
  #---- Generating IDs, sex, U ----
  obs <- tibble("id" = seq(from = 1, to = num_obs, by = 1), 
                "sex" = rbinom(num_obs, size = 1, prob = psex), 
                "U" = rnorm(num_obs, mean = 0, sd = 1))
  
  #---- Generating age data ----
  #Creating ages at each timepoint j
  ages <- as_tibble(matrix(NA, nrow = num_obs, ncol = length(age_varnames)))
  for(j in 1:length(age_varnames)){
    if(j == 1){
      ages[, i] = seq(from = 1, to = num_obs, by = 1) #Creates column of ids
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
    Cij <- as.data.frame(cog_func(obs)) %>% 
      cbind("id" = seq(from = 1, to = num_obs, by = 1), .) #Creating column of ids
    colnames(Cij) <- Cij_varnames
    
  #Function returns a list of 
      #(1) matrix of observations 
      #(2) matrix of Cij
  results_list <- list("Cij" = Cij, "obs" = obs) 
  return(results_list)
}

#---- Running the simulation----
#Storing the results of the simulation
sim_results <- replicate(2, sex_dem_sim())

#---- Find mean Cij by sex ----
obs <- as_tibble(sim_results$obs)
Cij <- as_tibble(sim_results$Cij) %>% mutate("sex" = obs$sex) %>% 
  mutate_at("sex", as.factor) 

mean_Cij <- Cij %>% group_by(sex) %>% summarise_all(mean)

#---- Creating plot data ----
#Defining mean Cij plot data for females
female_meanCij <- mean_Cij %>% filter(sex == 1) %>% dplyr::select(-id, -sex) %>% 
  t() %>% as.data.frame() %>% mutate("female" = V1) %>% dplyr::select(-V1) %>% 
  cbind(., "t" = visit_times) %>% melt(., id.vars = "t")

#Defining mean Cij plot data for males
male_meanCij <- mean_Cij %>% filter(sex == 0) %>% dplyr::select(-id, -sex) %>% 
  t() %>% as.data.frame() %>% mutate("male" = V1) %>% dplyr::select(-V1) %>% 
  cbind(., "t" = visit_times) %>% melt(., id.vars = "t")

#Combine all plot data into one dataframe (includes random sample of Cij)
samp_Cij <- sample_n(Cij, 10) %>% dplyr::select(-id, -sex) %>% t() %>%
  cbind(., "t" = visit_times) %>% as.data.frame() %>% 
  melt(., id.vars = "t") %>% rbind(., female_meanCij) %>% rbind(., male_meanCij)

#---- Plot a sample of Cij ----
#Creating a plot with random sample in the background
Cij_plot_samp <- ggplot(samp_Cij, aes(t, value)) + 
  geom_line(data = subset(samp_Cij, variable != "female" & variable != "male"), 
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

ggsave(filename = paste("mean_Cij_plot", str_extract(par_file, ".\\."), "jpeg", 
                        sep = ""), width = 10, height = 7, plot = Cij_plot)

