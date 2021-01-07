#---- Package loading, options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "here", "magrittr", "MASS", "matrixcalc", "Matrix", 
       "survival")

options(scipen = 999) #Standard Notation
options(digits = 6)   #Round to 6 decimal places
options(warn = -1)    #Suppress warnings

#---- Source files ----
#Choose the simulation scenario that you want to input here
source(here("RScripts", "scenario_A_pars.R"))
source(here("RScripts", "data_gen.R"))
source(here("RScripts", "dem_inc_ACT.R"))
source(here("RScripts", "life_table_calcs.R"))

#---- Set values ----
#You set the values to test in this section
opt_linear_term <- 0.0425
opt_quadratic_term <- -0.002875

#Fixed values ------------------
opt_baseline_var <- 0.05
opt_ci_cov02 <- 0
opt_ci_cov12 <- 0
#-------------------------------

opt_ci_cov01 <- -0.00009
opt_linear_var <- 0.001 
opt_quadratic_var <- 0.000009 

opt_dem_cut <- -6.5
 
#---- Check PD matrix ----
#PD is necessary for calculation in data_gen.R script
quad_coeff_cov_test <- matrix(c(opt_baseline_var, opt_ci_cov01, opt_ci_cov02,
                                opt_ci_cov01, opt_linear_var, opt_ci_cov12,
                                opt_ci_cov02, opt_ci_cov12,
                                opt_quadratic_var),
                              nrow = 3, byrow = TRUE)

if(!is.positive.definite(quad_coeff_cov_test)){
  quad_coeff_cov_test <-
    as.matrix(nearPD(quad_coeff_cov_test, corr = FALSE, keepDiag = FALSE)$mat)

  opt_baseline_var <- quad_coeff_cov_test[1, 1]
  opt_linear_var <- quad_coeff_cov_test[2, 2]
  opt_quadratic_var <- quad_coeff_cov_test[3, 3]
  opt_cij_cov01 <- quad_coeff_cov_test[1, 2]
  opt_cij_cov02 <- quad_coeff_cov_test[1, 3]
  opt_cij_cov12 <- quad_coeff_cov_test[2, 3]

  message("covariance parameters have changed:")

  quad_coeff_cov_test
}

#Assigns the test values to these variables so original code for data_gen can
#be used
b10 <- opt_linear_term
b20 <- opt_quadratic_term
cij_var0 <- opt_baseline_var
cij_var1 <- opt_linear_var
cij_var2 <- opt_quadratic_var
cij_cov01 <- opt_ci_cov01
cij_cov02 <- opt_ci_cov02
cij_cov12 <- opt_ci_cov12
dem_cut <- opt_dem_cut

#---- Generate the data ----
num_obs = 500000
obs <- data_gen(num_obs)

#---- Data by sex ----
male_data <- obs %>% filter(female == 0)
female_data <- obs %>% filter(female == 1)

#---- Compute incidence rates ----
all_inc_cases <- colSums(obs[, variable_names$dem_varnames])
all_PY <- colSums(obs[, variable_names$contributed_varnames[1:9]])
all_sim_inc_rates <- all_inc_cases[-1]/all_PY*1000
names(all_sim_inc_rates) <- variable_names$interval_ages[1:num_tests]

male_inc_cases <- colSums(male_data[, variable_names$dem_varnames])
male_PY <- colSums(male_data[, variable_names$contributed_varnames[1:9]])
male_sim_inc_rates <- male_inc_cases[-1]/male_PY*1000
names(male_sim_inc_rates) <- variable_names$interval_ages[1:num_tests]

female_inc_cases <- colSums(female_data[, variable_names$dem_varnames])
female_PY <- colSums(female_data[, variable_names$contributed_varnames[1:9]])
female_sim_inc_rates <- female_inc_cases[-1]/male_PY*1000
names(female_sim_inc_rates) <- variable_names$interval_ages[1:num_tests]

#---- Compute survival data ----

#---- Cohort size ----
num_obs <- nrow(obs)
num_males <- nrow(male_data)
num_females <- nrow(female_data)

#---- Cumulative survival probabilities ----
p_alive <- obs %>% 
  dplyr::select(variable_names$deathij_varnames[1:num_tests]) %>% 
  map_dbl(~sum(. == 0, na.rm = TRUE))/num_obs

#---- Cumulative survival probabilities by sex ----
p_alive_males <- male_data %>%
  dplyr::select(variable_names$deathij_varnames[1:num_tests]) %>% 
  map_dbl(~sum(. == 0, na.rm = TRUE))/num_males
 
p_alive_females <- female_data %>%
  dplyr::select(variable_names$deathij_varnames[1:num_tests]) %>% 
  map_dbl(~sum(. == 0, na.rm = TRUE))/num_females

#---- Conditional survival probabilities ----
cond_p_alive <- vector(length = 
                         length(na.omit(variable_names$deathij_varnames)))
alive <- obs %>% 
  dplyr::select(variable_names$deathij_varnames[1:num_tests])

for(i in 1:length(cond_p_alive)){
  if(i == 1){
    cond_p_alive[i] = 
      sum(alive[, variable_names$deathij_varnames[i]] == 0)/nrow(alive)
  } else{
    alive <- alive[alive[, variable_names$deathij_varnames[i - 1]] == 0, ]
    cond_p_alive[i] = 
      sum(alive[, variable_names$deathij_varnames[i]] == 0)/nrow(alive)
  }
}

#---- Conditional survival probabilities by sex ----
cond_p_alive_males <- vector(length = 
                         length(na.omit(variable_names$deathij_varnames)))
alive_males <- male_data %>% 
  dplyr::select(variable_names$deathij_varnames[1:num_tests])

for(i in 1:length(cond_p_alive_males)){
  if(i == 1){
    cond_p_alive_males[i] = 
      sum(alive_males[, variable_names$deathij_varnames[i]] == 0)/
      nrow(alive_males)
  } else{
    alive_males <- 
      alive_males[alive_males[, variable_names$deathij_varnames[i - 1]] == 0, ]
    cond_p_alive_males[i] = 
      sum(alive_males[, variable_names$deathij_varnames[i]] == 0)/
      nrow(alive_males)
  }
}

cond_p_alive_females <- vector(length = 
                         length(na.omit(variable_names$deathij_varnames)))
alive_females <- female_data %>% 
  dplyr::select(variable_names$deathij_varnames[1:num_tests])

for(i in 1:length(cond_p_alive)){
  if(i == 1){
    cond_p_alive_females[i] = 
      sum(alive_females[, variable_names$deathij_varnames[i]] == 0)/
      nrow(alive_females)
  } else{
    alive_females <- 
      alive_females[
        alive_females[, variable_names$deathij_varnames[i - 1]] == 0, ]
    cond_p_alive_females[i] = 
      sum(alive_females[, variable_names$deathij_varnames[i]] == 0)/
      nrow(alive_females)
  }
}

#---- Visualize Cij data ----
Cij_check <- obs %>% 
  dplyr::select("female", variable_names$Ci_varnames)

#Getting Cij by sex
female_Cij <- Cij_check %>% filter(female == 1) 
male_Cij <- Cij_check %>% filter(female == 0)

#Getting mean Cij by sex
female_meanCij<- female_Cij %>% colMeans(., na.rm = TRUE)
male_mean_Cij <- male_Cij %>% colMeans(., na.rm = TRUE)

female_mean_plot <- tibble("Age" = c(seq(50, 95, by = 5)), 
                           "variable" = "Women", 
                           "value" = female_meanCij[-1])

male_mean_plot <- tibble("Age" = seq(50, 95, by = 5), 
                         "variable" = "Men", 
                         "value" = male_mean_Cij[-1])

Cij_plot_data <- obs %>% 
  dplyr::select("id", "female", variable_names$Ci_varnames, "survtime", 
                "last_Ci") %>% sample_n(100) %>% 
  mutate("Age" = survtime + 50) %>% 
  dplyr::select(-c(survtime, female)) %>% 
  set_colnames(c("id", variable_names$Ci_varnames, "Cij", "Age")) 


#Plot data
samp_Cij <- Cij_plot_data[, c("id", variable_names$Ci_varnames)] %>%
  set_colnames(c("id", seq(50, 95, by = 5))) %>% 
  gather(as.character(seq(50, 95, by = 5)), key = "Age", value = "Cij") %>% 
  mutate_at("Age", as.numeric) %>% 
  rbind(., Cij_plot_data[, c("id", "Age", "Cij")]) %>% 
  set_colnames(c("variable", "Age", "value")) %>% 
  rbind(., female_mean_plot) %>% 
  rbind(., male_mean_plot)

#---- Plot a sample of Cij ----
#Creating a plot with random sample in the background
ggplot(samp_Cij, aes(Age, value)) + 
  geom_line(data = 
              subset(samp_Cij, variable != "Women" & variable != "Men"), 
            aes(group = variable), color = "gray") +
  geom_line(data = subset(samp_Cij, variable == "Women"),
            aes(color = variable), size = 1.25) +
  geom_line(data = subset(samp_Cij, variable == "Men"),
            aes(color = variable), size = 1.25, alpha = 0.6) +
  #geom_point(data = ) +
  geom_hline(yintercept = dem_cut, size = 1.25) + 
  labs(y = "Cognitive Function", 
       x = "Age", 
       color = "Mean Cognitive \n Function") + 
  scale_x_continuous(breaks = seq(50, 95, 5)) + 
  ggtitle("Mean Cognitive Trajectories") +
  theme_minimal() +  
  #theme(text = element_text(size = 28)) + 
  coord_cartesian(ylim = c(-10, 10)) + 
  geom_hline(yintercept = dem_cut) + 
  guides(color = guide_legend(reverse = TRUE))

#---- Checking values ----
head(ACT_inc_rates$Male_All_Dementia_1000PY)
scenario_A_rates <- c(7.046, 8.864, 21.173, 45.525, 71.635, 96.946)
scenario_A_rates
male_sim_inc_rates

simulated_mortality_logHRs <- vector(length = num_tests)

for(i in 1:length(simulated_mortality_logHRs)){
  survtime_name <- paste0("survtime", i - 1, "-", i)
  death_indicator_name <- paste0("death", i - 1, "-", i)
  cox_model <- coxph(Surv(obs[, survtime_name], 
                          obs[, death_indicator_name]) ~ obs$female)
  simulated_mortality_logHRs[i] <- cox_model$coefficients
}

Hratio_US[-1, ]
exp(simulated_mortality_logHRs)

#conditional survival
male_life_US$CP[-1]
cond_p_alive_males

female_life_US$CP[-1]
cond_p_alive_females



