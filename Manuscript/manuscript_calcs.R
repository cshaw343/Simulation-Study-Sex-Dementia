#---- Seed setting + package loading ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "here")

#---- read in sim results ----
sim_A <- read_csv(here("Results", "Scenario_A_no_bias", 
                       "sim_A_1000_20191125.csv"))
sim_B1 <- read_csv(here("Results", "Scenario_B", "sim_B1_1000_20191128.csv"))
sim_B2 <- read_csv(here("Results", "Scenario_B", "sim_B2_1000_20191126.csv"))
sim_C1 <- read_csv(here("Results", "Scenario_C", "sim_C1_1000_20191129.csv"))
sim_C2 <- read_csv(here("Results", "Scenario_C", "sim_C2_1000_20191126.csv"))

#---- individual samples ----
set.seed(20200226)
#Create one sample of size 100,000 for each simulation scenario
sample_A <- data_gen(num_obs = 100000)

source(here("RScripts","scenario_B1_pars.R"))
sample_B1 <- data_gen(num_obs = 100000)

source(here("RScripts","scenario_B2_pars.R"))
sample_B2 <- data_gen(num_obs = 100000)

source(here("RScripts","scenario_C1_pars.R"))
sample_C1 <- data_gen(num_obs = 100000)

source(here("RScripts","scenario_C2_pars.R"))
sample_C2 <- data_gen(num_obs = 100000)

#---- median survival ----
#Can just do this for simA b/c all others similarly calibrated

survtime_dists <- function(){
  data <- data_gen(num_obs = 100000)
  female <- data %>% filter(female == 1)
  male <- data %>% filter(female == 0)
  
  return(c(median(female$survtime), median(male$survtime)))
}

meds <- replicate(100, survtime_dists())
female_med <- median(meds[1, ])
female_IQR <- summary(meds[1, ])
male_med <- median(meds[2, ])
male_IQR <- summary(meds[2, ])

sample_A_women <- sample_A %>% filter(female == 1)
sample_A_men <- sample_A %>% filter(female == 0)

summary(sample_A_women$survtime)
summary(sample_A_men$survtime)

#---- cumulative survival ----
#Can just do this for simA b/c all others similarly calibrated
1 - mean(sim_A$p_alive_95)
sd(sim_A$p_alive_95)



