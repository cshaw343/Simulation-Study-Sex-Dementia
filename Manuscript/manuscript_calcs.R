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

#---- median survival ----
#Can just do this for simA b/c all others similarly calibrated


#---- cumulative survival ----
#Can just do this for simA b/c all others similarly calibrated
1 - mean(sim_A$p_alive_95)
sd(sim_A$p_alive_95)

#---- U stuff ----
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