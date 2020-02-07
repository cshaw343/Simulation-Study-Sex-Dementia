#---- Package Loading, Options, and Seed ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("here", "tidyverse")

#---- Source scripts ----
source(here("RScripts", "dementia_incidence_ACT.R"))
source(here("RScripts", "quadratic_model", "variable_names_quad.R"))

#---- Load results data ----
results_A <- read_csv(here(
  "Results", "quadratic_model", "Scenario_A_no_bias", 
  "sim_A_male_AllDem_1000_20191125.csv")) %>% 
  results_65_plus()

results_B1 <- read_csv(here(
  "Results", "quadratic_model", "Scenario_B", 
  "sim_B_male_AllDem_1000_20191128.csv")) %>% 
  results_65_plus()

results_B2 <- read_csv(here(
  "Results", "quadratic_model", "Scenario_B", 
  "sim_B_highUonSurv_male_AllDem_1000_20191126.csv")) %>% 
  results_65_plus()

results_C1 <- read_csv(here(
  "Results", "quadratic_model", "Scenario_C", 
  "sim_C_male_AllDem_1000_20191129.csv")) %>% 
  results_65_plus()

results_C2 <- read_csv(here(
  "Results", "quadratic_model", "Scenario_C", 
  "sim_C_highUonSurv_male_AllDem_1000_20191126.csv")) %>% 
  results_65_plus()

#---- Calculate mean and sd of results ----
mean_results_A <- results_A %>% colMeans()
mean_results_B1 <- results_B1 %>% colMeans()
mean_results_B2 <- results_B2 %>% colMeans()
mean_results_C1 <- results_C1 %>% colMeans()
mean_results_C2 <- results_C2 %>% colMeans()

#---- Figure 2 ----
figure2_data <- 
  data.frame("Ages" = c(80, 85, 90), 
             "ACT" = ACT_inc_rates$`All_Dem_IRR_F:M`[
               which(ACT_inc_rates$Low_Age == 80):
                 which(ACT_inc_rates$Low_Age == 90)], 
             "A" = )






+ scale_x_continuous(name = "Age range", breaks = c(70,80,90,100), labels=c("65-74", "75-84","85-94","95-104"))+"
IRR_table <- 
  tibble("ACT" = " ", 
         "IRR (Women:Men)_ACT" = c(rep("NA", 3), 
                                   round(ACT_inc_rates$`All_Dem_IRR_F:M`, 3)),
         "Simulation Scenario A" = " ", 
         "IRR (Women:Men)" = exp(
           mean_sim_results[na.omit(variable_names$logIRR_varnames)])) %>% 
  t() %>% 
  set_colnames(age_labels$formatted_age_intervals) %>% 
  set_rownames(c("ACT", "IRR (Women:Men)", "Simulation Scenario A", 
                 "IRR (Women:Men)"))
