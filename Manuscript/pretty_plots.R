#---- Package Loading, Options, and Seed ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("here", "tidyverse", "latex2exp", "magrittr")

#---- Source scripts ----
source(here("RScripts", "dementia_incidence_ACT.R"))
source(here("RScripts", "variable_names_quad.R"))
source(here("RScripts", "sex-dementia_sim_data_gen_quad.R"))

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
#Defining the column names that I want to pull
IRR_ages <- c("logIRR_80-85", "logIRR_85-90", "logIRR_90-95")
IRR_upper <- c("logIRR_CI_95_Upper_80-85", "logIRR_CI_95_Upper_85-90", 
               "logIRR_CI_95_Upper_90-95")
IRR_lower <- c("logIRR_CI_95_Lower_80-85", "logIRR_CI_95_Lower_85-90", 
               "logIRR_CI_95_Lower_90-95")

#Making a table of point estimates
point_estimates <- 
  data.frame("Ages" = c(80, 85, 90), 
             "ACT" = ACT_inc_rates$`All_Dem_IRR_F:M`[
               which(ACT_inc_rates$Low_Age == 80):
                 which(ACT_inc_rates$Low_Age == 90)],
             "A" = exp(mean_results_A[IRR_ages]), 
             "B1" = exp(mean_results_B1[IRR_ages]), 
             "B2" = exp(mean_results_B2[IRR_ages]), 
             "C1" = exp(mean_results_C1[IRR_ages]), 
             "C2" = exp(mean_results_C2[IRR_ages])) %>% 
  pivot_longer(cols = c("ACT", "A", "B1", "B2", "C1", "C2"), 
               names_to = "Scenario", values_to = "IRR")
                            
#Making a table of the lower bounds for the 95%CI
LB <- data.frame("Ages" = c(80, 85, 90), 
                 "ACT" = ACT_inc_rates$`All_Dem_95CI_lower`[
                   which(ACT_inc_rates$Low_Age == 80):
                     which(ACT_inc_rates$Low_Age == 90)], 
                 "A" = exp(mean_results_A[IRR_lower]), 
                 "B1" = exp(mean_results_B1[IRR_lower]), 
                 "B2" = exp(mean_results_B2[IRR_lower]), 
                 "C1" = exp(mean_results_C1[IRR_lower]), 
                 "C2" = exp(mean_results_C2[IRR_lower])) %>% 
  pivot_longer(cols = c("ACT", "A", "B1", "B2", "C1", "C2"), 
               names_to = "Scenario", values_to = "LB")

#Making a table of the upper bounds for the 95%CI
UB <- data.frame("Ages" = c(80, 85, 90), 
                 "ACT" = ACT_inc_rates$`All_Dem_95CI_upper`[
                   which(ACT_inc_rates$Low_Age == 80):
                     which(ACT_inc_rates$Low_Age == 90)], 
                 "A" = exp(mean_results_A[IRR_upper]), 
                 "B1" = exp(mean_results_B1[IRR_upper]), 
                 "B2" = exp(mean_results_B2[IRR_upper]), 
                 "C1" = exp(mean_results_C1[IRR_upper]), 
                 "C2" = exp(mean_results_C2[IRR_upper])) %>% 
  pivot_longer(cols = c("ACT", "A", "B1", "B2", "C1", "C2"), 
               names_to = "Scenario", values_to = "UB")

#Data frame with error bar data for all simulation scenarios
figure2_data_all_error <- 
  left_join(point_estimates, LB, 
            by = c("Ages", "Scenario")) %>%
  left_join(., UB, by = c("Ages", "Scenario")) %>% 
  mutate_at("Scenario", as.factor)
figure2_data_all_error$Scenario <- 
  relevel(figure2_data_all_error$Scenario, "ACT")

#Data frame with error bar data for only the ACT study
figure2_data_ACT_error <- 
  left_join(point_estimates, LB, 
            by = c("Ages", "Scenario")) %>%
  left_join(., UB, by = c("Ages", "Scenario")) %>% 
  mutate_at("Scenario", as.factor) 
figure2_data_ACT_error$UB[figure2_data_ACT_error$Scenario != "ACT"] <- NA
figure2_data_ACT_error$Scenario <- 
  relevel(figure2_data_ACT_error$Scenario, "ACT")

#Plot with error bars for all scenarios
figure2_option1 <- ggplot(figure2_data_all_error, 
                          aes(x = Ages, y = IRR, colour = Scenario)) + 
  geom_errorbar(aes(ymin = LB, ymax = UB), width = 1, 
                position = position_dodge(width = 0.5)) +
  geom_line(position = position_dodge(width = 0.5)) + 
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_hline(yintercept = 1, linetype="dashed", color = "black") +
  theme_minimal() + 
  scale_x_continuous(name = "Age bands", breaks = c(80, 85, 90), 
                     labels = c("[80-85)", "[85-90)","[90-95)")) + 
  ylab(TeX("$\\widehat{\\bar{IRR}}_{women:men}$"))

#Plot with error bars only for the ACT study
figure2_option2 <- ggplot(figure2_data_ACT_error, 
                          aes(x = Ages, y = IRR, colour = Scenario)) + 
  geom_errorbar(aes(ymin = LB, ymax = UB), width = 0.5) +
  geom_line() + geom_point() + 
  geom_hline(yintercept = 1, linetype="dashed", color = "black") + 
  theme_minimal() + 
  scale_x_continuous(name = "Age bands", breaks = c(80, 85, 90), 
  labels = c("[80-85)", "[85-90)","[90-95)")) + 
  ylab(TeX("$\\widehat{\\bar{IRR}}_{women:men}$")) 

#Saving figures
ggsave(here("Manuscript", "figure2_option1.jpeg"), plot = figure2_option1,
       device = "jpeg", dpi = 300)

ggsave(here("Manuscript", "figure2_option2.jpeg"), plot = figure2_option2,
       device = "jpeg", dpi = 300)

#---- Figure 3 ----
#Create one sample of size 100,000 for each simulation scenario
source(here("RScripts","Scenario_A_pars.R"))
sample_A <- data_gen(num_obs = 100000)

source(here("RScripts","Scenario_B1_pars.R"))
sample_B1 <- data_gen(num_obs = 100000)

source(here("RScripts","Scenario_B2_pars.R"))
sample_B2 <- data_gen(num_obs = 100000)

source(here("RScripts","Scenario_C1_pars.R"))
sample_C1 <- data_gen(num_obs = 100000)

source(here("RScripts","Scenario_C2_pars.R"))
sample_C2 <- data_gen(num_obs = 100000)

#Format data for plotting
plot_data_A <- format_plot_data(sample_A, "A") 
plot_data_B1 <- format_plot_data(sample_B1, "B1") 
plot_data_B2 <- format_plot_data(sample_B2, "B2") 
plot_data_C1 <- format_plot_data(sample_C1, "C1") 
plot_data_C2 <- format_plot_data(sample_C2, "C2") 

figure3_data <- do.call("rbind", 
                        list(plot_data_A, plot_data_B1, plot_data_B2, 
                             plot_data_C1, plot_data_C2))

#Plot
figure3 <- ggplot(data = figure3_data, 
                  aes(x = `Age Band` , y = U, 
                      fill = `Sex/Gender`, color = `Sex/Gender`)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 0, linetype="dashed", color = "black") + 
  theme_minimal() + 
  theme(panel.spacing = unit(1, "lines")) +
  facet_grid(. ~Scenario) + 
  coord_flip()

#Saving figure
ggsave(here("Manuscript", "figure3.jpeg"), plot = figure3,
       device = "jpeg", dpi = 300)
  






