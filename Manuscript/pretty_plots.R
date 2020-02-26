#---- Package Loading, Options, and Seed ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("here", "tidyverse", "latex2exp", "magrittr")

#---- Source scripts ----
source(here("RScripts", "scenario_A_pars.R")) #Need for parameters
source(here("RScripts", "var_names.R"))
source(here("RScripts", "data_gen.R"))
source(here("RScripts", "dem_inc_ACT.R"))
source(here("RScripts", "format_plot_data.R"))
source(here("RScripts", "results_65_plus.R"))
source(here("RScripts", "life_table_calcs.R"))

#---- Load results data ----
results_A <- read_csv(here("Results", "Scenario_A_no_bias", 
                           "sim_A_1000_20191125.csv")) %>% results_65_plus()

results_B1 <- read_csv(here("Results", "Scenario_B", 
                            "sim_B1_1000_20191128.csv")) %>% results_65_plus()

results_B2 <- read_csv(here("Results", "Scenario_B", 
                            "sim_B2_1000_20191126.csv")) %>% results_65_plus()

results_C1 <- read_csv(here("Results", "Scenario_C", 
                            "sim_C1_1000_20191129.csv")) %>% results_65_plus()

results_C2 <- read_csv(here("Results", "Scenario_C", 
                            "sim_C2_1000_20191126.csv")) %>% results_65_plus()

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
sample_A <- data_gen(num_obs = 100000)

source(here("RScripts","scenario_B1_pars.R"))
sample_B1 <- data_gen(num_obs = 100000)

source(here("RScripts","scenario_B2_pars.R"))
sample_B2 <- data_gen(num_obs = 100000)

source(here("RScripts","scenario_C1_pars.R"))
sample_C1 <- data_gen(num_obs = 100000)

source(here("RScripts","scenario_C2_pars.R"))
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

#Manuscript calculation (Results section)
mean_U_last <- figure3_data %>% 
  filter(`Age Band` == "[90-95)") %>% 
  group_by(Scenario, `Sex/Gender`) %>% 
  summarise_at("U", mean)

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

#---- Figure e1 ----
#(a)
#Survival plot data
female_cp_survival <- 
  c(mean_results_A[variable_names$cp_alive_females_varnames[1:num_tests]])
male_cp_survival <- 
  c(mean_results_A[variable_names$cp_alive_males_varnames[1:num_tests]])

cp_surv_plot_data <- tibble("age" = seq(55, 95, by = 5), 
                            "Published_Women_survival" = female_life_US$CP[-1], 
                            "Published_Men_survival" = male_life_US$CP[-1], 
                            "Simulated_Women_survival" = female_cp_survival, 
                            "Simulated_Men_survival" = male_cp_survival) %>% 
  pivot_longer(cols = -age, 
               names_to = c("Data Type", "Gender"),
               names_pattern = "(.*)_(.*)?_survival",
               values_to = "prob")

figure_e1a <- ggplot(cp_surv_plot_data, 
                     aes(age, prob, group = `Data Type`, color = `Data Type`)) + 
  geom_line(size = 1.25) + 
  scale_x_continuous(breaks = seq(50, 95, 5)) + 
  labs(title = "", 
       y = "Conditional Survival Probability from age 50", 
       x = "Age") + 
  facet_grid(. ~Gender) + 
  theme_minimal() 

ggsave(here("Manuscript", "figure_e1a.jpeg"), plot = figure_e1a,
       device = "jpeg", dpi = 300)

#(b)
HR_plot_data <- 
  cbind(Hratio_US[-1, ], 
        exp(mean_results_A[na.omit(variable_names$mortality_logHR_varnames)]), 
        exp(mean_results_B1[na.omit(variable_names$mortality_logHR_varnames)]), 
        exp(mean_results_B2[na.omit(variable_names$mortality_logHR_varnames)]), 
        exp(mean_results_C1[na.omit(variable_names$mortality_logHR_varnames)]), 
        exp(mean_results_C2[na.omit(
          variable_names$mortality_logHR_varnames)])) %>% 
  set_colnames(c("Published", "A", "B1", "B2", "C1", "C2")) %>% 
  as.data.frame() %>%
  mutate("Age" = seq(50, 90, by = 5)) %>%
  pivot_longer(cols = -Age, 
               names_to = "Scenario", 
               values_to = "HR") %>% 
  mutate_at("Scenario", as.factor)
HR_plot_data$Scenario <- fct_relevel(HR_plot_data$Scenario,"Published", 
                                     after = 0)

figure_e1b <- ggplot(HR_plot_data, aes(Age, HR)) + 
  geom_point(aes(color = Scenario, group = Scenario), size = 1.75) + 
  geom_line(aes(color = Scenario, group = Scenario), size = 1.25, 
            alpha = 0.6) + 
  scale_x_continuous(name = "Age bands", breaks = seq(50, 90, 5), 
                     labels = c("[50-55)", "[55-60)","[60-65)", "[65-70)", 
                                "[70-75)","[75-80)", "[80-85)", "[85-90)",
                                "[90-95)")) + 
  #ylim(0, 1) +
  #scale_y_continuous(breaks = seq(0, 1, 0.05)) +
  labs(y = "Mortality Hazard Ratio (Women:Men)", x = "Age", 
       color = "") + theme_minimal() + 
  #theme(text = element_text(size = 40)) + 
  ggtitle("")

ggsave(here("Manuscript", "figure_e1b.jpeg"), plot = figure_e1b,
       device = "jpeg", dpi = 300)
        
        
  
  









