#---- Package Loading, Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("here", "tidyverse", "latex2exp", "magrittr", "harrypotter", "grid")
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

# #Plot with error bars for all scenarios
# figure2_option1 <- ggplot(figure2_data_all_error,
#                           aes(x = Ages, y = IRR, colour = Scenario)) +
#   geom_errorbar(aes(ymin = LB, ymax = UB), width = 1,
#                 position = position_dodge(width = 0.5)) +
#   geom_line(position = position_dodge(width = 0.5)) +
#   geom_point(position = position_dodge(width = 0.5)) +
#   geom_hline(yintercept = 1, linetype = "solid", color = "black") +
#   theme_minimal() +
#   scale_x_continuous(name = "Age bands", breaks = c(80, 85, 90),
#                      labels = c("[80-85)", "[85-90)","[90-95)")) +
#   ylab(TeX("$\\widehat{\\bar{IRR}}_{women:men}$"))

# #Plot with error bars only for the ACT study
# figure2_option2 <- ggplot(figure2_data_ACT_error, 
#                           aes(x = Ages, y = IRR, colour = Scenario)) + 
#   geom_errorbar(aes(ymin = LB, ymax = UB), width = 0.25, size = 1) +
#   geom_line(aes(color = Scenario), size = 1) + 
#   geom_point(aes(color = Scenario, shape = Scenario), size = 3) + 
#   geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
#   scale_x_continuous(name = "Age bands", breaks = c(80, 85, 90), 
#   labels = c("[80-85)", "[85-90)","[90-95)")) + 
#   scale_colour_manual(name = "Scenario",
#     values = c("darkgrey", "black", "#A2D5C6", "#A2D5C6","#039FBE", "#039FBE"), 
#     labels = levels(figure2_data_ACT_error$Scenario)) + 
#   scale_shape_manual(values = c(19, 8, 19, 17, 19, 17)) + 
#   # scale_linetype_manual(name = "Scenario", 
#   #                       values = c(rep("solid", 3), rep("twodash", 3)), 
#   #                       labels = levels(figure2_data_ACT_error$Scenario)) + 
#   ylab(TeX("$\\widehat{\\bar{IRR}}_{women:men}$")) +
#   theme_minimal() + 
#   theme(text = element_text(size = 14))

#edit data for this plot
#Age data
figure2_data_all_error$Ages <- as.character(figure2_data_all_error$Ages)
figure2_data_all_error[which(figure2_data_all_error$Ages == "80"), "Ages"] <- 
  "[80, 85)"
figure2_data_all_error[which(figure2_data_all_error$Ages == "85"), "Ages"] <- 
  "[85, 90)"
figure2_data_all_error[which(figure2_data_all_error$Ages == "90"), "Ages"] <- 
  "[90, 95)"
#Scenario data
figure2_data_all_error$Scenario <- as.character(figure2_data_all_error$Scenario)
figure2_data_all_error[which(figure2_data_all_error$Scenario == "A"), 
                       "Scenario"] <- "No Selection"
figure2_data_all_error[which(figure2_data_all_error$Scenario == "B1"), 
                       "Scenario"] <- "HOM1"
figure2_data_all_error[which(figure2_data_all_error$Scenario == "B2"), 
                       "Scenario"] <- "HOM2"
figure2_data_all_error[which(figure2_data_all_error$Scenario == "C1"), 
                       "Scenario"] <- "HET1"
figure2_data_all_error[which(figure2_data_all_error$Scenario == "C2"), 
                       "Scenario"] <- "HET2"
figure2_data_all_error$Scenario <- as.factor(figure2_data_all_error$Scenario)
figure2_data_all_error$Scenario <- 
  factor(figure2_data_all_error$Scenario, levels = 
          c("ACT", "No Selection", "HOM1", "HOM2", "HET1", "HET2"))


#plot
figure2_option3 <- ggplot(figure2_data_all_error, 
                          aes(x = Ages, y = IRR, fill = Scenario)) +
  stat_summary(fun.y = identity, geom = "bar", 
               position = position_dodge(width = 0.9), size = 3) +
  geom_errorbar(aes(ymin = LB, ymax = UB),
                width = .2, position = position_dodge(0.9)) + 
  scale_fill_hp_d(option = "LunaLovegood", begin = 0, end = 1) +
  theme_minimal() + ylab(TeX("$\\widehat{\\bar{IRR}}_{women:men}$")) + 
  theme(text = element_text(size = 14)) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")


#Saving figures
# ggsave(here("Manuscript", "figure2_option1.jpeg"), plot = figure2_option1,
#        device = "jpeg", dpi = 300, width = 6.5, height = 5.5, units = "in")

# ggsave(here("Manuscript", "figure2_option2.pdf"), plot = figure2_option2,
#        device = "pdf", dpi = 300)

ggsave(here("Manuscript", "figure2_option3.pdf"), plot = figure2_option3,
       device = "pdf", dpi = 300)

#---- Figure 3 ----
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

#Format data for plotting
plot_data_A <- format_plot_data(sample_A, "No Selection") 
plot_data_B1 <- format_plot_data(sample_B1, "HOM1") 
plot_data_B2 <- format_plot_data(sample_B2, "HOM2") 
plot_data_C1 <- format_plot_data(sample_C1, "HET1") 
plot_data_C2 <- format_plot_data(sample_C2, "HET2") 

figure3_data <- do.call("rbind", 
                        list(plot_data_A, plot_data_B1, plot_data_B2, 
                             plot_data_C1, plot_data_C2))

figure3_data$Scenario <- 
  factor(figure3_data$Scenario, levels = 
           c("ACT", "No Selection", "HOM1", "HOM2", "HET1", "HET2"))


#Manuscript calculation (Results section)
mean_U_last <- figure3_data %>% 
  filter(`Age Band` == "[90-95)") %>% 
  group_by(Scenario, `Sex/Gender`) %>% 
  summarise_at("U", mean)

#Plot
figure3 <- ggplot(data = figure3_data) + 
  geom_boxplot(aes(x = `Age Band` , y = U, 
                   fill = `Sex/Gender`, color = `Sex/Gender`)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  theme_minimal() + 
  theme(panel.spacing = unit(1, "lines")) +
  facet_grid(. ~Scenario) + 
  coord_flip() + 
  scale_fill_hp_d(option = "LunaLovegood", begin = 0.3, end = 1) +
  scale_color_hp_d(option = "LunaLovegood", begin = 0.3, end = 1) +
  theme(text = element_text(size = 14), 
        panel.spacing = unit(1.5, "lines")) + 
  guides(fill = guide_legend(reverse = TRUE)) + 
  guides(color = guide_legend(reverse = TRUE))

#Saving figure
ggsave(here("Manuscript", "figure3.pdf"), plot = figure3,
       device = "pdf", dpi = 300)

#---- eFigure 1 ----
#(a)
#Survival plot data
female_cp_survival <- 
  cbind(mean_results_A[variable_names$cp_alive_females_varnames[1:num_tests]], 
        mean_results_B1[variable_names$cp_alive_females_varnames[1:num_tests]], 
        mean_results_B2[variable_names$cp_alive_females_varnames[1:num_tests]], 
        mean_results_C1[variable_names$cp_alive_females_varnames[1:num_tests]], 
        mean_results_C2[variable_names$cp_alive_females_varnames[1:num_tests]],
        female_life_US$CP[-1]) %>% as.data.frame() %>%
  set_colnames(c("No Selection", "HOM1", "HOM2", "HET1", "HET2", "Lifetable")
               ) %>% mutate("Age" = seq(55, 95, by = 5)) %>% 
  pivot_longer(cols = -Age, names_to = "Scenario", values_to = "prob") %>% 
  mutate("Sex/Gender" = "Women")
  
male_cp_survival <- 
  cbind(mean_results_A[variable_names$cp_alive_males_varnames[1:num_tests]], 
        mean_results_B1[variable_names$cp_alive_males_varnames[1:num_tests]], 
        mean_results_B2[variable_names$cp_alive_males_varnames[1:num_tests]], 
        mean_results_C1[variable_names$cp_alive_males_varnames[1:num_tests]], 
        mean_results_C2[variable_names$cp_alive_males_varnames[1:num_tests]],
        male_life_US$CP[-1]) %>% as.data.frame() %>%
  set_colnames(c("No Selection", "HOM1", "HOM2", "HET1", "HET2", "Lifetable")
               ) %>% mutate("Age" = seq(55, 95, by = 5)) %>% 
  pivot_longer(cols = -Age, names_to = "Scenario", values_to = "prob") %>% 
  mutate("Sex/Gender" = "Men")

cp_surv_plot_data <- rbind(female_cp_survival, male_cp_survival) %>% 
  mutate_at("Scenario", as.factor)
cp_surv_plot_data$Scenario <- 
  factor(cp_surv_plot_data$Scenario, levels = 
           c("Lifetable", "No Selection", "HOM1", "HOM2", "HET1", "HET2"))


figure_e1a <- ggplot(cp_surv_plot_data, 
                     aes(Age, prob, group = Scenario, color = Scenario)) + 
  geom_line(size = 1.25) + 
  scale_x_continuous(breaks = seq(50, 95, 5)) + 
  scale_color_hp_d(option = "LunaLovegood", begin = 0, end = 1) +
  labs(title = "", 
       y = "Conditional Survival Probability from age 50", 
       x = "Age") + 
  facet_grid(. ~`Sex/Gender`) + 
  theme_minimal() + 
  theme(text = element_text(size = 14))

ggsave(here("Manuscript", "figure_e1a.pdf"), plot = figure_e1a,
       device = "pdf", dpi = 300)

#(b)
HR_plot_data <- 
  cbind(Hratio_US[-1, ], 
        exp(mean_results_A[na.omit(variable_names$mortality_logHR_varnames)]), 
        exp(mean_results_B1[na.omit(variable_names$mortality_logHR_varnames)]), 
        exp(mean_results_B2[na.omit(variable_names$mortality_logHR_varnames)]), 
        exp(mean_results_C1[na.omit(variable_names$mortality_logHR_varnames)]), 
        exp(mean_results_C2[na.omit(
          variable_names$mortality_logHR_varnames)])) %>% 
  set_colnames(c("Lifetable", "No Selection", "HOM1", "HOM2", "HET1", "HET2")
               ) %>% as.data.frame() %>%
  mutate("Age" = seq(50, 90, by = 5)) %>%
  pivot_longer(cols = -Age, 
               names_to = "Scenario", 
               values_to = "HR") %>% 
  mutate_at("Scenario", as.factor)
HR_plot_data$Scenario <- 
  factor(HR_plot_data$Scenario, levels = 
           c("Lifetable", "No Selection", "HOM1", "HOM2", "HET1", "HET2"))

figure_e1b <- ggplot(HR_plot_data, aes(Age, HR)) + 
  geom_point(aes(color = Scenario, group = Scenario), 
             size = 3, alpha = rep(c(1, rep(0.3, 5)), 9)) + 
  geom_line(aes(color = Scenario, group = Scenario), size = 1.25, 
            alpha = rep(c(1, rep(0.3, 5)), each = 9)) + 
  scale_x_continuous(name = "Age bands", breaks = seq(50, 90, 5), 
                     labels = c("[50-55)", "[55-60)","[60-65)", "[65-70)", 
                                "[70-75)","[75-80)", "[80-85)", "[85-90)",
                                "[90-95)")) + 
  #scale_shape_manual(values = c(19, 8, 19, 17, 19, 17)) + 
  scale_color_hp_d(option = "LunaLovegood", begin = 0, end = 1) + 
  labs(y = "Mortality Hazard Ratio (Women:Men)", x = "Age", 
       color = "") + theme_minimal() + 
  theme(text = element_text(size = 14)) 

ggsave(here("Manuscript", "figure_e1b.pdf"), plot = figure_e1b,
       device = "pdf", dpi = 300)


#---- eFigure 2 ----
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

#Create plot data
Ci_data <- 
  rbind(sample_A %>% 
          dplyr::select("id", "female", variable_names$Ci_varnames, 
                        "survtime", "last_Ci") %>% 
          mutate("Gender" = case_when(female == 0 ~ "Men", 
                                      TRUE ~ "Women"), 
                 "Scenario" = "No Selection") %>% 
          dplyr::select(-one_of("female")), 
        sample_B1 %>% 
          dplyr::select("id", "female", variable_names$Ci_varnames, 
                        "survtime", "last_Ci") %>% 
          mutate("Gender" = case_when(female == 0 ~ "Men", 
                                      TRUE ~ "Women"), 
                 "Scenario" = "HOM1") %>% 
          dplyr::select(-one_of("female")), 
        sample_B2 %>% 
          dplyr::select("id", "female", variable_names$Ci_varnames, 
                        "survtime", "last_Ci") %>% 
          mutate("Gender" = case_when(female == 0 ~ "Men", 
                                      TRUE ~ "Women"), 
                 "Scenario" = "HOM2") %>% 
          dplyr::select(-one_of("female")), 
        sample_C1 %>% 
          dplyr::select("id", "female", variable_names$Ci_varnames, 
                        "survtime", "last_Ci") %>% 
          mutate("Gender" = case_when(female == 0 ~ "Men", 
                                      TRUE ~ "Women"), 
                 "Scenario" = "HET1") %>% 
          dplyr::select(-one_of("female")), 
        sample_C2 %>% 
          dplyr::select("id", "female", variable_names$Ci_varnames, 
                        "survtime", "last_Ci") %>% 
          mutate("Gender" = case_when(female == 0 ~ "Men", 
                                      TRUE ~ "Women"), 
                 "Scenario" = "HET2") %>% 
          dplyr::select(-one_of("female")))

#Getting mean Ci by sex
mean_Ci_data <- Ci_data %>% group_by(Gender, Scenario) %>% 
  summarise_at(variable_names$Ci_varnames, ~mean(., na.rm = TRUE)) %>% 
  mutate("Cij" = NA, 
         "Age" = NA) %>%
  dplyr::select("Gender", variable_names$Ci_varnames, "Scenario", "Age") %>% 
  set_colnames(c("id", seq(50, 95, by = 5), "Scenario", "Age"))

#Sample data 
Ci_plot_data <- Ci_data %>% group_by(Scenario) %>% sample_n(100) %>% 
  mutate("Age" = survtime + 50) %>% 
  dplyr::select(-c(survtime, Gender)) %>% 
  set_colnames(c("id", seq(50, 95, by = 5), 
                 "Cij", "Scenario", "Age")) %>% 
  mutate_at("id", as.character) %>%
  rbind(mean_Ci_data) 


#Plot data
samp_Ci <- Ci_plot_data[, c("id", seq(50, 95, by = 5), "Scenario")] %>%
  gather(as.character(seq(50, 95, by = 5)), key = "Age", value = "Cij") %>% 
  mutate_at("Age", as.numeric) %>% 
  rbind(., Ci_plot_data[, c("id", "Scenario", "Age", "Cij")]) %>% 
  set_colnames(c("variable", "Scenario", "Age", "value"))

samp_Ci$Scenario <- 
  factor(samp_Ci$Scenario, levels = 
           c("ACT", "No Selection", "HOM1", "HOM2", "HET1", "HET2"))



#Creating a plot with random sample in the background
figure_e2 <- ggplot(samp_Ci, aes(Age, value)) + 
  geom_line(data = 
              subset(samp_Ci, variable != "Women" & variable != "Men"), 
            aes(group = variable), color = "gray") +
  geom_line(data = subset(samp_Ci, variable == "Women"), 
            aes(color = variable), size = 1.25) + 
  geom_line(data = subset(samp_Ci, variable == "Men"), 
            aes(color = variable), size = 1.25, alpha = 0.6) + 
  geom_hline(yintercept = dem_cut, size = 1.25) + 
  labs(y = "Cognitive Function", 
       x = "Age", 
       color = "Sex/Gender") + 
  scale_x_continuous(breaks = seq(50, 95, 5)) + 
  facet_grid(. ~Scenario) +
  theme_minimal() +  
  theme(text = element_text(size = 14)) + 
  coord_cartesian(ylim = c(-10, 10)) + 
  scale_color_hp_d(option = "LunaLovegood", begin = 0.3, end = 1) + 
  #guides(color = guide_legend(reverse = TRUE)) +
  geom_hline(yintercept = dem_cut) 

ggsave(here("Manuscript", "figure_e2.jpeg"), plot = figure_e2,
       device = "jpeg", dpi = 300)

#---- eFigure 3 ----
male_inc_data <- 
  data.frame("Ages" = seq(50, 90, by = 5), 
             "ACT" = c(rep(0, 3),ACT_inc_rates$Male_All_Dementia_1000PY),
             "No Selection" = 1000*
               mean_results_A[na.omit(
                 variable_names$inc_cases_males_varnames)]/
               mean_results_A[na.omit(variable_names$PY_males_varnames)], 
             "HOM1" = 1000*
               mean_results_B1[na.omit(
                 variable_names$inc_cases_males_varnames)]/
               mean_results_B1[na.omit(variable_names$PY_males_varnames)], 
             "HOM2" = 1000*
               mean_results_B2[na.omit(
                 variable_names$inc_cases_males_varnames)]/
               mean_results_B2[na.omit(variable_names$PY_males_varnames)], 
             "HET1" = 1000*
               mean_results_C1[na.omit(
                 variable_names$inc_cases_males_varnames)]/
               mean_results_C1[na.omit(variable_names$PY_males_varnames)], 
             "HET2" = 1000*
               mean_results_C2[na.omit(
                 variable_names$inc_cases_males_varnames)]/
               mean_results_C2[na.omit(variable_names$PY_males_varnames)]) %>% 
  pivot_longer(cols = c("ACT", "No.Selection", "HOM1", "HOM2", "HET1", "HET2"), 
               names_to = "Scenario", values_to = "rates") %>% 
  filter(Ages >= 65) %>% mutate_at("Scenario", as.factor)  

male_inc_data$Scenario <- 
  factor(male_inc_data$Scenario, levels = 
           c("ACT", "No.Selection", "HOM1", "HOM2", "HET1", "HET2"))


#Make the plot
figure_e3 <- ggplot(male_inc_data, 
                    aes(x = Ages, y = rates, colour = Scenario)) + 
  geom_line() + geom_point(aes(shape = Scenario)) + theme_minimal() + 
  scale_x_continuous(name = "Age bands", breaks = seq(65, 90, by = 5), 
                     labels = c("[65-70)", "[70-75)","[75-80)", 
                                "[80-85)", "[85-90)","[90-95)")) + 
  ylab(TeX("Dementia incidence rates per 1000 PY")) + 
  ggtitle("Dementia incidence rates for men") + 
  scale_colour_manual(name = "Scenario",
                      values = c("darkgrey", "black", "#A2D5C6", "#A2D5C6",
                                 "#039FBE", "#039FBE"), 
                      labels = levels(male_inc_data$Scenario)) + 
  scale_shape_manual(values = c(19, 8, 19, 17, 19, 17)) + 
  theme(text = element_text(size = 14)) 

ggsave(here("Manuscript", "figure_e3.jpeg"), plot = figure_e3,
       device = "jpeg", dpi = 300)

#---- eTable 4 ----
exp(mean_results_A[variable_names$logIRR_varnames])
exp(mean_results_A[variable_names$logIRR_95CI_Lower_varnames])
exp(mean_results_A[variable_names$logIRR_95CI_Upper_varnames])

exp(mean_results_B1[variable_names$logIRR_varnames])
exp(mean_results_B1[variable_names$logIRR_95CI_Lower_varnames])
exp(mean_results_B1[variable_names$logIRR_95CI_Upper_varnames])

exp(mean_results_B2[variable_names$logIRR_varnames])
exp(mean_results_B2[variable_names$logIRR_95CI_Lower_varnames])
exp(mean_results_B2[variable_names$logIRR_95CI_Upper_varnames])

exp(mean_results_C1[variable_names$logIRR_varnames])
exp(mean_results_C1[variable_names$logIRR_95CI_Lower_varnames])
exp(mean_results_C1[variable_names$logIRR_95CI_Upper_varnames])

exp(mean_results_C2[variable_names$logIRR_varnames])
exp(mean_results_C2[variable_names$logIRR_95CI_Lower_varnames])
exp(mean_results_C2[variable_names$logIRR_95CI_Upper_varnames])






