#---- Package loading, options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "here", "scales")

#---- Source files ----
source(here("RScripts", "dementia_incidence_ACT.R"))

#---- ACT AD by sex/gender ----
ACT_plot_data <- 
  tibble("Age" = seq(69, 94, by = 5), 
         "Men" = ACT_inc_rates$Male_AD_1000PY, 
         "Women" = ACT_inc_rates$Female_AD_1000PY) %>%
  gather(key = "Sex/Gender", value = "value", c("Men", "Women")) %>% 
  mutate_at("Sex/Gender", as.factor)
EURODEM_plot_data$`Sex/Gender` <- 
  fct_relevel(EURODEM_plot_data$`Sex/Gender`, "Men", after = 1)

ggplot(aes(Age, value), data = EURODEM_plot_data) + 
  geom_point(aes(colour = `Sex/Gender`)) + theme_minimal() + 
  geom_line(aes(color = `Sex/Gender`)) + 
  #ggtitle("Dementia Incidence Rates by Sex from EURODEM Pooled Analysis ") + 
  ylab("Dementia Incidence Rate per 1000 Person Years")

#---- IRR plot ----
IRR_data <- tibble("Age" = seq(69, 94, by = 5), 
                   "ACT" = ACT_inc_rates$`All_Dem_IRR_F:M`,
                   "A" = c(1.002, 1.001, 1.003, 1.002, 1.003, 0.999), 
                   "B1" = c(1, 1, 1, 1, 1, 1), 
                   "B2" = c(0.999, 1.007, 1.017, 1.012, 1.015, 1.002), 
                   "C1" = c(1, 1.08, 1.16, 1.16, 1.15, 1.17), 
                   "C2" = c(1.01, 1.12, 1.23, 1.21, 1.20, 1.22)) %>%
  gather(key = "Scenario", value = "value", 
         c("ACT", "A", "B1", "B2", "C1", "C2")) %>% 
  mutate_at("Scenario", as.factor)
IRR_data$`Scenario` <- relevel(IRR_data$`Scenario`, "ACT")

ggplot(aes(Age, value), data = IRR_data) + 
  geom_point(aes(colour = `Scenario`)) + theme_minimal() + 
  theme(text = element_text(size = 24)) +
  geom_line(aes(color = `Scenario`)) + 
  scale_color_manual(values = c("ACT" = "Black", "A" = "#F8766D", 
                                "B1" = "#C77CFF", "B2" = "Blue", 
                                "C1" = "#00BFC4", "C2" = "#7CAE00")) + 
  geom_hline(yintercept = 1, lty = 2, color = "black", size = 0.75)
  

#---- ACT All Dem incidence matching plot (men) ----
dem_match_plot_data <- 
  tibble("Age" = seq(69, 94, by = 5), 
         "ACT" = ACT_inc_rates$Male_All_Dementia_1000PY, 
         "A" = c(7.03, 8.87, 21.19, 45.48, 71.50, 96.96), 
         "B1" = c(7.04, 8.84, 21.15, 45.57, 71.76, 97.38), 
         "B2" = c(7.06, 8.92, 21.06, 45.65, 72.48, 99.38), 
         "C1" = c(7.07, 9.11, 21.70, 46.29, 73.20, 99.56), 
         "C2" = c(7.05, 8.84, 20.68, 44.75, 71.24, 96.85)) %>%
  gather(key = "Scenario", value = "value", 
         c("ACT", "A", "B1", "B2", "C1", "C2")) %>% 
  mutate_at("Scenario", as.factor)
dem_match_plot_data$`Scenario` <- 
  relevel(dem_match_plot_data$`Scenario`, "ACT")

ggplot(aes(Age, value), data = dem_match_plot_data) + 
  geom_point(aes(colour = `Scenario`)) + theme_minimal() + 
  #ggtitle("Dementia Incidence Rates by Sex from EURODEM Pooled Analysis ") + 
  #ylab("AD Incidence Rate per 1000 Person Years") +  
  theme(text = element_text(size = 24)) +
  geom_line(aes(color = `Scenario`)) + 
  scale_color_manual(values = c("ACT" = "Black", "A" = "#F8766D", 
                                "B1" = "#C77CFF", "B2" = "Blue", 
                                "C1" = "#00BFC4", "C2" = "#7CAE00"))
