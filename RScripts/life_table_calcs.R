#---- Package loading, options, seed ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "magrittr", "here")

#---- Source files ----
source(here("RScripts", "US_life_table_2017.R"))

#---- Selected life table subsets ----
#Add computations for cumulative survival conditional on survival to age 50
all_life_US <- white_life %>% 
  filter(`Birth cohort` == "1919-1921" & Age %in% seq(50, 95, by = 5)) 
all_life_US %<>% 
  mutate("cum_surv_cond50" = Survivors/all_life_US$Survivors[1])

female_life_US <- white_female_life %>% 
  filter(`Birth cohort` == "1919-1921" & Age %in% seq(50, 95, by = 5)) 
female_life_US %<>%
  mutate("cum_surv_cond50" = Survivors/female_life_US$Survivors[1])

male_life_US <- white_male_life %>% 
  filter(`Birth cohort` == "1919-1921" & Age %in% seq(50, 95, by = 5)) 
male_life_US %<>%
  mutate("cum_surv_cond50" = Survivors/male_life_US$Survivors[1])

#---- Hazard ratios ----
male_haz_US <- white_male_life %>% 
  filter(`Birth cohort` == "1919-1921" & Age %in% seq(50, 95, by = 5)) %>% 
  dplyr::select("Birth cohort", "Haz")

female_haz_US <- white_female_life %>% 
  filter(`Birth cohort` == "1919-1921" & Age %in% seq(50, 95, by = 5)) %>% 
  dplyr::select("Birth cohort", "Haz")

Hratio_US <- female_haz_US$Haz/male_haz_US$Haz %>%
  as.data.frame() %>% set_colnames("US") 

