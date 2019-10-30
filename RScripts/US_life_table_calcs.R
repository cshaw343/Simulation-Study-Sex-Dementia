#---- Package loading, options, seed ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "magrittr", "here")

#---- Read in data ----
US_total <- read_csv(here("Data", "US_cohort_table_white_MF_calcs.csv"))
US_F <- read_csv(here("Data", "US_cohort_table_white_F_calcs.csv"))
US_M <- read_csv(here("Data", "US_cohort_table_white_M_calcs.csv"))

#---- Selected life table subsets ----
#Add computations for cumulative survival conditional on survival to age 50
all_life_US <- US_total %>% 
  filter(`Birth cohort` == "1919-1921" & Age %in% seq(50, 95, by = 5)) 
all_life_US %<>% 
  mutate("cum_surv_cond50" = Survivors/all_life_US$Survivors[1])

female_life_US <- US_F %>% 
  filter(`Birth cohort` == "1919-1921" & Age %in% seq(50, 95, by = 5)) 
female_life_US %<>%
  mutate("cum_surv_cond50" = Survivors/female_life_US$Survivors[1])

male_life_US <- US_M %>% 
  filter(`Birth cohort` == "1919-1921" & Age %in% seq(50, 95, by = 5)) 
male_life_US %<>%
  mutate("cum_surv_cond50" = Survivors/male_life_US$Survivors[1])

#---- Hazard ratios ----
#US
male_haz_US <- US_M %>% 
  filter(`Birth cohort` == "1919-1921" & Age %in% seq(50, 95, by = 5)) %>% 
  dplyr::select("Birth cohort", "Haz")

female_haz_US <- US_F %>% 
  filter(`Birth cohort` == "1919-1921" & Age %in% seq(50, 95, by = 5)) %>% 
  dplyr::select("Birth cohort", "Haz")

Hratio_US <- female_haz_US$Haz/male_haz_US$Haz %>%
  as.data.frame() %>% set_colnames("US") 

