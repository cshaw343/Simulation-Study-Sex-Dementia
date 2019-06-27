#---- Package loading, options, seed ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "magrittr", "here")

#---- Read in data ----
Netherlands_total <- read_csv(here("Data", "Netherlands_cohort_MF_calcs"))
Netherlands_F <- read_csv(here("Data", "Netherlands_cohort_F_calcs"))
Netherlands_M <- read_csv(here("Data", "Netherlands_cohort_M_calcs"))

Denmark_total <- read_csv(here("Data", "Denmark_cohort_MF_calcs"))
Denmark_F <- read_csv(here("Data", "Denmark_cohort_F_calcs"))
Denmark_M <- read_csv(here("Data", "Denmark_cohort_M_calcs"))

France_total <- read_csv(here("Data", "France_cohort_MF_calcs"))
France_F <- read_csv(here("Data", "France_cohort_F_calcs"))
France_M <- read_csv(here("Data", "France_cohort_M_calcs"))

#---- Selected life table subsets ----
#Add computations for cumulative survival conditional on survival to age 50
all_life_netherlands <- Netherlands_total %>% 
  filter(Year == "1920-1925" & Age %in% seq(50, 95, by = 5)) 
all_life_netherlands %<>% 
  mutate("cum_surv_cond50" = lx/all_life_netherlands$lx[1])

female_life_netherlands <- Netherlands_F %>% 
  filter(Year == "1920-1925" & Age %in% seq(50, 95, by = 5)) 
female_life_netherlands %<>%
  mutate("cum_surv_cond50" = lx/female_life_netherlands$lx[1])

male_life_netherlands <- Netherlands_M %>% 
  filter(Year == "1920-1925" & Age %in% seq(50, 95, by = 5)) 
male_life_netherlands %<>%
  mutate("cum_surv_cond50" = lx/male_life_netherlands$lx[1])

#---- Hazard ratios ----
#Netherlands
male_haz_netherlands <- Netherlands_M %>% 
  filter(Year == "1920-1925" & Age %in% seq(0, 95, by = 5)) %>% 
  dplyr::select("Year", "Haz")

female_haz_netherlands <- Netherlands_F %>% 
  filter(Year == "1920-1925" & Age %in% seq(0, 95, by = 5)) %>% 
  dplyr::select("Year", "Haz")

Hratio_Netherlands <- female_haz_netherlands$Haz/male_haz_netherlands$Haz %>%
  as.data.frame() %>% set_colnames("Netherlands") 

#Denmark
male_haz_denmark <- Denmark_M %>% 
  filter(Year == "1920-1925" & Age %in% seq(0, 95, by = 5)) %>% 
  dplyr::select("Year", "Haz")

female_haz_denmark <- Denmark_F %>% 
  filter(Year == "1920-1925" & Age %in% seq(0, 95, by = 5)) %>% 
  dplyr::select("Year", "Haz")

Hratio_Denmark <- female_haz_denmark$Haz/male_haz_denmark$Haz %>% 
  as.data.frame() %>% set_colnames("Denmark")

#France
male_haz_france <- France_M %>% 
  filter(Year == "1920-1925" & Age %in% seq(0, 95, by = 5)) %>% 
  dplyr::select("Year", "Haz")

female_haz_france <- France_F %>% 
  filter(Year == "1920-1925" & Age %in% seq(0, 95, by = 5)) %>% 
  dplyr::select("Year", "Haz")

Hratio_France <- female_haz_france$Haz/male_haz_france$Haz %>%
  as.data.frame() %>% set_colnames("France") 
