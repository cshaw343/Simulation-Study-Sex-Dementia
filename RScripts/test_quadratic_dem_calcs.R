#---- Package loading, options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "MASS", "reshape", "magrittr", "here")

options(scipen = 999) #Standard Notation
options(digits = 6)   #Round to 6 decimal places
options(warn = -1)    #Suppress warnings

#---- Source files ----
source(here("RScripts", "sex-dementia_sim_parA_onedemcut_nodemkill_UonInt.R"))
source(here("RScripts", "sex-dementia_sim_data_gen.R")) 

#---- Create the data ----
test_data <- data_gen(100000)

#---- Create extra data columns ----
test_data %<>% mutate("a0" = NA, "a1" = NA, "a2" = NA) 
test_data %<>% mutate("quad_dem_onset" = NA)

#---- Fit the quadratics ----
Cij_data <- test_data %>% dplyr::select(variable_names$Cij_varnames) %>% 
  t() %>% as.data.frame() %>% 
  mutate("c_age" = seq(0, 45, by = 5), 
         "c_age2" = c_age^2) 

for(i in 1:nrow(test_data)){
  test_data[i, c("a0", "a1", "a2")] <- 
    lm(Cij_data[, i] ~ c_age + c_age2, data = Cij_data)$coefficients
}  

#---- Solve for dementia onset time ----
quad_coefficients <- test_data %>% dplyr::select(c("a0", "a1", "a2")) %>% 
  t() %>% as.data.frame() 

quad_coefficients[is.na(quad_coefficients)] <- 0

for(i in 1:nrow(test_data)){
  
  quad_function <- function(x){
    abs(quad_coefficients["a0", i] + quad_coefficients["a1", i]*x + 
      quad_coefficients["a2", i]*x^2 - dem_cut)
  }
  
  test_data[i, "quad_dem_onset"] <- 
    optimize(quad_function, lower = 0, upper = 50)
}

#---- Compare dementia onset times (plot) ----
ggplot(test_data) + 
  geom_point(aes(id, timetodem)) + 
  geom_point(aes(id, quad_dem_onset), color = "pink") + 
  theme_minimal()

#---- Count incident cases/year ----
test_data %<>% mutate("floor_quad_onset" = floor(quad_dem_onset)) 

annual_counts <- test_data %>% group_by(floor_quad_onset) %>% 
  tally()

plot(annual_counts)

ggplot(aes(, n), data = annual_counts) + 
  geom_point(color = "#00BFC4", size = 2) + theme_minimal() + 
  xlab("Year") + ylab("Incident Dementia Cases") + 
  ggtitle("Incident Dementia Cases by Year") + 
  scale_x_continuous(breaks = seq(0, 45, by = 5))

#---- Compute regular inc rates/1000 PY ----
#---- Dementia incidence rates ----
reg_sim_rates <- vector(length = num_tests) 

#Computing incidence rates
for(slot in 1:num_tests){
  if(slot == 1){
    dem_last_wave <- paste0("dem", (slot - 1))
    dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
    contributed <- paste0("contributed", (slot - 1), "-", slot)
  } else {
    dem_last_wave <- paste0("dem", (slot - 2), "-", (slot - 1))
    dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
    contributed <- paste0("contributed", (slot - 1), "-", slot)
  }
  PY_data <- test_data %>% 
    dplyr::select(dem_last_wave, dem_this_wave, contributed) %>% 
    filter(!! as.name(dem_last_wave) == 0) 
  
  reg_sim_rates[slot] = 
    round(1000*(sum(PY_data[, dem_this_wave], na.rm = TRUE)/
                  sum(PY_data[, contributed], na.rm = TRUE)), 3)
}

#---- Dementia incidence cases, rates, and PY by sex (5-year bands) ----
female_data <- test_data %>% filter(female == 1)
male_data <- test_data %>% filter(female == 0)

reg_sim_rates_females <- vector(length = num_tests) 
reg_sim_rates_males <- vector(length = num_tests) 

reg_inc_cases_females <- vector(length = num_tests) 
reg_inc_cases_males <- vector(length = num_tests) 

reg_PY_females <- vector(length = num_tests) 
reg_PY_males <- vector(length = num_tests) 

#Computing female incidence cases, rates, PY
for(slot in 1:num_tests){
  if(slot == 1){
    dem_last_wave <- paste0("dem", (slot - 1))
    dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
    contributed <- paste0("contributed", (slot - 1), "-", slot)
  } else {
    dem_last_wave <- paste0("dem", (slot - 2), "-", (slot - 1))
    dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
    contributed <- paste0("contributed", (slot - 1), "-", slot)
  }
  PY_data <- female_data %>% 
    dplyr::select(dem_last_wave, dem_this_wave, contributed) %>% 
    filter(!! as.name(dem_last_wave) == 0) 
  
  reg_inc_cases_females[slot] = sum(PY_data[, dem_this_wave], na.rm = TRUE)
  
  reg_PY_females[slot] = sum(PY_data[, contributed], na.rm = TRUE)
  
  reg_sim_rates_females[slot] = 
    round(1000*(sum(PY_data[, dem_this_wave], na.rm = TRUE)/
                  sum(PY_data[, contributed], na.rm = TRUE)), 3)
}

#Computing male incidence cases and rates
for(slot in 1:num_tests){
  if(slot == 1){
    dem_last_wave <- paste0("dem", (slot - 1))
    dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
    contributed <- paste0("contributed", (slot - 1), "-", slot)
  } else {
    dem_last_wave <- paste0("dem", (slot - 2), "-", (slot - 1))
    dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
    contributed <- paste0("contributed", (slot - 1), "-", slot)
  }
  PY_data <- male_data %>% 
    dplyr::select(dem_last_wave, dem_this_wave, contributed) %>% 
    filter(!! as.name(dem_last_wave) == 0)
  
  reg_inc_cases_males[slot] = sum(PY_data[, dem_this_wave], na.rm = TRUE)
  
  reg_PY_males[slot] = sum(PY_data[, contributed], na.rm = TRUE)
  
  reg_sim_rates_males[slot] = 
    round(1000*sum(PY_data[, dem_this_wave], na.rm = TRUE)/
            sum(PY_data[, contributed], na.rm = TRUE), 3)
}

#---- Dementia logIRRs (5-year bands) ----
IRRs <- reg_sim_rates_females/reg_sim_rates_males

#---- Dementia incidence cases, rates, and PY by sex (1-year bands) ----
reg_sim_rates_females_1year <- vector(length = num_tests*5) 
reg_sim_rates_males_1year <- vector(length = num_tests*5) 

reg_inc_cases_females_1year <- vector(length = num_tests*5)
reg_inc_cases_males_1year <- vector(length = num_tests*5)

reg_PY_females_1year <- vector(length = num_tests*5)
reg_PY_males_1year <- vector(length = num_tests*5)

#Computing female incidence cases, rates, PY
for(slot in 1:(num_tests*5)){
  if(slot == 1){
    dem_last_wave <- paste0("dem", (slot - 1))
    dem_this_wave <- paste0("dem", ((slot - 1) + 50), "-", (slot + 50))
    contributed <- paste0("contributed", ((slot - 1) + 50), "-", (slot + 50))
  } else {
    dem_last_wave <- paste0("dem", ((slot - 2) + 50), "-", ((slot - 1) + 50))
    dem_this_wave <- paste0("dem", ((slot - 1) + 50), "-", (slot + 50))
    contributed <- paste0("contributed", ((slot - 1) + 50), "-", (slot + 50))
  }
  PY_data <- female_data %>%
    dplyr::select(dem_last_wave, dem_this_wave, contributed) %>%
    filter(!! as.name(dem_last_wave) == 0)
  
  reg_inc_cases_females_1year[slot] = sum(PY_data[, dem_this_wave], 
                                          na.rm = TRUE)
  
  reg_PY_females_1year[slot] = sum(PY_data[, contributed], na.rm = TRUE)
}

reg_sim_rates_females_1year <-
  round(1000*(reg_inc_cases_females_1year/reg_PY_females_1year), 3)


#Computing male incidence cases, rates, PY
for(slot in 1:(num_tests*5)){
  if(slot == 1){
    dem_last_wave <- paste0("dem", (slot - 1))
    dem_this_wave <- paste0("dem", ((slot - 1) + 50), "-", (slot + 50))
    contributed <- paste0("contributed", ((slot - 1) + 50), "-", (slot + 50))
  } else {
    dem_last_wave <- paste0("dem", ((slot - 2) + 50), "-", ((slot - 1) + 50))
    dem_this_wave <- paste0("dem", ((slot - 1) + 50), "-", (slot + 50))
    contributed <- paste0("contributed", ((slot - 1) + 50), "-", (slot + 50))
  }
  PY_data <- male_data %>%
    dplyr::select(dem_last_wave, dem_this_wave, contributed) %>%
    filter(!! as.name(dem_last_wave) == 0)
  
  reg_inc_cases_males_1year[slot] = sum(PY_data[, dem_this_wave], na.rm = TRUE)
  
  reg_PY_males_1year[slot] = sum(PY_data[, contributed], na.rm = TRUE)
}

reg_sim_rates_males_1year <-
  round(1000*(reg_inc_cases_males_1year/reg_PY_males_1year), 3)

#---- Dementia logIRRs (1-year bands) ----
reg_IRRs_1year <- reg_sim_rates_females_1year/reg_sim_rates_males_1year

#---- Fix data to do quadratic inc rates/1000 PY Calculation ----
#---- Compute quad inc rates/1000 PY ----
#---- Dementia incidence rates ----
quad_sim_rates <- vector(length = num_tests) 

#Computing incidence rates
for(slot in 1:num_tests){
  if(slot == 1){
    dem_last_wave <- paste0("dem", (slot - 1))
    dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
    contributed <- paste0("contributed", (slot - 1), "-", slot)
  } else {
    dem_last_wave <- paste0("dem", (slot - 2), "-", (slot - 1))
    dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
    contributed <- paste0("contributed", (slot - 1), "-", slot)
  }
  PY_data <- test_data %>% 
    dplyr::select(dem_last_wave, dem_this_wave, contributed) %>% 
    filter(!! as.name(dem_last_wave) == 0) 
  
  reg_sim_rates[slot] = 
    round(1000*(sum(PY_data[, dem_this_wave], na.rm = TRUE)/
                  sum(PY_data[, contributed], na.rm = TRUE)), 3)
}

#---- Dementia incidence cases, rates, and PY by sex (5-year bands) ----
female_data <- test_data %>% filter(female == 1)
male_data <- test_data %>% filter(female == 0)

reg_sim_rates_females <- vector(length = num_tests) 
reg_sim_rates_males <- vector(length = num_tests) 

reg_inc_cases_females <- vector(length = num_tests) 
reg_inc_cases_males <- vector(length = num_tests) 

reg_PY_females <- vector(length = num_tests) 
reg_PY_males <- vector(length = num_tests) 

#Computing female incidence cases, rates, PY
for(slot in 1:num_tests){
  if(slot == 1){
    dem_last_wave <- paste0("dem", (slot - 1))
    dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
    contributed <- paste0("contributed", (slot - 1), "-", slot)
  } else {
    dem_last_wave <- paste0("dem", (slot - 2), "-", (slot - 1))
    dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
    contributed <- paste0("contributed", (slot - 1), "-", slot)
  }
  PY_data <- female_data %>% 
    dplyr::select(dem_last_wave, dem_this_wave, contributed) %>% 
    filter(!! as.name(dem_last_wave) == 0) 
  
  reg_inc_cases_females[slot] = sum(PY_data[, dem_this_wave], na.rm = TRUE)
  
  reg_PY_females[slot] = sum(PY_data[, contributed], na.rm = TRUE)
  
  reg_sim_rates_females[slot] = 
    round(1000*(sum(PY_data[, dem_this_wave], na.rm = TRUE)/
                  sum(PY_data[, contributed], na.rm = TRUE)), 3)
}

#Computing male incidence cases and rates
for(slot in 1:num_tests){
  if(slot == 1){
    dem_last_wave <- paste0("dem", (slot - 1))
    dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
    contributed <- paste0("contributed", (slot - 1), "-", slot)
  } else {
    dem_last_wave <- paste0("dem", (slot - 2), "-", (slot - 1))
    dem_this_wave <- paste0("dem", (slot - 1), "-", slot)
    contributed <- paste0("contributed", (slot - 1), "-", slot)
  }
  PY_data <- male_data %>% 
    dplyr::select(dem_last_wave, dem_this_wave, contributed) %>% 
    filter(!! as.name(dem_last_wave) == 0)
  
  reg_inc_cases_males[slot] = sum(PY_data[, dem_this_wave], na.rm = TRUE)
  
  reg_PY_males[slot] = sum(PY_data[, contributed], na.rm = TRUE)
  
  reg_sim_rates_males[slot] = 
    round(1000*sum(PY_data[, dem_this_wave], na.rm = TRUE)/
            sum(PY_data[, contributed], na.rm = TRUE), 3)
}

#---- Dementia logIRRs (5-year bands) ----
IRRs <- reg_sim_rates_females/reg_sim_rates_males

#---- Dementia incidence cases, rates, and PY by sex (1-year bands) ----
reg_sim_rates_females_1year <- vector(length = num_tests*5) 
reg_sim_rates_males_1year <- vector(length = num_tests*5) 

reg_inc_cases_females_1year <- vector(length = num_tests*5)
reg_inc_cases_males_1year <- vector(length = num_tests*5)

reg_PY_females_1year <- vector(length = num_tests*5)
reg_PY_males_1year <- vector(length = num_tests*5)

#Computing female incidence cases, rates, PY
for(slot in 1:(num_tests*5)){
  if(slot == 1){
    dem_last_wave <- paste0("dem", (slot - 1))
    dem_this_wave <- paste0("dem", ((slot - 1) + 50), "-", (slot + 50))
    contributed <- paste0("contributed", ((slot - 1) + 50), "-", (slot + 50))
  } else {
    dem_last_wave <- paste0("dem", ((slot - 2) + 50), "-", ((slot - 1) + 50))
    dem_this_wave <- paste0("dem", ((slot - 1) + 50), "-", (slot + 50))
    contributed <- paste0("contributed", ((slot - 1) + 50), "-", (slot + 50))
  }
  PY_data <- female_data %>%
    dplyr::select(dem_last_wave, dem_this_wave, contributed) %>%
    filter(!! as.name(dem_last_wave) == 0)
  
  reg_inc_cases_females_1year[slot] = sum(PY_data[, dem_this_wave], 
                                          na.rm = TRUE)
  
  reg_PY_females_1year[slot] = sum(PY_data[, contributed], na.rm = TRUE)
}

reg_sim_rates_females_1year <-
  round(1000*(reg_inc_cases_females_1year/reg_PY_females_1year), 3)


#Computing male incidence cases, rates, PY
for(slot in 1:(num_tests*5)){
  if(slot == 1){
    dem_last_wave <- paste0("dem", (slot - 1))
    dem_this_wave <- paste0("dem", ((slot - 1) + 50), "-", (slot + 50))
    contributed <- paste0("contributed", ((slot - 1) + 50), "-", (slot + 50))
  } else {
    dem_last_wave <- paste0("dem", ((slot - 2) + 50), "-", ((slot - 1) + 50))
    dem_this_wave <- paste0("dem", ((slot - 1) + 50), "-", (slot + 50))
    contributed <- paste0("contributed", ((slot - 1) + 50), "-", (slot + 50))
  }
  PY_data <- male_data %>%
    dplyr::select(dem_last_wave, dem_this_wave, contributed) %>%
    filter(!! as.name(dem_last_wave) == 0)
  
  reg_inc_cases_males_1year[slot] = sum(PY_data[, dem_this_wave], na.rm = TRUE)
  
  reg_PY_males_1year[slot] = sum(PY_data[, contributed], na.rm = TRUE)
}

reg_sim_rates_males_1year <-
  round(1000*(reg_inc_cases_males_1year/reg_PY_males_1year), 3)

#---- Dementia logIRRs (1-year bands) ----
IRRs_1year <- reg_sim_rates_females_1year/reg_sim_rates_males_1year








