#---- Specify source file ----
source("sex-demensia_sim_parA.R")

#---- Data generation script ----
#Creating IDs, baseline age
#Creating centered ages
obs <- tibble("id" = seq(from = 1, to = num_obs, by = 1), 
              "age0" = age0,
              "sex" = rbinom(num_obs, size = 1, prob = psex)) %>% 
  mutate("age_c" = age0 - mean(age0))

#Check the proportion of males
pmales <- mean(obs$sex)


