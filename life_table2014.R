#***************************************************************
# 2014 Life Table
# Source: National Vital Statistics Reports, Vol. 66, No. 4, 
# August 14, 2017 (pg 48-49)
# Birth cohort:  1919 - 1921
#***************************************************************

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "ggplot2")

#---- Conditional Probabilities Function ----
cond_prob <- function(x){
  probs <- vector(length = length(x))
  for(i in 2:length(probs)){
    probs[i] = x[i]/x[i - 1]
  }
  probs[1] <- 1
  return(probs)
}

#---- Hazard Function ----
haz <- function(age, logprobs){
  HZ <- vector(length = length(age))
  for(i in 2:length(HZ)){
    HZ[i] = -(logprobs[i] - logprobs[i - 1])/(age[i] - age[i - 1])
  }
  return(HZ)
}


#---- Life Table Data ----
#"Survivors" represents number surviving out of 100,000 born alive
ages <- seq(from = 45, to = 100, by = 5) #Start at 45 for hazard calcs

life <- tibble("Age" = ages, 
               "Survivors" = c(72036, 68429, 63947, 58079, 50560, 41090, 29729, 
                               18298, 8683, 2941, 646, 67)) %>% 
  mutate("Prob" = Survivors/100000,
         "logProb" = log(Prob),
         "CP" = cond_prob(Survivors), 
         "Haz" = haz(age = Age, logprobs = logProb))

male_life <- tibble("MAge" = ages, 
                    "MSurvivors" = c(71244, 67553, 62965, 56917, 49218, 39668, 
                                     28316, 17128, 7920, 2527, 556, 62)) %>% 
  mutate("MProb" = MSurvivors/100000,
         "MlogProb" = log(MProb),
         "MCP" = cond_prob(MSurvivors), 
         "MHaz" = haz(age = MAge, logprobs = MlogProb))

female_life <- tibble("FAge" = ages, 
                      "FSurvivors" = c(72954, 69452, 65099, 59438, 52126, 42741, 
                                       31344, 19613, 9515, 3314, 728, 72)) %>% 
  mutate("FProb" = FSurvivors/100000,
         "FlogProb" = log(FProb),
         "FCP" = cond_prob(FSurvivors), 
         "FHaz" = haz(age = FAge, logprobs = FlogProb))

Hratio <- male_life$MHaz/female_life$FHaz %>% as.data.frame()
colnames(Hratio) <- c("ratio")

#---- Hazard Plots ----
#Creating plot data
female_hazards <- female_life %>% dplyr::select(c("FAge", "FHaz")) %>%
  mutate("Age" = FAge) %>% dplyr::select(-FAge) %>% melt(., id.vars = "Age")

male_hazards <- male_life %>% dplyr::select(c("MAge", "MHaz")) %>%
  mutate("Age" = MAge) %>% dplyr::select(-MAge) %>% melt(., id.vars = "Age")

haz_ratios <- Hratio %>% cbind(ages) %>% 
  mutate("Age" = ages) %>% dplyr::select(-ages) %>% 
  melt(., id.vars = "Age")

haz_plot_data <- rbind(female_hazards, male_hazards, haz_ratios)

#Creating plot
hazard_plot<- ggplot(haz_plot_data, aes(Age, value), color = variable) + 
  geom_line(data = subset(haz_plot_data, variable == "FHaz"), 
            aes(color = variable), size = 1) + 
  geom_line(data = subset(haz_plot_data, variable == "MHaz"), 
            aes(color = variable), size = 1, alpha = 0.6) + 
  geom_line(data = subset(haz_plot_data, variable == "ratio"), 
            aes(color = variable), size = 1) + ylim(0, 1.25) + 
  labs(y = "Hazard", 
       x = "Age", 
       color = "Hazard") + 
  theme_minimal()

#Saving plot output
ggsave(filename = "hazard_plot_2014life.jpeg", width = 10, height = 7, 
       plot = hazard_plot)
