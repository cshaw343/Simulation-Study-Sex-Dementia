#*******************************************************************************
# Life tables from mortality.org
#*******************************************************************************

#---- Package loading, options, seed ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "magrittr")

#---- Source Files ----
source("RScripts/misc_custom_functions.R")
source("RScripts/life_table2014.R") #For the purposes of comparing Euro with US data

#---- Hazard Function ----
haz <- function(age, logprobs){
  HZ <- vector(length = length(age))
  for(i in 2:length(HZ)){
    HZ[i] = -(logprobs[i] - logprobs[i - 1])/(age[i] - age[i - 1])
  }
  return(HZ)
}

#---- Read in and format life tables ----
#There was no cohort data for the UK

Netherlands_total <- read.table("Data/Netherlands_cohort_table_MF", skip = 2, 
                                header = TRUE) %>% 
  mutate("Sex" = "both") %>%            #Add Sex Indicator
  mutate("Country" = "Netherlands") %>% #Add Country Indicator
  mutate_at("Age", as.character) %>%    #Make Age numeric
  mutate_at("Age", as.numeric) %>% 
  filter(Age %% 5 == 0)                 #Keep Ages in 5 year increments
  
Netherlands_M <- read.table("Data/Netherlands_cohort_table_M", skip = 2, 
                            header = TRUE) %>% 
  mutate("Sex" = "Male") %>%            #Add Sex Indicator
  mutate("Country" = "Netherlands") %>% #Add Country Indicator
  mutate_at("Age", as.character) %>%    #Make Age numeric
  mutate_at("Age", as.numeric) %>% 
  filter(Age %% 5 == 0)                 #Keep Ages in 5 year increments

Netherlands_F <- read.table("Data/Netherlands_cohort_table_F", skip = 2, 
                                             header = TRUE) %>% 
  mutate("Sex" = "Female") %>%          #Add Sex Indicator
  mutate("Country" = "Netherlands") %>% #Add Country Indicator
  mutate_at("Age", as.character) %>%    #Make Age numeric
  mutate_at("Age", as.numeric) %>% 
  filter(Age %% 5 == 0)                 #Keep Ages in 5 year increments

Denmark_total <- read.table("Data/Denmark_cohort_table_MF", skip = 2, 
                                header = TRUE) %>% 
  mutate("Sex" = "both") %>%          #Add Sex Indicator
  mutate("Country" = "Denmark") %>%   #Add Country Indicator
  mutate_at("Age", as.character) %>%  #Make Age numeric
  mutate_at("Age", as.numeric) %>% 
  filter(Age %% 5 == 0)               #Keep Ages in 5 year increments

Denmark_M <- read.table("Data/Denmark_cohort_table_M", skip = 2, 
                            header = TRUE) %>% 
  mutate("Sex" = "Male") %>%          #Add Sex Indicator
  mutate("Country" = "Denmark") %>%   #Add Country Indicator
  mutate_at("Age", as.character) %>%  #Make Age numeric
  mutate_at("Age", as.numeric) %>% 
  filter(Age %% 5 == 0)               #Keep Ages in 5 year increments

Denmark_F <- read.table("Data/Denmark_cohort_table_F", skip = 2, 
                            header = TRUE) %>% 
  mutate("Sex" = "Female") %>%          #Add Sex Indicator
  mutate("Country" = "Denmark") %>%     #Add Country Indicator
  mutate_at("Age", as.character) %>%    #Make Age numeric
  mutate_at("Age", as.numeric) %>% 
  filter(Age %% 5 == 0)                 #Keep Ages in 5 year increments

France_total <- read.table("Data/France_cohort_table_MF", skip = 2, 
                                header = TRUE) %>% 
  mutate("Sex" = "both") %>%          #Add Sex Indicator
  mutate("Country" = "France") %>%    #Add Country Indicator
  mutate_at("Age", as.character) %>%  #Make Age numeric
  mutate_at("Age", as.numeric) %>% 
  filter(Age %% 5 == 0)               #Keep Ages in 5 year increments

France_M <- read.table("Data/France_cohort_table_M", skip = 2, 
                            header = TRUE) %>% 
  mutate("Sex" = "Male") %>%          #Add Sex Indicator
  mutate("Country" = "France") %>%    #Add Country Indicator
  mutate_at("Age", as.character) %>%  #Make Age numeric
  mutate_at("Age", as.numeric) %>% 
  filter(Age %% 5 == 0)               #Keep Ages in 5 year increments

France_F <- read.table("Data/France_cohort_table_F", skip = 2, 
                            header = TRUE) %>% 
  mutate("Sex" = "Female") %>%          #Add Sex Indicator
  mutate("Country" = "France") %>%      #Add Country Indicator
  mutate_at("Age", as.character) %>%    #Make Age numeric
  mutate_at("Age", as.numeric) %>% 
  filter(Age %% 5 == 0)                 #Keep Ages in 5 year increments

#---- Cohorts Available ----
#All of them have 1910-1919, but only Netherlands and Denmark have 1920-1925
table(Netherlands_total$Year)
table(Denmark_total$Year)
table(France_total$Year)

#---- Add Computations ----
Netherlands_total %<>% group_by(Year) %>%
  mutate("Prob" = lx/100000, 
         "logProb" = log(Prob), 
         "CP" = cond_prob(lx), 
         "Haz" = haz(age = Age, logprobs = logProb))

Netherlands_M %<>% group_by(Year) %>%
  mutate("Prob" = lx/100000, 
         "logProb" = log(Prob), 
         "CP" = cond_prob(lx), 
         "Haz" = haz(age = Age, logprobs = logProb))

Netherlands_F %<>% group_by(Year) %>%
  mutate("Prob" = lx/100000, 
         "logProb" = log(Prob), 
         "CP" = cond_prob(lx), 
         "Haz" = haz(age = Age, logprobs = logProb))

Denmark_total %<>% group_by(Year) %>%
  mutate("Prob" = lx/100000, 
         "logProb" = log(Prob), 
         "CP" = cond_prob(lx), 
         "Haz" = haz(age = Age, logprobs = logProb))

Denmark_M %<>% group_by(Year) %>%
  mutate("Prob" = lx/100000, 
         "logProb" = log(Prob), 
         "CP" = cond_prob(lx), 
         "Haz" = haz(age = Age, logprobs = logProb))

Denmark_F %<>% group_by(Year) %>%
  mutate("Prob" = lx/100000, 
         "logProb" = log(Prob), 
         "CP" = cond_prob(lx), 
         "Haz" = haz(age = Age, logprobs = logProb))
  
France_total %<>% group_by(Year) %>%
  mutate("Prob" = lx/100000, 
         "logProb" = log(Prob), 
         "CP" = cond_prob(lx), 
         "Haz" = haz(age = Age, logprobs = logProb))

France_M %<>% group_by(Year) %>%
  mutate("Prob" = lx/100000, 
         "logProb" = log(Prob), 
         "CP" = cond_prob(lx), 
         "Haz" = haz(age = Age, logprobs = logProb))

France_F %<>% group_by(Year) %>%
  mutate("Prob" = lx/100000, 
         "logProb" = log(Prob), 
         "CP" = cond_prob(lx), 
         "Haz" = haz(age = Age, logprobs = logProb))

#---- Combine tables into one dataset ----
all_tables <- rbind(Netherlands_total, Netherlands_M, Netherlands_F, 
                    Denmark_total, Denmark_M, Denmark_F, 
                    France_total, France_M, France_F)

#---- Making Plots ----
#Comparing countries for 1910-1919 birth cohort, facet by sex
all_countries_1910_1919 <- ggplot(all_tables, aes(Age, CP)) + 
  geom_line(data = subset(all_tables, Year == "1910-1919" & 
                            Age %in% seq(45, 100, by = 5)), 
            aes(color = Country, group = Country), size = 1.25,
            alpha = 0.6) + ylim(0, 1) + 
  labs(y = "Conditional Probability of Survival", x = "Age", 
       color = "") + theme_minimal() + facet_wrap(~Sex) + 
  ggtitle("1910-1919 Birth Cohort")

ggsave(filename = "Plots/1910-1919_birth_cohort.jpeg", width = 10, height = 7, 
       plot = all_countries_1910_1919)

#Comparing countries for 1920-1925 birth cohort, facet by sex
all_countries_1920_1925 <- ggplot(all_tables, aes(Age, CP)) + 
  geom_line(data = subset(all_tables, Year == "1920-1925" & 
                            Age %in% seq(45, 100, by = 5)), 
            aes(color = Country, group = Country), size = 1.25,
            alpha = 0.6) + ylim(0, 1) + 
  labs(y = "Conditional Probability of Survival", x = "Age", 
       color = "") + theme_minimal() + facet_wrap(~Sex) + 
  ggtitle("1920-1925 Birth Cohort")

ggsave(filename = "Plots/1920-1925_birth_cohort.jpeg", width = 10, height = 7, 
       plot = all_countries_1920_1925)

#Comparing countries for 1910-1919 birth cohort
all_countries_1920_1925_combined <- ggplot(all_tables, aes(Age, CP)) + 
  geom_line(data = subset(all_tables, Year == "1920-1925" & Sex == "Male"), 
            aes(color = Country, group = Country), 
            size = 1.25, linetype = "longdash", alpha = 0.6) + 
  geom_line(data = subset(all_tables, Year == "1910-1919" & Sex == "Female"), 
            aes(color = Country, group = Country), size = 1.25,
            alpha = 0.6) + 
  ylim(0, 1) + 
  labs(y = "Conditional Probability of Survival", x = "Age", 
       color = "") + theme_minimal() + 
  ggtitle("1920-1925 Birth Cohort")

ggsave(filename = "Plots/1920-1925_birth_cohort_combined.jpeg", width = 10, 
       height = 7, plot = all_countries_1920_1925_combined)

#---- Selected life table subsets ----
all_life_netherlands <- Netherlands_total %>% 
  filter(Year == "1920-1925" & Age %in% seq(50, 100, by = 5))

female_life_netherlands <- Netherlands_F %>% 
  filter(Year == "1920-1925" & Age %in% seq(50, 100, by = 5))

male_life_netherlands <- Netherlands_M %>% 
  filter(Year == "1920-1925" & Age %in% seq(50, 100, by = 5))

#---- Hazard ratios ----
#Netherlands
male_haz_netherlands <- Netherlands_M %>% 
  filter(Year == "1920-1925" & Age %in% seq(0, 100, by = 5)) %>% 
  dplyr::select("Haz")

female_haz_netherlands <- Netherlands_F %>% 
  filter(Year == "1920-1925" & Age %in% seq(0, 100, by = 5)) %>% 
  dplyr::select("Haz")

Hratio_Netherlands <- male_haz_netherlands$Haz/female_haz_netherlands$Haz %>%
  as.data.frame() %>% set_colnames("Netherlands_HR") 

#Denmark
male_haz_denmark <- Denmark_M %>% 
  filter(Year == "1920-1925" & Age %in% seq(0, 100, by = 5)) %>% 
  dplyr::select("Haz")

female_haz_denmark <- Denmark_F %>% 
  filter(Year == "1920-1925" & Age %in% seq(0, 100, by = 5)) %>% 
  dplyr::select("Haz")

Hratio_Denmark <- male_haz_denmark$Haz/female_haz_denmark$Haz %>% 
  as.data.frame() %>% set_colnames("Denmark_HR")

#France
male_haz_france <- France_M %>% 
  filter(Year == "1920-1925" & Age %in% seq(0, 100, by = 5)) %>% 
  dplyr::select("Haz")

female_haz_france <- France_F %>% 
  filter(Year == "1920-1925" & Age %in% seq(0, 100, by = 5)) %>% 
  dplyr::select("Haz")

Hratio_France <- male_haz_france$Haz/female_haz_france$Haz %>%
  as.data.frame() %>% set_colnames("France_HR") 

#Format US HR data
Hratio_US <- Hratio$ratio %>% as.data.frame() %>% set_colnames("US_HR")

#Combine hazard ratios
HR_plot_data <- 
  cbind(Hratio_Netherlands, Hratio_Denmark, Hratio_France, Hratio_US) %>% 
  mutate("Age" = seq(0, 100, by = 5)) %>%
  gather(contains("HR"), key = "Country", value = "HR")
  
#Plot Hazard Ratios 
all_countries_1920_1925_HR <- ggplot(HR_plot_data, aes(Age, HR)) + 
  geom_line(aes(color = Country, group = Country), size = 1.25, alpha = 0.6) +
  labs(y = "Hazard Ratio (Male:Female)", x = "Age", 
       color = "") + theme_minimal() + 
  ggtitle("1920-1925 Birth Cohort")

ggsave(filename = "Plots/1920-1925_birth_cohort_HR.jpeg", width = 10, 
       height = 7, plot = all_countries_1920_1925_HR)









