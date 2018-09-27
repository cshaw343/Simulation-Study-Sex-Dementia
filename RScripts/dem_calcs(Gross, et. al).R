#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse")

#Suppress warnings
options(warn = -1)

#Only consider the General cognitive performance data
#Relevant individuals aged 70-85
cog_decline <- tibble("WHICAP" = c(-0.4, -0.5, -0.4), 
                      "SENAS" = c(-0.8, -0.4, -0.4), 
                      "DUKE" = c(-0.4, -0.4, -0.45), 
                      "NCODE" = c(0, 0, 0.1)) %>% 
  map_df(.f = ~(.)/10) #Change scale to SD = 1 
rownames(cog_decline) <- c("White", "Black", "Hispanic")

#Average all data 
avg_decline <- cog_decline %>% map_dbl(mean) %>% mean(.)
