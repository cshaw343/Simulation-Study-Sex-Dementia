format_plot_data <- function(sample, scenario){
  plot_data <- sample %>% 
    dplyr::select(c("female", "U", "death0", 
                    head(variable_names$deathij_varnames, -1))) %>%
    mutate("Sex" = if_else(female == 0, "Men", "Women")) %>%
    mutate("Scenario" = scenario) %>%
    dplyr::select(-one_of("female")) %>% 
    dplyr::select("U", "death0", everything()) %>% 
    set_colnames(c("U", "Baseline", "[50-55)", "[55-60)","[60-65)", 
                   "[65-70)", "[70-75)","[75-80)", 
                   "[80-85)", "[85-90)","[90-95)", 
                   "Sex/Gender", "Scenario")) %>% 
    gather(-contains(c("U", "Sex/Gender", "Scenario")), key = "Age Band", 
           value = "death_indicator") %>%
    filter(death_indicator == 0) %>% 
    mutate_at("Sex/Gender", as.factor) %>% 
    mutate_at("Age Band", as.factor)
  plot_data$`Age Band` <- relevel(plot_data$`Age Band`, "Baseline")
  plot_data$`Age Band` <- fct_relevel(plot_data$`Age Band`, 
                                      rev(levels(plot_data$`Age Band`)))
  
  return(plot_data)
}
