format_plot_data <- function(sample, scenario){
  plot_data <- sample %>% 
    dplyr::select(c("female", "U", "death0", 
                    head(variable_names$deathij_varnames, -1))) %>%
    mutate("Sex" = if_else(female == 0, "Men", "Women")) %>%
    mutate("Scenario" = scenario) %>%
    dplyr::select(-one_of("female")) %>% 
    dplyr::select("U", "death0", everything()) %>% 
    set_colnames(c("U", seq(1:10), "Sex/Gender", "Scenario")) %>% 
    gather(-contains(c("U", "Sex/Gender", "Scenario")), key = "Visit", 
           value = "death_indicator") %>%
    filter(death_indicator == 0) %>% 
    mutate_at("Sex/Gender", as.factor) %>% 
    mutate_at("Visit", as.factor)
  plot_data$Visit <- fct_relevel(plot_data$Visit, "10", after = Inf)
  
  return(plot_data)
}
