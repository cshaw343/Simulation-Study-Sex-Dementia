results_65_plus <- function(results){
  results %<>% 
    dplyr::select(-one_of(variable_names$logIRR_varnames[1:3])) %>%
    dplyr::select(-one_of(variable_names$logIRR_95CI_Lower_varnames[1:3])) %>% 
    dplyr::select(-one_of(variable_names$logIRR_95CI_Upper_varnames[1:3])) %>% 
    dplyr::select(-one_of(variable_names$logIRR_95CI_Coverage_varnames[1:3])) %>% 
    dplyr::select(-one_of(variable_names$logIRR_SE_varnames[1:3])) %>% 
    dplyr::select(-one_of(variable_names$prop_dem_random_by_age[1:3])) %>%
    dplyr::select(-one_of(variable_names$prop_dem_Cij_by_age[1:3])) %>%
    dplyr::select(-one_of(variable_names$prop_dem_both_by_age[1:3])) %>%
    dplyr::select(-one_of(variable_names$prop_dem_random_W_by_age[1:3])) %>%
    dplyr::select(-one_of(variable_names$prop_dem_Cij_W_by_age[1:3])) %>%
    dplyr::select(-one_of(variable_names$prop_dem_both_W_by_age[1:3])) %>%
    dplyr::select(-one_of(variable_names$prop_dem_random_M_by_age[1:3])) %>%
    dplyr::select(-one_of(variable_names$prop_dem_Cij_M_by_age[1:3])) %>%
    dplyr::select(-one_of(variable_names$prop_dem_both_M_by_age[1:3])) %>%
    dplyr::select(-one_of(variable_names_1year$logIRR_varnames[1:15])) %>%
    dplyr::select(
      -one_of(variable_names_1year$logIRR_95CI_Lower_varnames[1:15])) %>% 
    dplyr::select(
      -one_of(variable_names_1year$logIRR_95CI_Upper_varnames[1:15])) %>% 
    dplyr::select(
      -one_of(variable_names_1year$logIRR_95CI_Coverage_varnames[1:15])) %>% 
    dplyr::select(-one_of(variable_names_1year$logIRR_SE_varnames[1:15]))
  
  return(results)
}
