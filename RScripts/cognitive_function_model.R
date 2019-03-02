#---- Model for Cognitive Function ----
cog_func <- function(knots_ages, slopes, obs_matrix){
  extend_slopes <- c(slopes[1], rep(0, num_tests - 1))
  ages <- visit_times + 50
  for(k in 1:length(knots_ages)){
    extend_slopes[which(ages == knots_ages[k])] <- slopes[k + 1]
  }
  mid_visits <- visit_times[c(-1, -length(visit_times))]
  test_nums = seq(from = 0, to = num_tests, by = 1)
  testXslope = (-1)*extend_slopes[-1]*mid_visits
  Cij <- matrix(nrow = nrow(obs_matrix), ncol = length(visit_times))
  for(j in 1:ncol(Cij)){
    t = visit_times[j]
    test_num = test_nums[j]
    eps_name <- paste("eps", test_num, sep = "")
    z0_name <- paste0("z0_", (j - 1), "i")
    z1_name <- paste0("z1_", (j - 1), "i")
    if(ages[j] <= knots_ages[1]){
      Cij[, j] = b00 + obs_matrix[, z0_name] + 
        b01*obs_matrix[, "sex"] + b02*obs_matrix[, "age0_c50"] + 
        b03*obs_matrix[, "U"] + obs_matrix[, eps_name] +
        (extend_slopes[1] + obs_matrix[, z1_name] + 
           b11*obs_matrix[, "sex"] + b12*obs_matrix[, "age0_c50"] + 
           b13*obs_matrix[, "U"])*t
    } else{
      Cij[, j] = b00 + obs_matrix[, z0_name] + 
        b01*obs_matrix[, "sex"] + b02*obs_matrix[, "age0_c50"] + 
        b03*obs_matrix[, "U"] + obs_matrix[, eps_name] +
        sum(testXslope[1:(test_num - 1)]) +
        (sum(extend_slopes[1:test_num]) + obs_matrix[, z1_name] + 
           b11*obs_matrix[, "sex"] +
           b12*obs_matrix[, "age0_c50"] + b13*obs_matrix[, "U"])*t
    }
  }
  
  slopes <- matrix(NA, nrow = nrow(obs_matrix), ncol= (length(visit_times) - 1))
  
  for(i in 1:ncol(slopes)){
    slopes[, i] <- (Cij[, i + 1] - Cij[, i])/int_time
  }
  
  return(list("Cij" = Cij, "slopes" = slopes))
}
