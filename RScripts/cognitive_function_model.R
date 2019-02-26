#---- Model for Cognitive Function ----
cog_func <- function(knots_ages, slopes, obs, eps, noise, age_c){
  extend_slopes <- c(slopes[1], rep(0, num_tests - 1))
  ages <- visit_times + 50
  for(k in 1:length(knots_ages)){
    extend_slopes[which(ages == knots_ages[k])] <- slopes[k + 1]
  }
  mid_visits <- visit_times[c(-1, -length(visit_times))]
  test_nums = seq(from = 0, to = num_tests, by = 1)
  testXslope = (-1)*extend_slopes[-1]*mid_visits
  Cij <- vector(length = length(visit_times))
  for(j in 1:length(visit_times)){
    t = visit_times[j]
    test_num = test_nums[j]
    eps_name <- paste("eps", test_num, sep = "")
    z0_name <- paste0("z0_", (j - 1), "i")
    z1_name <- paste0("z1_", (j - 1), "i")
    if(ages[j] <= knots_ages[1]){
      Cij[j] = b00 + noise[, z0_name] + b01*obs[, "sex"] + 
        b02*age_c[, "age0_c50"] + b03*obs[, "U"] + eps[, eps_name] +
        (extend_slopes[1] + noise[, z1_name] + b11*obs[, "sex"] +
           b12*age_c[, "age0_c50"] + b13*obs[, "U"])*t
    } else{
      Cij[j] = b00 + noise[, z0_name] + b01*obs[, "sex"] +
        b02*age_c[, "age0_c50"] + b03*obs[, "U"] + eps[, eps_name] +
        sum(testXslope[1:(test_num - 1)]) +
        (sum(extend_slopes[1:test_num]) + noise[, z1_name] + b11*obs[, "sex"] +
           b12*age_c[, "age0_c50"] + b13*obs[, "U"])*t
    }
  }
  slopes <- matrix(NA, nrow = nrow(obs), ncol= (length(visit_times) - 1))
  #Cij is stored as a list in this function environment so use list indexing
  for(j in 1:ncol(slopes)){
    b <- Cij[[j + 1]]
    a <- Cij[[j]]
    slopes[, j] = (b-a)/int_time
  }
  return(list("Cij" = Cij, "slopes" = slopes))
}

