#---- Model for Cognitive Function ----
cog_func <- function(knots_ages, slopes, obs){
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
    eps <- paste("eps", test_num, sep = "")
    if(ages[j] <= knots_ages[1]){
      Cij[j] = b00 + obs[, "z0i"] + b01*obs[, "sex"] + 
        b02*obs[, "age0_c50"] + b03*obs[, "U"] + obs[, eps] +
        (extend_slopes[1] + obs[, "z1i"] + b11*obs[, "sex"] +
           b12*obs[, "age0_c50"] + b13*obs[, "U"])*t
    } else{
      Cij[j] = b00 + obs[, "z0i"] + b01*obs[, "sex"] +
        b02*obs[, "age0_c50"] + b03*obs[, "U"] + obs[, eps] +
        sum(testXslope[1:(test_num - 1)]) +
        (sum(extend_slopes[1:test_num]) + obs[, "z1i"] + b11*obs[, "sex"] +
           b12*obs[, "age0_c50"] + b13*obs[, "U"])*t
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

