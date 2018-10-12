#---- Model for functional ability ----
func_ability <- function(knots_ages, slopes, obs){
  extend_slopes <- c(slopes[1], rep(0, num_tests - 1))
  ages <- visit_times + 50
  for(k in 1:length(knots_ages)){
    extend_slopes[which(ages == knots_ages[k])] <- slopes[k + 1]
  }
  mid_visits <- visit_times[c(-1, -length(visit_times))]
  test_nums = seq(from = 0, to = num_tests, by = 1)
  testXslope = (-1)*extend_slopes[-1]*mid_visits
  Fij <- vector(length = length(visit_times))
  for(j in 1:length(visit_times)){
    t = visit_times[j]
    test_num = test_nums[j]
    delta <- paste("delta", test_num, sep = "")
    if(ages[j] <= knots_ages[1]){
      Fij[j] = a00 + obs[, "w0i"] + a01*obs[, "age0_c50"] + 
        a02*obs[, "U"] + obs[, delta] +
        (extend_slopes[1] + obs[, "w1i"] + a11*obs[, "age0_c50"] + 
           a12*obs[, "U"])*t
    } else{
      Fij[j] = a00 + obs[, "w0i"] + a01*obs[, "age0_c50"] + 
        a02*obs[, "U"] + obs[, delta] +
        sum(testXslope[1:(test_num - 1)]) +
        (sum(extend_slopes[1:test_num]) + obs[, "w1i"] + 
           a11*obs[, "age0_c50"] + a12*obs[, "U"])*t
    }
  }
  slopes <- matrix(NA, nrow = num_obs, ncol= (length(visit_times) - 1))
  #Fij is stored as a list in this function environment so use list indexing
  for(j in 1:ncol(slopes)){
    b <- Fij[[j + 1]]
    a <- Fij[[j]]
    slopes[, j] = (b-a)/int_time
  }
  return(list("Fij" = Fij, "slopes" = slopes))
}