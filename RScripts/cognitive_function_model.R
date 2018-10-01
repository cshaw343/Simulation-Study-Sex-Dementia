#---- Model for Cognitive Function ----
cog_func <- function(slopes, obs){
  knots = visit_times
  mid_visits <- knots[c(-1, -length(knots))]
  test_nums = seq(from = 0, to = num_tests, by = 1)
  testXslope = (-1)*slopes[-1]*mid_visits
  Cij <- vector(length = length(visit_times))
  for(j in 1:length(Cij)){
    t = knots[j]
    test_num = test_nums[j]
    eps <- paste("eps", test_num, sep = "")
    if(t <= 5){
      Cij[j] = b00 + obs[, "z0i"] + b01*obs[, "sex"] +
        b02*obs[, "age0_c50"] + b03*obs[, "U"] + obs[, eps] +
        (slopes[1] + obs[, "z1i"] + b11*obs[, "sex"] +
           b12*obs[, "age0_c50"] + b13*obs[, "U"])*t
    } else{
      Cij[j] = b00 + obs[, "z0i"] + b01*obs[, "sex"] +
        b02*obs[, "age0_c50"] + b03*obs[, "U"] + obs[, eps] +
        sum(testXslope[1:(test_num - 1)]) +
        (sum(slopes[1:test_num]) + obs[, "z1i"] + b11*obs[, "sex"] +
           b12*obs[, "age0_c50"] + b13*obs[, "U"])*t
    }
  }
  slopes <- matrix(NA, nrow = num_obs, ncol= (length(visit_times) - 1))
  #Cij is stored as a list in this function environment so use list indexing
  for(j in 1:ncol(slopes)){
    b <- Cij[[j + 1]]
    a <- Cij[[j]]
    slopes[, j] = (b-a)/int_time
  }
  return(list("Cij" = Cij, "slopes" = slopes))
}

#Old cog_func Code
# knots = c(0, 20, 35)
# Cij <- vector(length = length(visit_times))
# for(j in 1:length(Cij)){
#   t = visit_times[j]
#   test_num = j - 1
#   eps <- paste("eps", test_num, sep = "")
#   Cij[j] = b00 + obs[, "z0i"] + b01*obs[, "sex"] + b02*obs[, "age0_c50"] + 
#     b03*obs[, "U"] + obs[, eps] + 
#     (b10a - b10b)*knots[2]*(t >= knots[2]) + 
#     (b10b - b10c)*knots[3]*(t >= knots[3]) + 
#     (obs[, "z1i"] + b11*obs[, "sex"] + b12*obs[, "age0_c50"] + 
#        b13*obs[, "U"] + b10a*(t >= knots[1] & t< knots[2]) + 
#        b10b*(t >= knots[2] & t< knots[3]) +
#        b10c*(t >= knots[3]))*t
#   }
#   slopes <- matrix(NA, nrow = num_obs, ncol= (length(visit_times) - 1))
#   #Cij is stored as a list in this function environment so use list indexing
#   for(j in 1:ncol(slopes)){
#     b <- Cij[[j + 1]] 
#     a <- Cij[[j]]
#     slopes[, j] = (b-a)/int_time
#   }
#   return(list("Cij" = Cij, "slopes" = slopes))
# }

