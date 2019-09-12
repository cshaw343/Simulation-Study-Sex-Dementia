calc_coeff <- function(obs){
  a0 <- b00 + obs[, "z_0i"] + b01*obs[, "female"] + b02*obs[, "age0_c50"] +
    b03*obs[, "U"] + obs[, "eps_ij"]
  a1 <- b10 + obs[, "z_1i"] + b11*obs[, "female"] + b12*obs[, "age0_c50"] + 
    b13*obs[, "U"] 
  a2 <- b20 + obs[, "z_2i"] + b21*obs[, "female"] + b22*obs[, "age0_c50"] + 
    b23*obs[, "U"] 
  
  return(cbind(a0, a1, a2))
}
