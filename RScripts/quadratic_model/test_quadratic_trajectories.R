#---- Package loading, options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "MASS", "reshape", "magrittr", "here")

options(scipen = 999) #Standard Notation
options(digits = 6)   #Round to 6 decimal places
options(warn = -1)    #Suppress warnings

#---- Source files ----
source(here("RScripts", "sex-dementia_sim_parA_onedemcut_nodemkill_UonInt.R"))
source(here("RScripts", "sex-dementia_sim_data_gen.R")) 

#---- Generate data ----
test_data <- data_gen(1000)

#---- Generate data before survival ----
#set size for sample
num_obs = 1000

obs <- matrix(NA, nrow = num_obs, ncol = length(column_names)) %>% 
  as.data.frame() %>% set_colnames(column_names)

#---- Generating IDs, female, U ----
obs[, "id"] <- seq(from = 1, to = num_obs, by = 1)
obs[, "female"] <- rbinom(num_obs, size = 1, prob = pfemale)
obs[, "U"] <- rnorm(num_obs, mean = 0, sd = 1)

#---- Generating age data ----
#Creating ages at each timepoint j
ages = matrix(seq(50, 95, by = 5), nrow = 1)
ones = matrix(1, nrow = num_obs, ncol = 1)
obs[, variable_names$age_varnames] <- ones %*% ages

#---- Generating centered age data ----
#Creating baseline-mean-centered ages at each timepoint j
obs[, variable_names$agec_varnames] <- 
  obs[, variable_names$age_varnames] - mean(age0)

#---- Generating "true" cognitive function Cij ----
#Refer to Manuscript/manuscript_equations.pdf for equation

#Generating random terms for slope and intercept
#Covariance matrices for random slope and intercept terms
cij_slope_int_cov <- lapply(1:(num_tests + 1), 
                            function(x) matrix(NA, nrow = 2, ncol = 2))
for(i in 1:(num_tests + 1)){
  cij_slope_int_cov[[i]] <- matrix(c(cij_var0, cij_cov, 
                                     cij_cov, cij_var1[i]), 
                                   nrow = 2, byrow = TRUE)
}

#Generate random terms for each individual
for(i in 1:(num_tests + 1)){
  noise <- mvrnorm(n = num_obs, mu = rep(0, 2), 
                   Sigma = cij_slope_int_cov[[i]]) 
  if(i == 1){
    obs[, c("z0_i", paste0("z1_", (i - 1), "i"))] <- noise
  } else{
    obs[, paste0("z1_", (i - 1), "i")] <- noise[, 2]
  }
}

#Generating noise term (unexplained variance in Cij) for each visit
#Creating AR(1) correlation matrix
num_visits = num_tests + 1
powers <- abs(outer(1:(num_visits), 1:(num_visits), "-")) #Exponents for autoregressive terms
corr <- sqrt(cij_var3)*(cij_r1^powers)                    #Correlation matrix
S <- diag(rep(sqrt(cij_var3)), nrow(corr))                #Diagonal matrix of SDs
cij_cov_mat <- S%*%corr%*%S                               #Covariance matrix

#Generating noise terms
obs[, variable_names$eps_varnames] <- 
  mvrnorm(n = num_obs, mu = rep(0, num_visits), Sigma = cij_cov_mat)

#Calculating Cij for each individual
#Store Cij values and slope values for each assessment
compute_Cij <- cog_func(cij_knots, cij_slopes, obs)
obs[, variable_names$Cij_varnames] <- compute_Cij[["Cij"]]
obs[, variable_names$cij_slopeij_varnames[1:num_tests]] <- 
  compute_Cij[["slopes"]]

#---- Create a competing risk outcome ----
dem_cuts_mat <- matrix(dem_cut, nrow = nrow(obs), 
                       ncol = length(variable_names$Cij_varnames), 
                       byrow = TRUE)

obs[, variable_names$dem_varnames] <- 
  (obs[, variable_names$Cij_varnames] <= dem_cuts_mat)*1

obs %<>% filter(dem0 == 0)

#---- Fit quadratic trajectories ----
sample_index <- sample(seq(1, nrow(test_data), by = 1), size = 100)

model_data <- as.data.frame(matrix(nrow = 10, ncol = 3))
colnames(model_data) <- c("Cij", "age", "age2")
model_data$age <- seq(50, 95, by = 5)
model_data$age2 <- (model_data$age)^2

for(i in 1:1){
  obs_num <- sample_index[i]
  model_data$Cij <- as.numeric(test_data[obs_num, variable_names$Cij_varnames])
  quad_model <- lm(Cij ~ age + age2, data = model_data)
}

#---- Sample plots; live people ----
plot_data <- sample_n(size = 10, test_data) %>% 
  dplyr::select(c("id", variable_names$Cij_varnames)) %>% 
  set_colnames(c("id", seq(50, 95, by = 5))) 
plot_data %<>% gather(key = "age", value = "Cij", colnames(plot_data)[-1])
plot_data$id <- as.character(plot_data$id)
plot_data$age <- as.numeric(plot_data$age)

ggplot(plot_data, aes(age, Cij, color = id)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se = FALSE) + 
  theme_minimal()

#---- Sample plots; full trajectory ----
plot_data_full <- sample_n(size = 10, obs) %>% 
  dplyr::select(c("id", variable_names$Cij_varnames)) %>% 
  set_colnames(c("id", seq(0, 45, by = 5))) 
plot_data_full %<>% 
  gather(key = "centered_age", value = "Cij", colnames(plot_data_full)[-1])
plot_data_full$id <- as.character(plot_data_full$id)
plot_data_full$centered_age <- as.numeric(plot_data_full$centered_age)

ggplot(plot_data_full, aes(centered_age, Cij, color = id)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se = FALSE) + 
  theme_minimal()

#---- Mean quadratic trajectory ----
all_obs <- obs %>% sample_n(size = 100000) %>% 
  dplyr::select(variable_names$Cij_varnames) %>% 
  set_colnames(seq(0, 45, by = 5)) %>% 
  gather(key = "centered_age", value = "Cij") %>%
  mutate_at("centered_age", as.numeric) %>% 
  mutate("centered_age2" = centered_age^2)

ggplot(all_obs, aes(centered_age, Cij)) + geom_point() + 
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se = FALSE) + 
  theme_minimal()

#---- Correlations and covariances in quadratic trajectories ----
obs <- data_gen(100000)

#Mean model
all_obs <- obs %>% 
  dplyr::select(variable_names$Cij_varnames) %>% 
  set_colnames(seq(0, 45, by = 5)) %>% 
  gather(key = "centered_age", value = "Cij") %>%
  mutate_at("centered_age", as.numeric) %>% 
  mutate("centered_age2" = centered_age^2)

fit_mean_quad <- lm(Cij ~ centered_age + centered_age2, data = all_obs)

#Individual trajectories
Cij_data <- obs %>% sample_n(5000) %>% 
  dplyr::select(variable_names$Cij_varnames) %>% 
  t() %>% as.data.frame() %>% mutate("centered_age" = seq(0, 45, by = 5), 
                                     "centered_age2" = centered_age^2)

quad_coeff <- as.data.frame(matrix(nrow = 3, ncol = 5000)) %>% 
  set_rownames(c("a0", "a1", "a2"))


for(i in 1:ncol(quad_coeff)){
  quad_coeff[, i] <- lm(Cij_data[, i] ~ centered_age + centered_age2, 
       data = Cij_data)$coefficients
}

quad_coeff %<>% t() %>% as.data.frame()

var(quad_coeff$a0)
var(quad_coeff$a1, na.rm = TRUE)
var(quad_coeff$a2, na.rm = TRUE)

cov(quad_coeff$a0, quad_coeff$a1, use = "complete.obs")
cov(quad_coeff$a0, quad_coeff$a2, use = "complete.obs")
cov(quad_coeff$a1, quad_coeff$a2, use = "complete.obs")
