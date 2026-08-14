############################################################
# MBCUSUM scheme, n = 1
# Control-limit calibration by bisection method
############################################################

rm(list = ls())
gc()

library(MASS)

time1 <- proc.time()

# ==========================================================
# 1. Basic simulation settings
# ==========================================================

q_dim <- 2L                 # process dimension q
subgroup_size <- 1L         # subgroup size n
sampling_interval <- 1      # sampling interval d
tau0 <- 50L                 # change-point setting tau
num_rep <- 10000L           # number of Monte Carlo replications
K_max <- 8000L              # maximum monitoring length

# ==========================================================
# 2. IC process parameters
# ==========================================================

mu0 <- c(0, 0)

Sigma0 <- matrix(
  c(1, 0,
    0, 1),
  nrow = q_dim,
  byrow = TRUE
)

Sigma0_inv <- solve(Sigma0)

# Since Y_i is the subgroup mean vector:
# Y_i ~ MVN(mu, Sigma0 / n)
# For n = 1, Sigma_Y = Sigma0.

Sigma_Y <- Sigma0 / subgroup_size

# ==========================================================
# 3. Actual process mean for IC calibration
# ==========================================================
# For IC control-limit calibration, the process remains IC.
# Therefore, the actual post-tau mean is still mu0.

mu_actual <- mu0

# ==========================================================
# 4. Reference OOC mean vector for MBCUSUM statistic
# ==========================================================
# In the MBCUSUM statistic, mu1 is the pre-specified
# reference OOC mean vector used to construct the CUSUM chart.

mu1 <- c(0.1, 0.1)

# ==========================================================
# 5. Prior parameters
# ==========================================================

m0 <- c(0.1, 0.1)

Z0 <- matrix(
  c(1, 0,
    0, 1),
  nrow = q_dim,
  byrow = TRUE
)

Z0_inv <- solve(Z0)

# ==========================================================
# 6. Posterior covariance matrix
# ==========================================================
# Z1 = (Z0^{-1} + n Sigma0^{-1})^{-1}

Z1_inv <- Z0_inv + subgroup_size * Sigma0_inv
Z1 <- solve(Z1_inv)

# Precompute terms for m_i
# m_i = Z1 (Z0^{-1} m0 + n Sigma0^{-1} Y_i)

A_vec <- as.numeric(Z1 %*% (Z0_inv %*% m0))
B_mat <- Z1 %*% (subgroup_size * Sigma0_inv)

# Precompute terms for the MBCUSUM statistic
# W_i = (mu1 - mu0)' Z1^{-1} m_i
# K   = 1/2 (mu1' Z1^{-1} mu1 - mu0' Z1^{-1} mu0)

w_vec <- as.numeric(Z1_inv %*% (mu1 - mu0))

K_const <- as.numeric(
  0.5 * (
    t(mu1) %*% Z1_inv %*% mu1 -
      t(mu0) %*% Z1_inv %*% mu0
  )
)

# ==========================================================
# 7. Function: simulate run lengths for MBCUSUM
# ==========================================================

simulate_rls_MBCUSUM <- function(H_MBCUSUM) {
  
  rls <- rep(K_max, num_rep)
  active <- rep(TRUE, num_rep)
  
  # C_i for each Monte Carlo replication
  C_values <- rep(0, num_rep)
  
  for (k in seq_len(K_max)) {
    
    active_index <- which(active)
    n_active <- length(active_index)
    
    if (n_active == 0L) {
      break
    }
    
    # IC calibration: data are generated from the IC process.
    # Before and after tau0, the mean is still mu0.
    
    if (k <= tau0) {
      current_mu <- mu0
    } else {
      current_mu <- mu_actual
    }
    
    # Generate subgroup mean vectors Y_i
    # Y_i ~ MVN(mu, Sigma0 / n)
    
    Y_mat <- MASS::mvrnorm(
      n = n_active,
      mu = current_mu,
      Sigma = Sigma_Y
    )
    
    if (n_active == 1L) {
      Y_mat <- matrix(Y_mat, nrow = 1)
    }
    
    # Calculate posterior mean vectors m_i
    # m_i = Z1 (Z0^{-1} m0 + n Sigma0^{-1} Y_i)
    
    mi_mat <- sweep(
      Y_mat %*% t(B_mat),
      2,
      A_vec,
      "+"
    )
    
    # Calculate W_i for active replications
    # W_i = (mu1 - mu0)' Z1^{-1} m_i
    
    W_values <- as.numeric(mi_mat %*% w_vec)
    
    # Recursive MBCUSUM statistic:
    # C_i = max(0, C_{i-1} + W_i - K)
    
    C_values[active_index] <- pmax(
      0,
      C_values[active_index] + W_values - K_const
    )
    
    signal_local <- which(C_values[active_index] > H_MBCUSUM)
    
    if (length(signal_local) > 0L) {
      signal_index <- active_index[signal_local]
      rls[signal_index] <- k
      active[signal_index] <- FALSE
    }
  }
  
  return(rls)
}

# ==========================================================
# 8. Function: calculate IC SSATS for a given H_MBCUSUM
# ==========================================================

ARL_cal <- function(H_MBCUSUM) {
  
  rls <- simulate_rls_MBCUSUM(H_MBCUSUM)
  
  # Remove false alarms before tau0
  rls <- rls[rls > tau0]
  
  if (length(rls) == 0) {
    return(NA_real_)
  }
  
  # Since the change point is uniformly distributed on [tau, tau+1],
  # subtract tau + 0.5.
  # Since d = 1, multiply by d to obtain ATS scale.
  
  SSATS_values <- (rls - tau0 - 0.5) * sampling_interval
  
  SSATS <- mean(SSATS_values)
  
  return(SSATS)
}

# ==========================================================
# 9. Test one given control limit
# ==========================================================
# You can replace test_H with the value reported in your table.

#test_H <- 2.7150

#test_SSATS <- ARL_cal(test_H)

#print("Test result for MBCUSUM:")
#print(c(H_MBCUSUM = test_H, SSATS0 = test_SSATS))

# ==========================================================
# 10. Find H_MBCUSUM using the bisection method
# ==========================================================

H_low <-2.7
H_high <- 3

target_SSATS0 <- 370
tol <- 1

time_bisect1 <- proc.time()

f_low <- ARL_cal(H_low) - target_SSATS0
f_high <- ARL_cal(H_high) - target_SSATS0

if (is.na(f_low) || is.na(f_high)) {
  
  print("NA appears in the initial interval. Please increase K_max or adjust the interval.")
  
} else if (f_low * f_high > 0) {
  
  print("Wrong interval: f_low * f_high > 0")
  
  print(c(
    H_low, f_low + target_SSATS0,
    H_high, f_high + target_SSATS0
  ))
  
} else if (f_low == 0 || f_high == 0) {
  
  print("One endpoint is already the calibrated control limit")
  
  print(c(
    H_low, f_low + target_SSATS0,
    H_high, f_high + target_SSATS0
  ))
  
} else {
  
  for (iter in 1:50) {
    
    H_mid <- (H_low + H_high) / 2
    
    f_mid <- ARL_cal(H_mid) - target_SSATS0
    
    if (is.na(f_mid)) {
      
      print("NA appears during bisection. Please increase K_max or adjust the interval.")
      
      break
    }
    
    if (f_mid > -tol && f_mid < tol) {
      
      print("OK")
      
      H_MBCUSUM <- H_mid
      
      print("Calibrated control limit and estimated SSATS0:")
      print(c(H_MBCUSUM = H_MBCUSUM, SSATS0 = f_mid + target_SSATS0))
      
      print("Number of bisection iterations:")
      print(iter)
      
      break
    }
    
    if (f_low * f_mid < 0) {
      
      H_high <- H_mid
      f_high <- f_mid
      
    } else {
      
      H_low <- H_mid
      f_low <- f_mid
    }
    
    if (iter == 50) {
      
      print("The maximum number of iterations is reached")
      
      print(c(
        H_low, f_low + target_SSATS0,
        H_high, f_high + target_SSATS0
      ))
    }
  }
}

time_bisect2 <- proc.time()

print("Bisection running time:")
print(time_bisect2 - time_bisect1)

time2 <- proc.time()

print("Total running time:")
print(time2 - time1)