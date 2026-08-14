############################################################
# 2MBCUSUM scheme, n = 1
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
# 4. Reference OOC mean vectors for 2MBCUSUM statistic
# ==========================================================
# In the 2MBCUSUM statistic, mu11 and mu12 are the
# pre-specified reference OOC mean vectors.

mu11 <- c(0.1, 0.1)
mu12 <- c(0.3, 0.3)

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

# ==========================================================
# 7. Precompute terms for the 2MBCUSUM statistics
# ==========================================================
# W1_i = (mu11 - mu0)' Z1^{-1} m_i
# K1   = 1/2 (mu11' Z1^{-1} mu11 - mu0' Z1^{-1} mu0)
#
# W2_i = (mu12 - mu0)' Z1^{-1} m_i
# K2   = 1/2 (mu12' Z1^{-1} mu12 - mu0' Z1^{-1} mu0)

w1_vec <- as.numeric(Z1_inv %*% (mu11 - mu0))
w2_vec <- as.numeric(Z1_inv %*% (mu12 - mu0))

K1_const <- as.numeric(
  0.5 * (
    t(mu11) %*% Z1_inv %*% mu11 -
      t(mu0) %*% Z1_inv %*% mu0
  )
)

K2_const <- as.numeric(
  0.5 * (
    t(mu12) %*% Z1_inv %*% mu12 -
      t(mu0) %*% Z1_inv %*% mu0
  )
)

# ==========================================================
# 8. Function: simulate run lengths for 2MBCUSUM
# ==========================================================

simulate_rls_2MBCUSUM <- function(H_2MBCUSUM) {
  
  rls <- rep(K_max, num_rep)
  active <- rep(TRUE, num_rep)
  
  # C1_i and C2_i for each Monte Carlo replication
  C1_values <- rep(0, num_rep)
  C2_values <- rep(0, num_rep)
  
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
    
    # Calculate W1_i and W2_i for active replications
    
    W1_values <- as.numeric(mi_mat %*% w1_vec)
    W2_values <- as.numeric(mi_mat %*% w2_vec)
    
    # Recursive 2MBCUSUM statistics:
    # C1_i = max(0, C1_{i-1} + W1_i - K1)
    # C2_i = max(0, C2_{i-1} + W2_i - K2)
    
    C1_values[active_index] <- pmax(
      0,
      C1_values[active_index] + W1_values - K1_const
    )
    
    C2_values[active_index] <- pmax(
      0,
      C2_values[active_index] + W2_values - K2_const
    )
    
    # Signal rule:
    # max(C1_i, C2_i) > H_2MBCUSUM
    
    signal_local <- which(
      pmax(C1_values[active_index], C2_values[active_index]) > H_2MBCUSUM
    )
    
    if (length(signal_local) > 0L) {
      signal_index <- active_index[signal_local]
      rls[signal_index] <- k
      active[signal_index] <- FALSE
    }
  }
  
  return(rls)
}

# ==========================================================
# 9. Function: calculate IC SSATS for a given H_2MBCUSUM
# ==========================================================

ARL_cal <- function(H_2MBCUSUM) {
  
  rls <- simulate_rls_2MBCUSUM(H_2MBCUSUM)
  
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
# 10. Test one given control limit
# ==========================================================
# Your original code used H = 3.125.

#test_H <- 3.0801

#test_SSATS <- ARL_cal(test_H)

#print("Test result for 2MBCUSUM:")
#print(c(H_2MBCUSUM = test_H, SSATS0 = test_SSATS))

# ==========================================================
# 11. Find H_2MBCUSUM using the bisection method
# ==========================================================

H_low <- 2.9
H_high <- 3.1

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
      
      H_2MBCUSUM <- H_mid
      
      print("Calibrated control limit and estimated SSATS0:")
      print(c(H_2MBCUSUM = H_2MBCUSUM, SSATS0 = f_mid + target_SSATS0))
      
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