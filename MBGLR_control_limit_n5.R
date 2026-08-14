############################################################
# MBGLR scheme, n = 5, m = 1
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
subgroup_size <- 5L         # subgroup size n
sampling_interval <- 5      # sampling interval d
window_size <- 1L           # moving-window size m
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
Sigma_Y <- Sigma0 / subgroup_size

# ==========================================================
# 3. OOC mean vector
# ==========================================================
# For IC control-limit calibration, set mu1 = mu0.

mu1 <- mu0

# Example for OOC performance:
# Delta <- c(0.1, 0)
# mu1 <- mu0 + Delta

# ==========================================================
# 4. Prior parameters
# ==========================================================
# For this setting:
# n = 5, d = 5, m = 1, m0 = (0.03, 0.03)'

m0 <- c(0.03, 0.03)

Z0 <- matrix(
  c(1, 0,
    0, 1),
  nrow = q_dim,
  byrow = TRUE
)

Z0_inv <- solve(Z0)

# ==========================================================
# 5. Posterior covariance matrix
# ==========================================================
# Z1 = (Z0^{-1} + n Sigma0^{-1})^{-1}

Z1_inv <- Z0_inv + subgroup_size * Sigma0_inv
Z1 <- solve(Z1_inv)

# Precompute terms for m_i
# m_i = Z1 (Z0^{-1} m0 + n Sigma0^{-1} Y_i)

A_vec <- as.numeric(Z1 %*% (Z0_inv %*% m0))
B_mat <- Z1 %*% (subgroup_size * Sigma0_inv)

# ==========================================================
# 6. Function: simulate run lengths for m = 1
# ==========================================================
# When m = 1, the candidate change point is tau = k - 1.
# Then the MAP estimate is mu1_hat = m_i, and S1 = 0.
# Therefore:
# R_{k,1} = 0.5 * (mu0 - m_i)' Z1^{-1} (mu0 - m_i)

simulate_rls_MBGLR_m1 <- function(H_MBGLR) {
  
  rls <- rep(K_max, num_rep)
  active <- rep(TRUE, num_rep)
  
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
      current_mu <- mu1
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
    
    # Calculate R_{k,1}
    
    diff_mat <- sweep(
      mi_mat,
      2,
      mu0,
      "-"
    )
    
    R_values <- 0.5 * rowSums(
      (diff_mat %*% Z1_inv) * diff_mat
    )
    
    signal_local <- which(R_values > H_MBGLR)
    
    if (length(signal_local) > 0L) {
      signal_index <- active_index[signal_local]
      rls[signal_index] <- k
      active[signal_index] <- FALSE
    }
  }
  
  return(rls)
}

# ==========================================================
# 7. Function: calculate IC SSATS for a given H_MBGLR
# ==========================================================

ARL_cal <- function(H_MBGLR) {
  
  rls <- simulate_rls_MBGLR_m1(H_MBGLR)
  
  # Remove false alarms before tau0
  rls <- rls[rls > tau0]
  
  if (length(rls) == 0) {
    return(NA_real_)
  }
  
  # Since the change point is uniformly distributed on [tau, tau+1],
  # subtract tau + 0.5.
  # Since d = 5, multiply by d to obtain ATS scale.
  
  SSATS_values <- (rls - tau0 - 0.5) * sampling_interval
  
  SSATS <- mean(SSATS_values)
  
  return(SSATS)
}

# ==========================================================
# 8. Test one given control limit
# ==========================================================
# For n = 5, d = 5, m = 1, m0 = (0.03, 0.03)'.
# The manuscript reports H_MBGLR = 3.5860.

#test_SSATS <- ARL_cal(3.5860)

#print("Test result for H_MBGLR = 3.5860:")
#print(test_SSATS)

# ==========================================================
# 9. Find H_MBGLR using the bisection method
# ==========================================================

H_low <- 3.5
H_high <- 3.7

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
      
      H_MBGLR <- H_mid
      
      print("Calibrated control limit and estimated SSATS0:")
      print(c(H_MBGLR = H_MBGLR, SSATS0 = f_mid + target_SSATS0))
      
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