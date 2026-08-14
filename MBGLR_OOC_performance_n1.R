############################################################
# MBGLR scheme, n = 1, m = 1
# OOC performance evaluation
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
# 3. Prior parameters
# ==========================================================

m0 <- c(0, 0)

Z0 <- matrix(
  c(1, 0,
    0, 1),
  nrow = q_dim,
  byrow = TRUE
)

Z0_inv <- solve(Z0)

# ==========================================================
# 4. Posterior covariance matrix
# ==========================================================
# Z1 = (Z0^{-1} + n Sigma0^{-1})^{-1}

Z1_inv <- Z0_inv + subgroup_size * Sigma0_inv
Z1 <- solve(Z1_inv)

# Precompute terms for m_i
# m_i = Z1 (Z0^{-1} m0 + n Sigma0^{-1} Y_i)

A_vec <- as.numeric(Z1 %*% (Z0_inv %*% m0))
B_mat <- Z1 %*% (subgroup_size * Sigma0_inv)

# ==========================================================
# 5. Control limit
# ==========================================================
# Use the calibrated control limit from the manuscript/control-limit code.
# For n = 1, d = 1, m = 1, m0 = (0,0)':

H_MBGLR <- 2.9601

# ==========================================================
# 6. Shift settings
# ==========================================================
# The shift size is
# delta = sqrt((mu1 - mu0)' Sigma0^{-1} (mu1 - mu0)).
#
# Since Sigma0 = I, choosing Delta = c(delta, 0)
# gives shift size delta.

delta_values <- c(
  0,
  0.1, 0.2, 0.3, 0.4, 0.5,
  0.6, 0.7, 0.8, 0.9, 1.0,
  1.2, 1.4, 1.6, 1.8, 2.0,
  2.5, 3.0, 3.5, 4.0
)

# ==========================================================
# 7. Function: simulate run lengths for m = 1
# ==========================================================
# When m = 1, the candidate change point is tau = k - 1.
# Then the MAP estimate is mu1_hat = m_i, and S1 = 0.
# Therefore:
# R_{k,1} = 0.5 * (mu0 - m_i)' Z1^{-1} (mu0 - m_i)

simulate_rls_MBGLR_m1_OOC <- function(delta) {
  
  Delta <- c(delta, 0)
  mu1 <- mu0 + Delta
  
  rls <- rep(K_max, num_rep)
  active <- rep(TRUE, num_rep)
  
  for (k in seq_len(K_max)) {
    
    active_index <- which(active)
    n_active <- length(active_index)
    
    if (n_active == 0L) {
      break
    }
    
    # Before tau0, the process is IC.
    # After tau0, the process shifts to mu1.
    
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
# 8. Function: calculate OOC SSATS1 and SSSDTS
# ==========================================================

OOC_cal <- function(delta) {
  
  rls <- simulate_rls_MBGLR_m1_OOC(delta)
  
  # Remove false alarms before tau0
  rls <- rls[rls > tau0]
  
  if (length(rls) == 0) {
    return(
      data.frame(
        delta = delta,
        SSATS1 = NA_real_,
        SSSDTS = NA_real_,
        N_effective = 0
      )
    )
  }
  
  # Since the change point is uniformly distributed on [tau, tau+1],
  # subtract tau + 0.5.
  # Since d = 1, multiply by d to obtain ATS scale.
  
  SSATS_values <- (rls - tau0 - 0.5) * sampling_interval
  
  SSATS1 <- mean(SSATS_values)
  SSSDTS <- sd(SSATS_values)
  
  output <- data.frame(
    delta = delta,
    SSATS1 = SSATS1,
    SSSDTS = SSSDTS,
    N_effective = length(rls)
  )
  
  return(output)
}

# ==========================================================
# 9. Run OOC performance evaluation
# ==========================================================

results_list <- lapply(delta_values, OOC_cal)

results <- do.call(rbind, results_list)

print(results)

# ==========================================================
# 10. Save results
# ==========================================================

write.csv(
  results,
  file = "MBGLR_OOC_performance_n1_results.csv",
  row.names = FALSE
)

# ==========================================================
# 11. Running time
# ==========================================================

time2 <- proc.time()

print("Total running time:")
print(time2 - time1)