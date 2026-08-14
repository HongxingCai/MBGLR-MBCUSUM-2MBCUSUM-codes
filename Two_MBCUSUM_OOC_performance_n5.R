############################################################
# 2MBCUSUM scheme, n = 5
# OOC performance evaluation for Table 7
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
# For n = 5, Table 7 uses m0 = (0.3, 0.3)'.

m0 <- c(0.3, 0.3)

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
# 5. Actual OOC shift vectors in Table 7
# ==========================================================

actual_shift_labels <- c(
  "(0,0)",
  "(0.1,0)",
  "(0.2,0)",
  "(0.3,0)",
  "(0.4,0)",
  "(0.5,0)",
  "(0.8,0)",
  "(1.0,0)",
  "(0.6,0.6)",
  "(0.7,0.7)",
  "(0.8,0.8)",
  "(0.9,0.9)",
  "(1.0,1.0)",
  "(1.5,1.5)",
  "(2.0,2.0)"
)

actual_shift_mat <- matrix(
  c(
    0,   0,
    0.1, 0,
    0.2, 0,
    0.3, 0,
    0.4, 0,
    0.5, 0,
    0.8, 0,
    1.0, 0,
    0.6, 0.6,
    0.7, 0.7,
    0.8, 0.8,
    0.9, 0.9,
    1.0, 1.0,
    1.5, 1.5,
    2.0, 2.0
  ),
  ncol = 2,
  byrow = TRUE
)

rownames(actual_shift_mat) <- actual_shift_labels
colnames(actual_shift_mat) <- c("delta_1", "delta_2")

# ==========================================================
# 6. Reference shift-vector combinations in Table 7
# ==========================================================

reference_pair_labels <- c(
  "d11=(0.1,0.1), d12=(0.3,0.3)",
  "d11=(0.1,0.1), d12=(0.5,0.5)",
  "d11=(0.1,0.1), d12=(0.8,0.8)",
  "d11=(0.2,0.2), d12=(0.3,0.3)",
  "d11=(0.2,0.2), d12=(0.5,0.5)",
  "d11=(0.2,0.2), d12=(0.8,0.8)"
)

delta11_mat <- matrix(
  c(
    0.1, 0.1,
    0.1, 0.1,
    0.1, 0.1,
    0.2, 0.2,
    0.2, 0.2,
    0.2, 0.2
  ),
  ncol = 2,
  byrow = TRUE
)

delta12_mat <- matrix(
  c(
    0.3, 0.3,
    0.5, 0.5,
    0.8, 0.8,
    0.3, 0.3,
    0.5, 0.5,
    0.8, 0.8
  ),
  ncol = 2,
  byrow = TRUE
)

rownames(delta11_mat) <- reference_pair_labels
rownames(delta12_mat) <- reference_pair_labels

# Control limits for n = 5 in Table 7

H_values <- c(
  3.2998,
  3.3870,
  3.2309,
  3.2005,
  3.3550,
  3.2605
)

names(H_values) <- reference_pair_labels

# ==========================================================
# 7. Function: simulate run lengths for 2MBCUSUM under OOC
# ==========================================================

simulate_rls_2MBCUSUM_OOC <- function(delta_actual, delta11, delta12, H_2MBCUSUM) {
  
  # Actual OOC mean vector
  mu_actual <- mu0 + delta_actual
  
  # Reference OOC mean vectors used in the 2MBCUSUM statistic
  mu11 <- mu0 + delta11
  mu12 <- mu0 + delta12
  
  # W1_i = (mu11 - mu0)' Z1^{-1} m_i
  # K1   = 1/2 (mu11' Z1^{-1} mu11 - mu0' Z1^{-1} mu0)
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
    
    # Before tau0, the process is IC.
    # After tau0, the process shifts to mu_actual.
    
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
    
    # Calculate W1_i and W2_i
    
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
# 8. Function: calculate SSATS1 and SSSDTS
# ==========================================================

OOC_cal <- function(delta_actual, delta11, delta12, H_2MBCUSUM) {
  
  rls <- simulate_rls_2MBCUSUM_OOC(
    delta_actual = delta_actual,
    delta11 = delta11,
    delta12 = delta12,
    H_2MBCUSUM = H_2MBCUSUM
  )
  
  # Remove false alarms before tau0
  
  rls <- rls[rls > tau0]
  
  if (length(rls) == 0) {
    return(
      c(
        SSATS1 = NA_real_,
        SSSDTS = NA_real_,
        N_effective = 0
      )
    )
  }
  
  # Since the change point is uniformly distributed on [tau, tau+1],
  # subtract tau + 0.5.
  # Since d = 5, multiply by d to obtain ATS scale.
  
  SSATS_values <- (rls - tau0 - 0.5) * sampling_interval
  
  output <- c(
    SSATS1 = mean(SSATS_values),
    SSSDTS = sd(SSATS_values),
    N_effective = length(rls)
  )
  
  return(output)
}

# ==========================================================
# 9. Run all combinations
# ==========================================================

results_long <- data.frame()

SSATS_mat <- matrix(
  NA_real_,
  nrow = nrow(actual_shift_mat),
  ncol = length(reference_pair_labels)
)

rownames(SSATS_mat) <- rownames(actual_shift_mat)
colnames(SSATS_mat) <- reference_pair_labels

SSSDTS_mat <- SSATS_mat

for (ref_id in seq_along(reference_pair_labels)) {
  
  pair_label <- reference_pair_labels[ref_id]
  delta11 <- delta11_mat[ref_id, ]
  delta12 <- delta12_mat[ref_id, ]
  H_2MBCUSUM <- H_values[ref_id]
  
  for (shift_id in seq_len(nrow(actual_shift_mat))) {
    
    delta_actual <- actual_shift_mat[shift_id, ]
    shift_label <- rownames(actual_shift_mat)[shift_id]
    
    cat(
      "Running: n = 5, actual shift =", shift_label,
      ",", pair_label,
      ", H =", H_2MBCUSUM,
      "\n"
    )
    
    temp <- OOC_cal(
      delta_actual = delta_actual,
      delta11 = delta11,
      delta12 = delta12,
      H_2MBCUSUM = H_2MBCUSUM
    )
    
    SSATS_mat[shift_label, pair_label] <- temp["SSATS1"]
    SSSDTS_mat[shift_label, pair_label] <- temp["SSSDTS"]
    
    results_long <- rbind(
      results_long,
      data.frame(
        n = subgroup_size,
        d = sampling_interval,
        actual_shift = shift_label,
        actual_delta_1 = delta_actual[1],
        actual_delta_2 = delta_actual[2],
        reference_pair = pair_label,
        delta11_1 = delta11[1],
        delta11_2 = delta11[2],
        delta12_1 = delta12[1],
        delta12_2 = delta12[2],
        H_2MBCUSUM = H_2MBCUSUM,
        SSATS1 = temp["SSATS1"],
        SSSDTS = temp["SSSDTS"],
        N_effective = temp["N_effective"]
      )
    )
  }
}

# ==========================================================
# 10. Print and save results
# ==========================================================

print("SSATS1 matrix for 2MBCUSUM, n = 5:")
print(round(SSATS_mat, 2))

print("SSSDTS matrix for 2MBCUSUM, n = 5:")
print(round(SSSDTS_mat, 2))

write.csv(
  results_long,
  file = "Two_MBCUSUM_OOC_performance_n5_long_results.csv",
  row.names = FALSE
)

write.csv(
  SSATS_mat,
  file = "Two_MBCUSUM_OOC_performance_n5_SSATS1_matrix.csv"
)

write.csv(
  SSSDTS_mat,
  file = "Two_MBCUSUM_OOC_performance_n5_SSSDTS_matrix.csv"
)

# ==========================================================
# 11. Running time
# ==========================================================

time2 <- proc.time()

print("Total running time:")
print(time2 - time1)