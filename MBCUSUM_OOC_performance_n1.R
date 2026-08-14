############################################################
# MBCUSUM scheme, n = 1
# OOC performance evaluation for Table 6
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

Sigma_Y <- Sigma0 / subgroup_size

# ==========================================================
# 3. Prior parameters
# ==========================================================
# For n = 1, Table 6 uses m0 = (0.1, 0.1)'.

m0 <- c(0.1, 0.1)

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
# 5. Actual OOC shift vectors in Table 6
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
# 6. Reference shift vectors and control limits in Table 6
# ==========================================================

reference_shift_labels <- c(
  "(0.1,0.1)",
  "(0.2,0.2)",
  "(0.3,0.3)",
  "(0.5,0.5)",
  "(0.8,0.8)"
)

reference_shift_mat <- matrix(
  c(
    0.1, 0.1,
    0.2, 0.2,
    0.3, 0.3,
    0.5, 0.5,
    0.8, 0.8
  ),
  ncol = 2,
  byrow = TRUE
)

rownames(reference_shift_mat) <- reference_shift_labels
colnames(reference_shift_mat) <- c("delta_ref_1", "delta_ref_2")

# Control limits for n = 1 in Table 6
H_values <- c(
  2.7150,
  2.6551,
  2.6502,
  2.6251,
  2.4869
)

names(H_values) <- reference_shift_labels

# ==========================================================
# 7. Function: simulate run lengths for MBCUSUM under OOC
# ==========================================================

simulate_rls_MBCUSUM_OOC <- function(delta_actual, delta_ref, H_MBCUSUM) {
  
  # Actual OOC mean vector
  mu_actual <- mu0 + delta_actual
  
  # Reference OOC mean vector used in the MBCUSUM statistic
  mu_ref <- mu0 + delta_ref
  
  # W_i = (mu_ref - mu0)' Z1^{-1} m_i
  # K   = 1/2 (mu_ref' Z1^{-1} mu_ref - mu0' Z1^{-1} mu0)
  
  w_vec <- as.numeric(Z1_inv %*% (mu_ref - mu0))
  
  K_const <- as.numeric(
    0.5 * (
      t(mu_ref) %*% Z1_inv %*% mu_ref -
        t(mu0) %*% Z1_inv %*% mu0
    )
  )
  
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
    
    # Calculate W_i
    
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
# 8. Function: calculate SSATS1 and SSSDTS
# ==========================================================

OOC_cal <- function(delta_actual, delta_ref, H_MBCUSUM) {
  
  rls <- simulate_rls_MBCUSUM_OOC(
    delta_actual = delta_actual,
    delta_ref = delta_ref,
    H_MBCUSUM = H_MBCUSUM
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
  # Since d = 1, multiply by d to obtain ATS scale.
  
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
  ncol = nrow(reference_shift_mat)
)

rownames(SSATS_mat) <- rownames(actual_shift_mat)
colnames(SSATS_mat) <- rownames(reference_shift_mat)

SSSDTS_mat <- SSATS_mat

for (ref_id in seq_len(nrow(reference_shift_mat))) {
  
  delta_ref <- reference_shift_mat[ref_id, ]
  ref_label <- rownames(reference_shift_mat)[ref_id]
  H_MBCUSUM <- H_values[ref_id]
  
  for (shift_id in seq_len(nrow(actual_shift_mat))) {
    
    delta_actual <- actual_shift_mat[shift_id, ]
    shift_label <- rownames(actual_shift_mat)[shift_id]
    
    cat(
      "Running: n = 1, actual shift =", shift_label,
      ", reference shift =", ref_label,
      ", H =", H_MBCUSUM,
      "\n"
    )
    
    temp <- OOC_cal(
      delta_actual = delta_actual,
      delta_ref = delta_ref,
      H_MBCUSUM = H_MBCUSUM
    )
    
    SSATS_mat[shift_label, ref_label] <- temp["SSATS1"]
    SSSDTS_mat[shift_label, ref_label] <- temp["SSSDTS"]
    
    results_long <- rbind(
      results_long,
      data.frame(
        n = subgroup_size,
        d = sampling_interval,
        actual_shift = shift_label,
        actual_delta_1 = delta_actual[1],
        actual_delta_2 = delta_actual[2],
        reference_shift = ref_label,
        reference_delta_1 = delta_ref[1],
        reference_delta_2 = delta_ref[2],
        H_MBCUSUM = H_MBCUSUM,
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

print("SSATS1 matrix for MBCUSUM, n = 1:")
print(round(SSATS_mat, 2))

print("SSSDTS matrix for MBCUSUM, n = 1:")
print(round(SSSDTS_mat, 2))

write.csv(
  results_long,
  file = "MBCUSUM_OOC_performance_n1_long_results.csv",
  row.names = FALSE
)

write.csv(
  SSATS_mat,
  file = "MBCUSUM_OOC_performance_n1_SSATS1_matrix.csv"
)

write.csv(
  SSSDTS_mat,
  file = "MBCUSUM_OOC_performance_n1_SSSDTS_matrix.csv"
)

# ==========================================================
# 11. Running time
# ==========================================================

time2 <- proc.time()

print("Total running time:")
print(time2 - time1)