################################################################################
# Lumber manufacturing example
# MBGLR chart for p = 2 and n = 5
# Calculate RL1 using whitened subgroup means
################################################################################

rm(list = ls())

# ==========================================================
# 1. IC parameters in the original scale
# ==========================================================

Sigma0_raw <- matrix(
  c(
    100, 66,
    66, 121
  ),
  nrow = 2,
  byrow = TRUE
)

mu0_raw <- c(265, 470)

# Calculate Sigma0^{-1/2}.
eig <- eigen(Sigma0_raw, symmetric = TRUE)

if (any(eig$values <= 0)) {
  stop("Sigma0_raw is not positive definite and cannot be used for whitening.")
}

Sigma0_inv_sqrt <- eig$vectors %*%
  diag(1 / sqrt(eig$values)) %*%
  t(eig$vectors)

# ==========================================================
# 2. Phase-II subgroup mean vectors in the original scale
# ==========================================================

Ybar_raw <- matrix(
  c(
    266.12, 473.17,
    271.88, 476.22,
    260.35, 462.72,
    266.88, 473.84,
    268.55, 470.86,
    266.80, 471.77,
    265.82, 465.83,
    264.05, 473.65,
    262.33, 460.25,
    270.25, 471.00,
    269.32, 479.75,
    266.00, 476.56,
    276.96, 479.36,
    263.78, 465.75,
    269.28, 476.12,
    261.78, 467.85,
    265.45, 475.85,
    265.19, 476.39,
    273.92, 475.66,
    272.23, 471.40,
    273.03, 474.40,
    266.42, 472.67,
    268.97, 472.69,
    279.58, 475.18
  ),
  nrow = 24,
  byrow = TRUE
)

colnames(Ybar_raw) <- c("Ybar1_raw", "Ybar2_raw")

# ==========================================================
# 3. Whiten the subgroup mean vectors
#    Y_i = Sigma0^{-1/2} (Ybar_i - mu0)
# ==========================================================

Y <- t(
  apply(
    Ybar_raw,
    1,
    function(ybar_i) {
      Sigma0_inv_sqrt %*% (ybar_i - mu0_raw)
    }
  )
)

colnames(Y) <- c("Y1", "Y2")

cat("\nWhitened subgroup mean vectors Y_i:\n")
print(Y)

# ==========================================================
# 4. MBGLR parameter setting in the whitened space
# ==========================================================

p <- 2
n <- 5

# Moving-window size m.
m_window <- 80

# Calibrated control limit H_MBGLR.
H_MBGLR <- 4.9980

# In the whitened space, Sigma0 is the identity matrix.
Sigma0 <- diag(p)
Sigma0_inv <- solve(Sigma0)

# IC mean vector in the whitened space.
mu0 <- rep(0, p)

# Prior mean vector m0.
m0 <- c(0.3, 0.3)

# Prior covariance matrix Z0.
Z0 <- diag(p)
Z0_inv <- solve(Z0)

# Posterior covariance matrix:
# Z1 = (Z0^{-1} + n Sigma0^{-1})^{-1}.
Z1_inv <- Z0_inv + n * Sigma0_inv
Z1 <- solve(Z1_inv)

# ==========================================================
# 5. Calculate the MBGLR statistic R_{k,m}
# ==========================================================

Kmax <- nrow(Y)

mi_mat <- matrix(0, nrow = p, ncol = Kmax)

s0 <- numeric(Kmax)

R_km_vec <- numeric(Kmax)

RL1 <- NA

for (k in 1:Kmax) {
  
  yk <- Y[k, ]
  
  # Posterior mean vector:
  # m_i = Z1 (Z0^{-1} m0 + n Sigma0^{-1} Y_i).
  
  mi_mat[, k] <- Z1 %*% Z0_inv %*% m0 +
    n * Z1 %*% Sigma0_inv %*% yk
  
  # Calculate:
  # (mu0 - m_i)' Z1^{-1} (mu0 - m_i).
  
  diff0 <- mu0 - mi_mat[, k]
  
  s0[k] <- as.numeric(
    t(diff0) %*% Z1_inv %*% diff0
  )
  
  # Candidate change-point set in the moving window.
  
  tau_start <- max(0, k - m_window)
  
  R_values <- numeric(k - tau_start)
  
  idx <- 1
  
  for (tau in tau_start:(k - 1)) {
    
    # MAP estimator of the OOC mean vector:
    # mu1_hat = sum_{i=tau+1}^k m_i / (k - tau).
    
    if (tau == k - 1) {
      
      mu1_hat <- mi_mat[, k]
      
    } else {
      
      mu1_hat <- rowMeans(
        mi_mat[, (tau + 1):k, drop = FALSE]
      )
    }
    
    # Calculate:
    # S0_part = sum_{i=tau+1}^k
    # (mu0 - m_i)' Z1^{-1} (mu0 - m_i).
    
    S0_part <- sum(s0[(tau + 1):k])
    
    # Calculate:
    # S1_part = sum_{i=tau+1}^k
    # (mu1_hat - m_i)' Z1^{-1} (mu1_hat - m_i).
    
    S1_part <- 0
    
    for (i in (tau + 1):k) {
      
      diff1 <- mu1_hat - mi_mat[, i]
      
      S1_part <- S1_part + as.numeric(
        t(diff1) %*% Z1_inv %*% diff1
      )
    }
    
    # MBGLR statistic for a given candidate tau:
    # R_tau = 0.5 * (S0_part - S1_part).
    
    R_tau <- 0.5 * (S0_part - S1_part)
    
    R_values[idx] <- R_tau
    
    idx <- idx + 1
  }
  
  # Moving-window MBGLR statistic:
  # R_{k,m} = max R_tau.
  
  R_km_vec[k] <- max(R_values)
  
  cat("k =", k, "MBGLR statistic R_km =", R_km_vec[k], "\n")
  
  # Signal rule:
  # R_{k,m} > H_MBGLR.
  
  if (is.na(RL1) && R_km_vec[k] > H_MBGLR) {
    
    RL1 <- k
    
    cat("Signal at k =", k, "\n")
  }
}

# ==========================================================
# 6. Display RL1
# ==========================================================

if (is.na(RL1)) {
  
  cat("No signal in", Kmax, "observations, so RL1 >", Kmax, "\n")
  
} else {
  
  cat("RL1 =", RL1, "\n")
}

# ==========================================================
# 7. Save and display monitoring statistics
# ==========================================================

stat_result <- data.frame(
  subgroup = 1:Kmax,
  Ybar1_raw = Ybar_raw[, 1],
  Ybar2_raw = Ybar_raw[, 2],
  Y1 = Y[, 1],
  Y2 = Y[, 2],
  posterior_m1 = mi_mat[1, ],
  posterior_m2 = mi_mat[2, ],
  MBGLR_statistic = R_km_vec,
  UCL = H_MBGLR,
  signal = R_km_vec > H_MBGLR
)

print(stat_result)

# Display the MBGLR statistic sequence.
R_km_vec


























################################################################################
# Lumber manufacturing example
# MBCUSUM and 2MBCUSUM charts for p = 2 and n = 5
# Calculate RL1 using whitened subgroup means
################################################################################

rm(list = ls())

# ==========================================================
# 1. IC parameters in the original scale
# ==========================================================

Sigma0_raw <- matrix(
  c(
    100, 66,
    66, 121
  ),
  nrow = 2,
  byrow = TRUE
)

mu0_raw <- c(265, 470)

# Calculate Sigma0^{-1/2}.
eig <- eigen(Sigma0_raw, symmetric = TRUE)

if (any(eig$values <= 0)) {
  stop("Sigma0_raw is not positive definite and cannot be used for whitening.")
}

Sigma0_inv_sqrt <- eig$vectors %*%
  diag(1 / sqrt(eig$values)) %*%
  t(eig$vectors)

cat("Check the whitening matrix. The result should be close to the identity matrix:\n")
print(Sigma0_inv_sqrt %*% Sigma0_raw %*% Sigma0_inv_sqrt)

# ==========================================================
# 2. Phase-II subgroup mean vectors in the original scale
# ==========================================================

Ybar_raw <- matrix(
  c(
    266.12, 473.17,
    271.88, 476.22,
    260.35, 462.72,
    266.88, 473.84,
    268.55, 470.86,
    266.80, 471.77,
    265.82, 465.83,
    264.05, 473.65,
    262.33, 460.25,
    270.25, 471.00,
    269.32, 479.75,
    266.00, 476.56,
    276.96, 479.36,
    263.78, 465.75,
    269.28, 476.12,
    261.78, 467.85,
    265.45, 475.85,
    265.19, 476.39,
    273.92, 475.66,
    272.23, 471.40,
    273.03, 474.40,
    266.42, 472.67,
    268.97, 472.69,
    279.58, 475.18
  ),
  nrow = 24,
  byrow = TRUE
)

colnames(Ybar_raw) <- c("Ybar1_raw", "Ybar2_raw")

Kmax <- nrow(Ybar_raw)

cat("\nNumber of Phase-II subgroup means =", Kmax, "\n")

# ==========================================================
# 3. Whiten the subgroup mean vectors
#    Y_i = Sigma0^{-1/2} (Ybar_i - mu0)
# ==========================================================

Y <- t(
  apply(
    Ybar_raw,
    1,
    function(ybar_i) {
      Sigma0_inv_sqrt %*% (ybar_i - mu0_raw)
    }
  )
)

colnames(Y) <- c("Y1", "Y2")

cat("\nWhitened subgroup mean vectors Y_i:\n")
print(Y)

# ==========================================================
# 4. Common Bayesian parameter settings in the whitened space
# ==========================================================

p <- 2
n <- 5

# In the whitened space:
# Y_i ~ MVN(mu0, Sigma0 / n) under the IC state,
# where mu0 = 0 and Sigma0 = I.

Sigma0 <- diag(p)
Sigma0_inv <- solve(Sigma0)

mu0 <- rep(0, p)

# Recommended prior mean vector for n = 5.
m0 <- c(0.3, 0.3)

# Prior covariance matrix Z0.
Z0 <- diag(p)
Z0_inv <- solve(Z0)

# Posterior covariance matrix:
# Z1 = (Z0^{-1} + n Sigma0^{-1})^{-1}.
Z1_inv <- Z0_inv + n * Sigma0_inv
Z1 <- solve(Z1_inv)

cat("\nPosterior covariance matrix Z1:\n")
print(Z1)

# ==========================================================
# 5. MBCUSUM and 2MBCUSUM parameter settings
# ==========================================================

# ----------------------------------------------------------
# MBCUSUM
# ----------------------------------------------------------
# mu1 is the reference OOC mean vector used to construct
# the MBCUSUM statistic. It is not necessarily the actual
# shift observed in the Phase-II data.

mu1 <- c(0.1, 0.1)

# Calibrated control limit for the MBCUSUM scheme.
H_MBCUSUM <- 2.6473

# ----------------------------------------------------------
# 2MBCUSUM
# ----------------------------------------------------------
# mu11 and mu12 are the two reference OOC mean vectors used
# to construct the two CUSUM components.

mu11 <- c(0.1, 0.1)
mu12 <- c(0.3, 0.3)

# Calibrated control limit for the 2MBCUSUM scheme.
H_2MBCUSUM <- 3.2998

# ==========================================================
# 6. Calculate posterior mean vectors m_i
# ==========================================================

mi_mat <- matrix(0, nrow = p, ncol = Kmax)

for (k in 1:Kmax) {
  
  yk <- Y[k, ]
  
  # Posterior mean vector:
  # m_i = Z1 (Z0^{-1} m0 + n Sigma0^{-1} Y_i).
  
  mi_mat[, k] <- Z1 %*% Z0_inv %*% m0 +
    n * Z1 %*% Sigma0_inv %*% yk
}

rownames(mi_mat) <- paste0("dim", 1:p)
colnames(mi_mat) <- paste0("k", 1:Kmax)

cat("\nPosterior mean vectors m_i:\n")
print(mi_mat)

# ==========================================================
# 7. MBCUSUM monitoring
# ==========================================================

C_i_vec <- numeric(Kmax)
W_i_vec <- numeric(Kmax)
K_vec <- numeric(Kmax)

C_prev <- 0

# For a fixed reference vector mu1, K is constant:
# K = 1/2 (mu1' Z1^{-1} mu1 - mu0' Z1^{-1} mu0).

K_const <- 0.5 * as.numeric(
  t(mu1) %*% Z1_inv %*% mu1 -
    t(mu0) %*% Z1_inv %*% mu0
)

for (k in 1:Kmax) {
  
  mk <- mi_mat[, k]
  
  # W_i = (mu1 - mu0)' Z1^{-1} m_i.
  
  W_i <- as.numeric(
    t(mu1 - mu0) %*% Z1_inv %*% mk
  )
  
  # Recursive MBCUSUM statistic:
  # C_i = max(0, C_{i-1} + W_i - K).
  
  C_i <- max(0, C_prev + W_i - K_const)
  
  W_i_vec[k] <- W_i
  K_vec[k] <- K_const
  C_i_vec[k] <- C_i
  
  C_prev <- C_i
  
  cat("k =", k, "MBCUSUM statistic C_i =", C_i, "\n")
}

# ==========================================================
# 8. 2MBCUSUM monitoring
# ==========================================================

C1_i_vec <- numeric(Kmax)
C2_i_vec <- numeric(Kmax)
Cmax_i_vec <- numeric(Kmax)

W1_i_vec <- numeric(Kmax)
W2_i_vec <- numeric(Kmax)

K1_vec <- numeric(Kmax)
K2_vec <- numeric(Kmax)

C1_prev <- 0
C2_prev <- 0

# For fixed reference vectors mu11 and mu12, K1 and K2 are constants.

K1_const <- 0.5 * as.numeric(
  t(mu11) %*% Z1_inv %*% mu11 -
    t(mu0) %*% Z1_inv %*% mu0
)

K2_const <- 0.5 * as.numeric(
  t(mu12) %*% Z1_inv %*% mu12 -
    t(mu0) %*% Z1_inv %*% mu0
)

for (k in 1:Kmax) {
  
  mk <- mi_mat[, k]
  
  # --------------------------------------------------------
  # First MBCUSUM component
  # --------------------------------------------------------
  
  # W1_i = (mu11 - mu0)' Z1^{-1} m_i.
  
  W1_i <- as.numeric(
    t(mu11 - mu0) %*% Z1_inv %*% mk
  )
  
  # C1_i = max(0, C1_{i-1} + W1_i - K1).
  
  C1_i <- max(0, C1_prev + W1_i - K1_const)
  
  # --------------------------------------------------------
  # Second MBCUSUM component
  # --------------------------------------------------------
  
  # W2_i = (mu12 - mu0)' Z1^{-1} m_i.
  
  W2_i <- as.numeric(
    t(mu12 - mu0) %*% Z1_inv %*% mk
  )
  
  # C2_i = max(0, C2_{i-1} + W2_i - K2).
  
  C2_i <- max(0, C2_prev + W2_i - K2_const)
  
  # Store the monitoring statistics.
  
  W1_i_vec[k] <- W1_i
  W2_i_vec[k] <- W2_i
  
  K1_vec[k] <- K1_const
  K2_vec[k] <- K2_const
  
  C1_i_vec[k] <- C1_i
  C2_i_vec[k] <- C2_i
  Cmax_i_vec[k] <- max(C1_i, C2_i)
  
  C1_prev <- C1_i
  C2_prev <- C2_i
  
  cat(
    "k =", k,
    "2MBCUSUM C1_i =", C1_i,
    "C2_i =", C2_i,
    "max(C1_i, C2_i) =", Cmax_i_vec[k],
    "\n"
  )
}

# ==========================================================
# 9. Calculate RL1
# ==========================================================

# ----------------------------------------------------------
# MBCUSUM RL1
# ----------------------------------------------------------
# Signal rule:
# C_i > H_MBCUSUM.

signal_MBCUSUM <- C_i_vec > H_MBCUSUM

signal_index_MBCUSUM <- which(signal_MBCUSUM)

if (length(signal_index_MBCUSUM) == 0) {
  
  RL1_MBCUSUM <- NA
  
  cat(
    "\nMBCUSUM: No signal in",
    Kmax,
    "observations, so RL1 >",
    Kmax,
    "\n"
  )
  
} else {
  
  RL1_MBCUSUM <- signal_index_MBCUSUM[1]
  
  cat("\nMBCUSUM RL1 =", RL1_MBCUSUM, "\n")
}

# ----------------------------------------------------------
# 2MBCUSUM RL1
# ----------------------------------------------------------
# Signal rule:
# max(C1_i, C2_i) > H_2MBCUSUM.

signal_2MBCUSUM <- Cmax_i_vec > H_2MBCUSUM

signal_index_2MBCUSUM <- which(signal_2MBCUSUM)

if (length(signal_index_2MBCUSUM) == 0) {
  
  RL1_2MBCUSUM <- NA
  
  cat(
    "\n2MBCUSUM: No signal in",
    Kmax,
    "observations, so RL1 >",
    Kmax,
    "\n"
  )
  
} else {
  
  RL1_2MBCUSUM <- signal_index_2MBCUSUM[1]
  
  cat("\n2MBCUSUM RL1 =", RL1_2MBCUSUM, "\n")
}

# ==========================================================
# 10. Output all monitoring statistics
# ==========================================================

stat_result <- data.frame(
  subgroup = 1:Kmax,
  
  Ybar1_raw = Ybar_raw[, 1],
  Ybar2_raw = Ybar_raw[, 2],
  
  Y1 = Y[, 1],
  Y2 = Y[, 2],
  
  posterior_m1 = mi_mat[1, ],
  posterior_m2 = mi_mat[2, ],
  
  # MBCUSUM statistics
  MBCUSUM_W_i = W_i_vec,
  MBCUSUM_K = K_vec,
  MBCUSUM_C_i = C_i_vec,
  MBCUSUM_UCL = H_MBCUSUM,
  MBCUSUM_signal = signal_MBCUSUM,
  
  # 2MBCUSUM statistics
  Two_MBCUSUM_W1_i = W1_i_vec,
  Two_MBCUSUM_K1 = K1_vec,
  Two_MBCUSUM_C1_i = C1_i_vec,
  
  Two_MBCUSUM_W2_i = W2_i_vec,
  Two_MBCUSUM_K2 = K2_vec,
  Two_MBCUSUM_C2_i = C2_i_vec,
  
  Two_MBCUSUM_max_C = Cmax_i_vec,
  Two_MBCUSUM_UCL = H_2MBCUSUM,
  Two_MBCUSUM_signal = signal_2MBCUSUM
)

cat("\nAll monitoring statistics:\n")
print(stat_result)

# ==========================================================
# 11. Save results
# ==========================================================

write.csv(
  stat_result,
  file = "lumber_MBCUSUM_2MBCUSUM_statistics.csv",
  row.names = FALSE
)

cat("\nResults saved to: lumber_MBCUSUM_2MBCUSUM_statistics.csv\n")

# ==========================================================
# 12. Summary
# ==========================================================

cat("\nSummary:\n")

cat("MBCUSUM reference vector mu1 =", paste(mu1, collapse = ", "), "\n")
cat("MBCUSUM H_MBCUSUM =", H_MBCUSUM, "\n")
cat(
  "MBCUSUM RL1 =",
  ifelse(is.na(RL1_MBCUSUM), paste0(">", Kmax), RL1_MBCUSUM),
  "\n\n"
)

cat("2MBCUSUM reference vector mu11 =", paste(mu11, collapse = ", "), "\n")
cat("2MBCUSUM reference vector mu12 =", paste(mu12, collapse = ", "), "\n")
cat("2MBCUSUM H_2MBCUSUM =", H_2MBCUSUM, "\n")
cat(
  "2MBCUSUM RL1 =",
  ifelse(is.na(RL1_2MBCUSUM), paste0(">", Kmax), RL1_2MBCUSUM),
  "\n"
)









