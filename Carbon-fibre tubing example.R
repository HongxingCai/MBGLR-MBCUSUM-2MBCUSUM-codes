############################################################
# Carbon-fibre tubing Phase-II data
# MBGLR chart for p = 3 and n = 8
# Calculate RL1 using whitened subgroup means
############################################################

rm(list = ls())

# ==========================================================
# 1. Read Phase-II data
# ==========================================================

file_path <- "D:/example/carbon2.csv"

carbon2 <- read.csv(file_path, header = TRUE)

print(head(carbon2))
print(dim(carbon2))
print(names(carbon2))

# If all 24 columns in carbon2.csv are observations, use this line.
X_wide <- as.matrix(carbon2[, 1:24])

# ==========================================================
# 2. IC parameter setting in the original scale
# ==========================================================

p <- 3
n <- 8

mu0_raw <- c(0.99, 1.04, 49.98)

Sigma0_raw <- matrix(
  c(
    0.0025, 0.0036, 0.0067,
    0.0036, 0.0145, 0.0102,
    0.0067, 0.0102, 0.0592
  ),
  nrow = p,
  byrow = TRUE
)

# ==========================================================
# 3. Check the number of data columns
# ==========================================================

if (ncol(X_wide) != p * n) {
  stop("The number of data columns is not 3*8=24. Please check whether the file contains subgroup indices or other non-data columns.")
}

Kmax <- nrow(X_wide)

cat("Number of Phase-II subgroups =", Kmax, "\n")

# ==========================================================
# 4. Recover each subgroup as an 8 by 3 matrix
#    and calculate subgroup mean vectors
# ==========================================================

Ybar_raw <- matrix(NA, nrow = Kmax, ncol = p)

for (i in 1:Kmax) {
  
  xi <- as.numeric(X_wide[i, ])
  
  Xi_mat <- matrix(
    xi,
    nrow = n,
    ncol = p,
    byrow = TRUE
  )
  
  Ybar_raw[i, ] <- colMeans(Xi_mat)
}

colnames(Ybar_raw) <- c("Ybar1_raw", "Ybar2_raw", "Ybar3_raw")

cat("\nSubgroup mean vectors in the original scale:\n")
print(Ybar_raw)

# ==========================================================
# 5. Calculate Sigma0^{-1/2} and whiten subgroup means
#    Y_i = Sigma0^{-1/2} (Ybar_i - mu0)
# ==========================================================

eig <- eigen(Sigma0_raw, symmetric = TRUE)

if (any(eig$values <= 0)) {
  stop("Sigma0_raw is not positive definite and cannot be used for whitening.")
}

Sigma0_inv_sqrt <- eig$vectors %*%
  diag(1 / sqrt(eig$values)) %*%
  t(eig$vectors)

cat("\nCheck the whitening matrix. The result should be close to the identity matrix:\n")
print(Sigma0_inv_sqrt %*% Sigma0_raw %*% Sigma0_inv_sqrt)

Y <- t(
  apply(
    Ybar_raw,
    1,
    function(ybar_i) {
      Sigma0_inv_sqrt %*% (ybar_i - mu0_raw)
    }
  )
)

colnames(Y) <- c("Y1", "Y2", "Y3")

cat("\nWhitened subgroup mean vectors Y_i:\n")
print(Y)

# ==========================================================
# 6. MBGLR parameter setting in the whitened space
# ==========================================================

# Calibrated control limit.
H_MBGLR <- 5.7123

# Moving-window size m.
m_window <- 80

# In the whitened space, Sigma0 is the identity matrix.
Sigma0 <- diag(p)
Sigma0_inv <- solve(Sigma0)

# IC mean vector in the whitened space.
mu0 <- rep(0, p)

# Prior mean vector m0.
m0 <- rep(0.3, p)

# Prior covariance matrix Z0.
Z0 <- diag(p)
Z0_inv <- solve(Z0)

# Posterior covariance matrix:
# Z1 = (Z0^{-1} + n Sigma0^{-1})^{-1}.
Z1_inv <- Z0_inv + n * Sigma0_inv
Z1 <- solve(Z1_inv)

# ==========================================================
# 7. Calculate the MBGLR statistic for all Phase-II subgroups
# ==========================================================

mi_mat <- matrix(0, nrow = p, ncol = Kmax)

s0 <- numeric(Kmax)

R_km_vec <- numeric(Kmax)

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
    # S0_part = sum_{i=tau+1}^k (mu0 - m_i)' Z1^{-1} (mu0 - m_i).
    
    S0_part <- sum(s0[(tau + 1):k])
    
    # Calculate:
    # S1_part = sum_{i=tau+1}^k (mu1_hat - m_i)' Z1^{-1} (mu1_hat - m_i).
    
    S1_part <- 0
    
    for (j in (tau + 1):k) {
      
      diff1 <- mu1_hat - mi_mat[, j]
      
      S1_part <- S1_part + as.numeric(
        t(diff1) %*% Z1_inv %*% diff1
      )
    }
    
    # MBGLR statistic for a given candidate tau:
    # R = 0.5 * (S0_part - S1_part).
    
    R_tau <- 0.5 * (S0_part - S1_part)
    
    R_values[idx] <- R_tau
    
    idx <- idx + 1
  }
  
  # Moving-window MBGLR statistic:
  # R_{k,m} = max R_tau.
  
  R_km_vec[k] <- max(R_values)
  
  cat("k =", k, "MBGLR statistic R_km =", R_km_vec[k], "\n")
}

# ==========================================================
# 8. Calculate RL1
# ==========================================================

signal_index <- which(R_km_vec > H_MBGLR)

if (length(signal_index) == 0) {
  
  RL1 <- NA
  
  cat(
    "\nNo signal in",
    Kmax,
    "Phase-II subgroups. Therefore, RL1 >",
    Kmax,
    "\n"
  )
  
} else {
  
  RL1 <- signal_index[1]
  
  cat("\nRL1 =", RL1, "\n")
}

# ==========================================================
# 9. Save all monitoring statistics
# ==========================================================

stat_result <- data.frame(
  subgroup = 1:Kmax,
  Ybar1_raw = Ybar_raw[, 1],
  Ybar2_raw = Ybar_raw[, 2],
  Ybar3_raw = Ybar_raw[, 3],
  Y1 = Y[, 1],
  Y2 = Y[, 2],
  Y3 = Y[, 3],
  MBGLR_statistic = R_km_vec,
  UCL = H_MBGLR,
  signal = R_km_vec > H_MBGLR
)

print(stat_result)

write.csv(
  stat_result,
  file = "D:/example/carbon2_MBGLR_p3_n8_all_statistics.csv",
  row.names = FALSE
)

cat(
  "\nAll MBGLR statistics have been saved to: D:/example/carbon2_MBGLR_p3_n8_all_statistics.csv\n"
)

# Display the MBGLR statistic sequence.
stat_result$MBGLR_statistic














############################################################
# Carbon-fibre tubing Phase-II data
# MBCUSUM and 2MBCUSUM charts for p = 3 and n = 8
# Calculate RL1 using whitened subgroup means
############################################################

rm(list = ls())

# ==========================================================
# 1. Read Phase-II data
# ==========================================================

file_path <- "D:/example/carbon2.csv"

carbon2 <- read.csv(file_path, header = TRUE)

print(head(carbon2))
print(dim(carbon2))
print(names(carbon2))


X_wide <- as.matrix(carbon2[, 1:24])

# ==========================================================
# 2. IC parameter setting in the original scale
# ==========================================================

p <- 3
n <- 8

mu0_raw <- c(0.99, 1.04, 49.98)

Sigma0_raw <- matrix(
  c(
    0.0025, 0.0036, 0.0067,
    0.0036, 0.0145, 0.0102,
    0.0067, 0.0102, 0.0592
  ),
  nrow = p,
  byrow = TRUE
)

# ==========================================================
# 3. Check the number of data columns
# ==========================================================

if (ncol(X_wide) != p * n) {
  stop("The number of data columns is not 3*8=24. Please check whether the file contains subgroup indices or other non-data columns.")
}

Kmax <- nrow(X_wide)

cat("Number of Phase-II subgroups =", Kmax, "\n")

# ==========================================================
# 4. Recover each subgroup as an 8 by 3 matrix
#    and calculate subgroup mean vectors
# ==========================================================

Ybar_raw <- matrix(NA, nrow = Kmax, ncol = p)

for (i in 1:Kmax) {
  
  xi <- as.numeric(X_wide[i, ])
  
  Xi_mat <- matrix(
    xi,
    nrow = n,
    ncol = p,
    byrow = TRUE
  )
  
  Ybar_raw[i, ] <- colMeans(Xi_mat)
}

colnames(Ybar_raw) <- c("Ybar1_raw", "Ybar2_raw", "Ybar3_raw")

cat("\nSubgroup mean vectors in the original scale:\n")
print(Ybar_raw)

# ==========================================================
# 5. Calculate Sigma0^{-1/2} and whiten subgroup means
#    Y_i = Sigma0^{-1/2} (Ybar_i - mu0)
# ==========================================================

eig <- eigen(Sigma0_raw, symmetric = TRUE)

if (any(eig$values <= 0)) {
  stop("Sigma0_raw is not positive definite and cannot be used for whitening.")
}

Sigma0_inv_sqrt <- eig$vectors %*%
  diag(1 / sqrt(eig$values)) %*%
  t(eig$vectors)

cat("\nCheck the whitening matrix. The result should be close to the identity matrix:\n")
print(Sigma0_inv_sqrt %*% Sigma0_raw %*% Sigma0_inv_sqrt)

Y <- t(
  apply(
    Ybar_raw,
    1,
    function(ybar_i) {
      Sigma0_inv_sqrt %*% (ybar_i - mu0_raw)
    }
  )
)

colnames(Y) <- c("Y1", "Y2", "Y3")

cat("\nWhitened subgroup mean vectors Y_i:\n")
print(Y)

# ==========================================================
# 6. MBCUSUM and 2MBCUSUM parameter settings
# ==========================================================

# Calibrated control limits.
H_MBCUSUM <- 2.35
H_2MBCUSUM <- 2.945

# In the whitened space, Sigma0 is the identity matrix.
Sigma0 <- diag(p)
Sigma0_inv <- solve(Sigma0)

# IC mean vector in the whitened space.
mu0 <- rep(0, p)

# Prior mean vector m0.
m0 <- rep(0.3, p)

# Prior covariance matrix Z0.
Z0 <- diag(p)
Z0_inv <- solve(Z0)

# Posterior covariance matrix:
# Z1 = (Z0^{-1} + n Sigma0^{-1})^{-1}.
Z1_inv <- Z0_inv + n * Sigma0_inv
Z1 <- solve(Z1_inv)

# ----------------------------------------------------------
# Reference OOC mean vector for the MBCUSUM statistic
# ----------------------------------------------------------
# This vector is not the actual shift observed in the data.
# It is the design reference vector used to construct the MBCUSUM chart.

mu1 <- c(0.1, 0.1, 0.1)

# ----------------------------------------------------------
# Reference OOC mean vectors for the 2MBCUSUM statistic
# ----------------------------------------------------------
# mu11 is usually used for detecting small shifts, while
# mu12 is used for detecting moderate or relatively large shifts.

mu11 <- c(0.1, 0.1, 0.1)
mu12 <- c(0.3, 0.3, 0.3)

# ==========================================================
# 7. Calculate posterior mean vectors m_i
# ==========================================================

mi_mat <- matrix(0, nrow = p, ncol = Kmax)

for (k in 1:Kmax) {
  
  yk <- Y[k, ]
  
  mi_mat[, k] <- Z1 %*% Z0_inv %*% m0 +
    n * Z1 %*% Sigma0_inv %*% yk
}

colnames(mi_mat) <- paste0("m", 1:Kmax)

# ==========================================================
# 8. MBCUSUM monitoring
# ==========================================================

C_i_vec <- numeric(Kmax)
W_i_vec <- numeric(Kmax)
K_vec <- numeric(Kmax)

C_prev <- 0

K_const <- 0.5 * as.numeric(
  t(mu1) %*% Z1_inv %*% mu1 -
    t(mu0) %*% Z1_inv %*% mu0
)

for (k in 1:Kmax) {
  
  mk <- mi_mat[, k]
  
  
  W_i <- as.numeric(
    t(mu1 - mu0) %*% Z1_inv %*% mk
  )
  
  C_i <- max(0, C_prev + W_i - K_const)
  
  W_i_vec[k] <- W_i
  K_vec[k] <- K_const
  C_i_vec[k] <- C_i
  
  C_prev <- C_i
  
  cat("k =", k, "MBCUSUM statistic C_i =", C_i, "\n")
}

# ==========================================================
# 9. 2MBCUSUM monitoring
# ==========================================================

C1_i_vec <- numeric(Kmax)
C2_i_vec <- numeric(Kmax)

W1_i_vec <- numeric(Kmax)
W2_i_vec <- numeric(Kmax)

K1_vec <- numeric(Kmax)
K2_vec <- numeric(Kmax)

C1_prev <- 0
C2_prev <- 0


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
  
  C1_prev <- C1_i
  C2_prev <- C2_i
  
  cat(
    "k =", k,
    "2MBCUSUM C1_i =", C1_i,
    "C2_i =", C2_i,
    "max(C1_i, C2_i) =", max(C1_i, C2_i),
    "\n"
  )
}

# ==========================================================
# 10. Calculate RL1
# ==========================================================

# Signal rule for the MBCUSUM chart:
# C_i > H_MBCUSUM.

signal_MBCUSUM <- C_i_vec > H_MBCUSUM

signal_index_MBCUSUM <- which(signal_MBCUSUM)

if (length(signal_index_MBCUSUM) == 0) {
  
  RL1_MBCUSUM <- NA
  
  cat(
    "\nMBCUSUM: No signal in",
    Kmax,
    "Phase-II subgroups. Therefore, RL1 >",
    Kmax,
    "\n"
  )
  
} else {
  
  RL1_MBCUSUM <- signal_index_MBCUSUM[1]
  
  cat("\nMBCUSUM RL1 =", RL1_MBCUSUM, "\n")
}

# Signal rule for the 2MBCUSUM chart:
# max(C1_i, C2_i) > H_2MBCUSUM.

signal_2MBCUSUM <- pmax(C1_i_vec, C2_i_vec) > H_2MBCUSUM

signal_index_2MBCUSUM <- which(signal_2MBCUSUM)

if (length(signal_index_2MBCUSUM) == 0) {
  
  RL1_2MBCUSUM <- NA
  
  cat(
    "\n2MBCUSUM: No signal in",
    Kmax,
    "Phase-II subgroups. Therefore, RL1 >",
    Kmax,
    "\n"
  )
  
} else {
  
  RL1_2MBCUSUM <- signal_index_2MBCUSUM[1]
  
  cat("\n2MBCUSUM RL1 =", RL1_2MBCUSUM, "\n")
}

# ==========================================================
# 11. Save all monitoring statistics
# ==========================================================

stat_result <- data.frame(
  subgroup = 1:Kmax,
  
  Ybar1_raw = Ybar_raw[, 1],
  Ybar2_raw = Ybar_raw[, 2],
  Ybar3_raw = Ybar_raw[, 3],
  
  Y1 = Y[, 1],
  Y2 = Y[, 2],
  Y3 = Y[, 3],
  
  posterior_m1 = mi_mat[1, ],
  posterior_m2 = mi_mat[2, ],
  posterior_m3 = mi_mat[3, ],
  
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
  
  Two_MBCUSUM_max_C = pmax(C1_i_vec, C2_i_vec),
  Two_MBCUSUM_UCL = H_2MBCUSUM,
  Two_MBCUSUM_signal = signal_2MBCUSUM
)

print(stat_result)

write.csv(
  stat_result,
  file = "D:/example/carbon2_MBCUSUM_2MBCUSUM_p3_n8_all_statistics.csv",
  row.names = FALSE
)

cat("\nAll MBCUSUM and 2MBCUSUM statistics have been saved to:\n")
cat("D:/example/carbon2_MBCUSUM_2MBCUSUM_p3_n8_all_statistics.csv\n")

# ==========================================================
# 12. Summary
# ==========================================================

cat("\nSummary:\n")

cat("MBCUSUM H_MBCUSUM =", H_MBCUSUM, "\n")
cat("MBCUSUM reference vector mu1 =", paste(mu1, collapse = ", "), "\n")
cat(
  "MBCUSUM RL1 =",
  ifelse(is.na(RL1_MBCUSUM), paste0(">", Kmax), RL1_MBCUSUM),
  "\n\n"
)

cat("2MBCUSUM H_2MBCUSUM =", H_2MBCUSUM, "\n")
cat("2MBCUSUM reference vector mu11 =", paste(mu11, collapse = ", "), "\n")
cat("2MBCUSUM reference vector mu12 =", paste(mu12, collapse = ", "), "\n")
cat(
  "2MBCUSUM RL1 =",
  ifelse(is.na(RL1_2MBCUSUM), paste0(">", Kmax), RL1_2MBCUSUM),
  "\n"
)









