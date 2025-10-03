# ===========================================
# Cluster interference: SNMM, basic sandwich
# ===========================================

suppressWarnings(suppressMessages({
  library(nleqslv)
  library(numDeriv)  # Jacobian for sandwich bread
  library(MASS)      # ginv fallback
}))

inv_logit <- function(x) exp(x)/(1+exp(x))

# ----- TRUE PARAMETERS -----
psi11_t <- c(1, .5)
psi12_t <- c(2, 1)
psi22_t <- c(.75, .25, .1)
psi_true <- c(psi11_t, psi12_t, psi22_t) # length 7

# ----- BLIP FUNCTIONS (unchanged) -----
blip111 <- function(a11,a21,psi11=c(1,.5))  psi11[1]*a11 + psi11[2]*a21
blip211 <- function(a11,a21,psi11=c(1,.5))  psi11[1]*a21 + psi11[2]*a11
blip112 <- function(a11,a21,psi12=c(2,1))   psi12[1]*a11 + psi12[2]*a21
blip212 <- function(a11,a21,psi12=c(2,1))   psi12[1]*a21 + psi12[2]*a11
blip122 <- function(a12,a22,a11,a21,psi22=c(.75,.25,.1)) psi22[1]*a12 + psi22[2]*a22 + psi22[3]*a12*(a21+a22)
blip222 <- function(a12,a22,a11,a21,psi22=c(.75,.25,.1)) psi22[1]*a22 + psi22[2]*a12 + psi22[3]*a22*(a11+a12)

# ----- ONE DATASET (N clusters) -----
gen_data <- function(N) {
  U   <- rbinom(N, 1, .5)
  Y10 <- rnorm(N, U, .1)
  Y20 <- rnorm(N, U, .1)

  A11 <- rbinom(N, 1, .3 + .2*U)
  A21 <- rbinom(N, 1, .3 + .2*U)

  Y11_0 <- rnorm(N, U, .1)
  Y21_0 <- rnorm(N, U, .1)
  Y11   <- rnorm(N, Y11_0 + sapply(1:N, function(i) blip111(A11[i],A21[i],psi11_t)), .1)
  Y21   <- rnorm(N, Y21_0 + sapply(1:N, function(i) blip211(A11[i],A21[i],psi11_t)), .1)

  A12 <- ifelse(A11==1, 0, rbinom(N, 1, .3 + .2*U))
  A22 <- ifelse(A21==1, 0, rbinom(N, 1, .3 + .2*U))

  Y12_0 <- rnorm(N, U, .1)
  Y22_0 <- rnorm(N, U, .1)

  Y12 <- rnorm(
    N,
    Y12_0 +
      sapply(1:N, function(i) blip112(A11[i],A21[i],psi12_t)) +
      sapply(1:N, function(i) blip122(A12[i],A22[i],A11[i],A21[i],psi22_t)),
    .1
  )

  Y22 <- rnorm(
    N,
    Y22_0 +
      sapply(1:N, function(i) blip212(A11[i],A21[i],psi12_t)) +
      sapply(1:N, function(i) blip222(A12[i],A22[i],A11[i],A21[i],psi22_t)),
    .1
  )


  data2 <- data.frame(A11,A21,A12,A22,Y10,Y20,Y11,Y21,Y12,Y22,Y11_0,Y21_0,Y12_0,Y22_0)
  data2$id <- 1:N

  # ----- Nuisance  -----
  data2$A11_hat <- mean(data2$A11)
  data2$A21_hat <- mean(data2$A21)

  A12_mod <- glm(A12 ~ A21, data = data2[data2$A11==0, ], family = binomial())
  data2$A12_hat <- 0; data2$A12_hat[data2$A11==0] <- predict(A12_mod, newdata = data2[data2$A11==0, ], type = "response")

  A22_mod <- lm(A22 ~ A11, data = data2[data2$A21==0, ])
  data2$A22_hat <- 0; data2$A22_hat[data2$A21==0] <- predict(A22_mod, newdata = data2[data2$A21==0, ])

  data2$A12sumA2     <- data2$A12 * (data2$A21 + data2$A22)
  data2$A12sumA2_hat <- 0
  A12sumA2_mod <- lm(A12sumA2 ~ A21, data = data2[data2$A11==0, ])
  data2$A12sumA2_hat[data2$A11==0] <- predict(A12sumA2_mod, newdata = data2[data2$A11==0, ])

  data2$A22sumA1     <- data2$A22 * (data2$A11 + data2$A12)
  data2$A22sumA1_hat <- 0
  A22sumA1_mod <- lm(A22sumA1 ~ A11, data = data2[data2$A21==0, ])
  data2$A22sumA1_hat[data2$A21==0] <- predict(A22sumA1_mod, newdata = data2[data2$A21==0, ])

  # Instruments (your qi... variables)
  data2$qi1m1k1_1     <- data2$A11; data2$qi1m1k1_1_hat <- data2$A11_hat
  data2$qi1m1k1_2     <- data2$A21; data2$qi1m1k1_2_hat <- data2$A21_hat
  data2$qi2m1k1_1     <- data2$A21; data2$qi2m1k1_1_hat <- data2$A21_hat
  data2$qi2m1k1_2     <- data2$A11; data2$qi2m1k1_2_hat <- data2$A11_hat

  data2$qi1m1k2_1     <- data2$A11; data2$qi1m1k2_1_hat <- data2$A11_hat
  data2$qi1m1k2_2     <- data2$A21; data2$qi1m1k2_2_hat <- data2$A21_hat
  data2$qi2m1k2_1     <- data2$A21; data2$qi2m1k2_1_hat <- data2$A21_hat
  data2$qi2m1k2_2     <- data2$A11; data2$qi2m1k2_2_hat <- data2$A11_hat

  data2$qi1m2k2_1     <- data2$A12; data2$qi1m2k2_1_hat <- data2$A12_hat
  data2$qi1m2k2_2     <- data2$A22; data2$qi1m2k2_2_hat <- data2$A22_hat
  data2$qi1m2k2_3     <- data2$A12 * (data2$A22 + data2$A21); data2$qi1m2k2_3_hat <- data2$A12sumA2_hat

  data2$qi2m2k2_1     <- data2$A22; data2$qi2m2k2_1_hat <- data2$A22_hat
  data2$qi2m2k2_2     <- data2$A12; data2$qi2m2k2_2_hat <- data2$A12_hat
  data2$qi2m2k2_3     <- data2$A22 * (data2$A12 + data2$A11); data2$qi2m2k2_3_hat <- data2$A22sumA1_hat

  data2
}

# ----- estimating equations  -----
compute_est_eq_psi <- function(psi, estmat) {
  psi11 <- psi[1:2]; psi12 <- psi[3:4]; psi22 <- psi[5:7]

  estmat$H110 <- estmat$Y10
  estmat$H210 <- estmat$Y20

  estmat$H111 <- estmat$Y11 - (psi11[1]*estmat$A11 + psi11[2]*estmat$A21)
  estmat$H211 <- estmat$Y21 - (psi11[1]*estmat$A21 + psi11[2]*estmat$A11)

  estmat$H112 <- estmat$Y12 - (psi22[1]*estmat$A12 + psi22[2]*estmat$A22 + psi22[3]*estmat$A12*(estmat$A21+estmat$A22)) -
    (psi12[1]*estmat$A11 + psi12[2]*estmat$A21)
  estmat$H212 <- estmat$Y22 - (psi22[1]*estmat$A22 + psi22[2]*estmat$A12 + psi22[3]*estmat$A22*(estmat$A11+estmat$A12)) -
    (psi12[1]*estmat$A21 + psi12[2]*estmat$A11)

  estmat$H122 <- estmat$Y12 - (psi22[1]*estmat$A12 + psi22[2]*estmat$A22 + psi22[3]*estmat$A12*(estmat$A21+estmat$A22))
  estmat$H222 <- estmat$Y22 - (psi22[1]*estmat$A22 + psi22[2]*estmat$A12 + psi22[3]*estmat$A22*(estmat$A11+estmat$A12))

  estmat$H121 <- estmat$Y11
  estmat$H221 <- estmat$Y21

  estmat$H_diff111 <- estmat$H111 - estmat$H110
  estmat$H_diff211 <- estmat$H211 - estmat$H210
  estmat$H_diff112 <- estmat$H112 - estmat$H111
  estmat$H_diff212 <- estmat$H212 - estmat$H211
  estmat$H_diff122 <- estmat$H122 - estmat$H121
  estmat$H_diff222 <- estmat$H222 - estmat$H221

  H_mod111 <- lm(H_diff111 ~ 1,           data = estmat)
  H_mod211 <- lm(H_diff211 ~ 1,           data = estmat)
  H_mod112 <- lm(H_diff112 ~ 1,           data = estmat)
  H_mod212 <- lm(H_diff212 ~ 1,           data = estmat)
  H_mod122 <- lm(H_diff122 ~ A11*A21,     data = estmat)
  H_mod222 <- lm(H_diff222 ~ A11*A21,     data = estmat)

  estmat$V111 <- predict(H_mod111, newdata = estmat)
  estmat$V211 <- predict(H_mod211, newdata = estmat)
  estmat$V112 <- predict(H_mod112, newdata = estmat)
  estmat$V212 <- predict(H_mod212, newdata = estmat)
  estmat$V122 <- predict(H_mod122, newdata = estmat)
  estmat$V222 <- predict(H_mod222, newdata = estmat)

  eqs <- cbind(
    c( (estmat$H_diff111 - estmat$V111)*(estmat$qi1m1k1_1 - estmat$qi1m1k1_1_hat),
       (estmat$H_diff211 - estmat$V211)*(estmat$qi2m1k1_1 - estmat$qi2m1k1_1_hat) ),
    c( (estmat$H_diff111 - estmat$V111)*(estmat$qi1m1k1_2 - estmat$qi1m1k1_2_hat),
       (estmat$H_diff211 - estmat$V211)*(estmat$qi2m1k1_2 - estmat$qi2m1k1_2_hat) ),
    c( (estmat$H_diff112 - estmat$V112)*(estmat$qi1m1k2_1 - estmat$qi1m1k2_1_hat),
       (estmat$H_diff212 - estmat$V212)*(estmat$qi2m1k2_1 - estmat$qi2m1k2_1_hat) ),
    c( (estmat$H_diff112 - estmat$V112)*(estmat$qi1m1k2_2 - estmat$qi1m1k2_2_hat),
       (estmat$H_diff212 - estmat$V212)*(estmat$qi2m1k2_2 - estmat$qi2m1k2_2_hat) ),
    c( (estmat$H_diff122 - estmat$V122)*(estmat$qi1m2k2_1 - estmat$qi1m2k2_1_hat),
       (estmat$H_diff222 - estmat$V222)*(estmat$qi2m2k2_1 - estmat$qi2m2k2_1_hat) ),
    c( (estmat$H_diff122 - estmat$V122)*(estmat$qi1m2k2_2 - estmat$qi1m2k2_2_hat),
       (estmat$H_diff222 - estmat$V222)*(estmat$qi2m2k2_2 - estmat$qi2m2k2_2_hat) ),
    c( (estmat$H_diff122 - estmat$V122)*(estmat$qi1m2k2_3 - estmat$qi1m2k2_3_hat),
       (estmat$H_diff222 - estmat$V222)*(estmat$qi2m2k2_3 - estmat$qi2m2k2_3_hat) )
  )
  colSums(eqs)
}

# ----- Per-cluster (row) score vector g_i(psi) for sandwich meat -----
cluster_scores <- function(psi, estmat) {
  psi11 <- psi[1:2]; psi12 <- psi[3:4]; psi22 <- psi[5:7]

  # Build the same residual pieces (vectorized)
  H110 <- estmat$Y10
  H210 <- estmat$Y20
  H111 <- estmat$Y11 - (psi11[1]*estmat$A11 + psi11[2]*estmat$A21)
  H211 <- estmat$Y21 - (psi11[1]*estmat$A21 + psi11[2]*estmat$A11)
  H112 <- estmat$Y12 - (psi22[1]*estmat$A12 + psi22[2]*estmat$A22 + psi22[3]*estmat$A12*(estmat$A21+estmat$A22)) -
    (psi12[1]*estmat$A11 + psi12[2]*estmat$A21)
  H212 <- estmat$Y22 - (psi22[1]*estmat$A22 + psi22[2]*estmat$A12 + psi22[3]*estmat$A22*(estmat$A11+estmat$A12)) -
    (psi12[1]*estmat$A21 + psi12[2]*estmat$A11)
  H122 <- estmat$Y12 - (psi22[1]*estmat$A12 + psi22[2]*estmat$A22 + psi22[3]*estmat$A12*(estmat$A21+estmat$A22))
  H222 <- estmat$Y22 - (psi22[1]*estmat$A22 + psi22[2]*estmat$A12 + psi22[3]*estmat$A22*(estmat$A11+estmat$A12))
  H121 <- estmat$Y11
  H221 <- estmat$Y21

  H_diff111 <- H111 - H110
  H_diff211 <- H211 - H210
  H_diff112 <- H112 - H111
  H_diff212 <- H212 - H211
  H_diff122 <- H122 - H121
  H_diff222 <- H222 - H221

  # Working models for V (same as compute_est_eq_psi)
  V111 <- as.numeric(predict(lm(H_diff111 ~ 1, data = estmat)))
  V211 <- as.numeric(predict(lm(H_diff211 ~ 1, data = estmat)))
  V112 <- as.numeric(predict(lm(H_diff112 ~ 1, data = estmat)))
  V212 <- as.numeric(predict(lm(H_diff212 ~ 1, data = estmat)))
  V122 <- as.numeric(predict(lm(H_diff122 ~ A11*A21, data = estmat), newdata = estmat))
  V222 <- as.numeric(predict(lm(H_diff222 ~ A11*A21, data = estmat), newdata = estmat))

  # Row-wise 7-vector:
  g1 <- (H_diff111 - V111) * (estmat$qi1m1k1_1 - estmat$qi1m1k1_1_hat)
  g2 <- (H_diff211 - V211) * (estmat$qi2m1k1_1 - estmat$qi2m1k1_1_hat)

  g3 <- (H_diff111 - V111) * (estmat$qi1m1k1_2 - estmat$qi1m1k1_2_hat)
  g4 <- (H_diff211 - V211) * (estmat$qi2m1k1_2 - estmat$qi2m1k1_2_hat)

  g5 <- (H_diff112 - V112) * (estmat$qi1m1k2_1 - estmat$qi1m1k2_1_hat)
  g6 <- (H_diff212 - V212) * (estmat$qi2m1k2_1 - estmat$qi2m1k2_1_hat)

  g7 <- (H_diff112 - V112) * (estmat$qi1m1k2_2 - estmat$qi1m1k2_2_hat)
  g8 <- (H_diff212 - V212) * (estmat$qi2m1k2_2 - estmat$qi2m1k2_2_hat)

  g9  <- (H_diff122 - V122) * (estmat$qi1m2k2_1 - estmat$qi1m2k2_1_hat)
  g10 <- (H_diff222 - V222) * (estmat$qi2m2k2_1 - estmat$qi2m2k2_1_hat)

  g11 <- (H_diff122 - V122) * (estmat$qi1m2k2_2 - estmat$qi1m2k2_2_hat)
  g12 <- (H_diff222 - V222) * (estmat$qi2m2k2_2 - estmat$qi2m2k2_2_hat)

  g13 <- (H_diff122 - V122) * (estmat$qi1m2k2_3 - estmat$qi1m2k2_3_hat)
  g14 <- (H_diff222 - V222) * (estmat$qi2m2k2_3 - estmat$qi2m2k2_3_hat)

  # Stack by parameter groups (same columns as colSums(eqs))
  cbind(g1 + g2,   # psi11[1]
        g3 + g4,   # psi11[2]
        g5 + g6,   # psi12[1]
        g7 + g8,   # psi12[2]
        g9 + g10,  # psi22[1]
        g11 + g12, # psi22[2]
        g13 + g14  # psi22[3]
  )  # returns N x 7
}

# ----- Basic i.i.d. sandwich at cluster level -----
sandwich_se <- function(psi_hat, estmat, ridge = 1e-8) {
  N  <- nrow(estmat)
  # meat Σ̂ = (1/N) Σ g_i g_i^T  with g_i evaluated at psi_hat
  Gi <- cluster_scores(psi_hat, estmat)   # N x 7
  Sigma_hat <- crossprod(Gi) / N

  # bread G_hat = ∂ ḡ(psi) / ∂psi^T at psi_hat, where ḡ = (1/N) Σ g_i
  g_bar <- function(th) {
    compute_est_eq_psi(th, estmat) / N
  }
  G_hat <- numDeriv::jacobian(g_bar, psi_hat)  # 7 x 7

  # robust inverse with small ridge + pseudoinverse fallback
  G_reg <- G_hat + diag(ridge, 7)
  inv_try <- try(solve(G_reg), silent = TRUE)
  Ginv <- if (inherits(inv_try, "try-error")) ginv(G_hat) else inv_try

  V_hat <- (Ginv %*% Sigma_hat %*% t(Ginv)) / N
  se    <- sqrt(pmax(diag(V_hat), 0))
  list(se = se, V = V_hat, G = G_hat, Sigma = Sigma_hat)
}

# ----- One replication: solve + SE + coverage -----
one_rep <- function(N) {
  dat <- gen_data(N)
  # Solve estimating equations (sum=0)
  root <- nleqslv(x = rep(0, 7), fn = function(th) compute_est_eq_psi(th, dat))
  psi_hat <- as.numeric(root$x)

  se_out <- sandwich_se(psi_hat, dat)
  se     <- se_out$se

  ci_lo <- psi_hat - 1.96 * se
  ci_hi <- psi_hat + 1.96 * se
  cover <- (psi_true >= ci_lo) & (psi_true <= ci_hi)

  list(psi_hat = psi_hat, se = se, cover = cover)
}

# ----- Simulation wrapper -----
run_sim <- function(nsims = 1000, N = 10000, verbose = TRUE) {
  p <- length(psi_true)
  est_mat <- matrix(NA_real_, nsims, p)
  se_mat  <- matrix(NA_real_, nsims, p)
  cov_mat <- matrix(NA, nsims, p)

  for (s in 1:nsims) {
    out <- one_rep(N)
    est_mat[s, ] <- out$psi_hat
    se_mat[s, ]  <- out$se
    cov_mat[s, ] <- out$cover
    if (verbose && s %% max(1, floor(nsims/10)) == 0) {
      message(sprintf("Sim %d / %d", s, nsims))
    }
  }

  colnames(est_mat) <- colnames(se_mat) <- colnames(cov_mat) <-
    c("psi11_1","psi11_2","psi12_1","psi12_2","psi22_1","psi22_2","psi22_3")

  list(
    coverage = colMeans(cov_mat, na.rm = TRUE),
    bias     = colMeans(est_mat, na.rm = TRUE) - psi_true,
    mean_se  = colMeans(se_mat, na.rm = TRUE),
    sd_est   = apply(est_mat, 2, sd, na.rm = TRUE),
    mean_est = colMeans(est_mat, na.rm = TRUE)
  )
}

# ===== Run (example) =====
set.seed(1)
res <- run_sim(nsims = 500, N = 5000, verbose = TRUE)  # increase nsims to 1000 for paper
print(res$coverage)
print(res$bias)
print(res$mean_se)
print(res$sd_est)
