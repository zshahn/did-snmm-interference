# ============================================================
# Time-varying network SNMM
#   - Moving Block Bootstrap (MBB) for SEs/CIs
#   - Monte-Carlo coverage summary
# ============================================================

suppressWarnings(suppressMessages({
  library(dplyr)
}))

set.seed(1)

# --------------------------
# True psi and "blip"
# --------------------------
blip <- function(m, k, l = 0, a, past_a, psi) {
  if (m == 1) {
    a[1]*psi[1] + a[2]*psi[2] + a[1]*(k-m)*psi[3] + a[2]*(k-m)*psi[4] +
      a[1]*a[2]*psi[5] + a[1]*a[2]*(k-m)*psi[6]
  } else {
    a[1]*psi[7] + a[2]*psi[8] + a[1]*a[2]*psi[9] +
      past_a[1]*a[2]*psi[10] + past_a[2]*a[1]*psi[11] +
      past_a[2]*a[2]*psi[12] + past_a[2]*a[1]*a[2]*psi[13]
  }
}

psi_true <- c(
  1, .5, -.1, -.1, -.2, -.05,    # m = 1 block
  1, .5, -.1, -.1, -.1, -.05, -0.05  # m = 2 block
)

# --------------------------
# Simulator (line network with neighbors i-1, i+1)
# --------------------------
simulate_one <- function(N = 5000, psi = psi_true) {
  U  <- rbinom(N, 1, .5)
  Y0 <- rnorm(N, U, .1)

  A1 <- rbinom(N, 1, .3 + .2*U)
  lag1_A1  <- c(0, A1[1:(N-1)])
  lead1_A1 <- c(A1[2:N], 0)
  int1 <- pmax(lag1_A1, lead1_A1)  # (fixed a past typo)
  phi_A1 <- cbind(A1 = A1, int1 = int1)

  Y1_0 <- rnorm(N, U, .1)
  Y1 <- rnorm(
    N,
    Y1_0 + sapply(1:N, function(i)
      blip(m = 1, k = 1, a = phi_A1[i, ], past_a = c(0, 0), psi = psi)
    ),
    .1
  )

  A2 <- ifelse(A1 == 1, 0, rbinom(N, 1, .3 + .2*U))  # absorbing
  lag1_A2  <- c(0, A2[1:(N-1)])
  lead1_A2 <- c(A2[2:N], 0)
  int2 <- pmax(lag1_A2, lead1_A2)
  phi_A2 <- cbind(A2 = A2, int2 = int2)

  Y2_0 <- rnorm(N, U, .1)
  Y2 <- rnorm(
    N,
    Y2_0 +
      sapply(1:N, function(i) blip(m = 1, k = 2, a = phi_A1[i, ], past_a = c(0, 0), psi = psi)) +
      sapply(1:N, function(i) blip(m = 2, k = 2, a = phi_A2[i, ], past_a = phi_A1[i, ], psi = psi)),
    .1
  )

  # Wide data2 we'll treat as "observed" for bootstrap
  data2 <- data.frame(A1, int1, A2, int2, Y0, Y1, Y2)
  data2$id <- 1:N
  list(data2 = data2)
}

# --------------------------
# Panel + nuisance from data2 (no recomputation of int1/int2)
# --------------------------
build_panel_from_data2 <- function(data2) {
  N <- nrow(data2)
  data2$id <- 1:N

  Hmk <- data.frame(
    id = rep(1:N, each = 5),
    m  = rep(c(1, 1, 1, 2, 2), times = N),
    k  = rep(c(0, 1, 2, 1, 2), times = N)
  )
  Hmk <- merge(Hmk, data2[, c("id","A1","int1","A2","int2")], by = "id", all.x = TRUE)
  Hmk <- merge(Hmk, data2[, c("id","Y0","Y1","Y2")], by = "id", all.x = TRUE)

  Hmk$Y <- ifelse(Hmk$k == 0, Hmk$Y0, ifelse(Hmk$k == 1, Hmk$Y1, Hmk$Y2))
  Hmk$A   <- ifelse(Hmk$m == 1, Hmk$A1,   Hmk$A2)
  Hmk$int <- ifelse(Hmk$m == 1, Hmk$int1, Hmk$int2)

  Hmk$past_A    <- ifelse(Hmk$m == 1, 0, Hmk$A1)
  Hmk$past_A1   <- 0
  Hmk$past_A2   <- Hmk$A1
  Hmk$past_int  <- ifelse(Hmk$m == 1, 0, Hmk$int1)
  Hmk$past_int1 <- 0
  Hmk$past_int2 <- Hmk$int1

  # S1..S13 (basis)
  Hmk$S1  <- Hmk$A*(Hmk$m==1)
  Hmk$S2  <- Hmk$int*(Hmk$m==1)
  Hmk$S3  <- Hmk$A*(Hmk$k-Hmk$m)*(Hmk$m==1)
  Hmk$S4  <- Hmk$int*(Hmk$k-Hmk$m)
  Hmk$S5  <- Hmk$A*Hmk$int*(Hmk$m==1)
  Hmk$S6  <- Hmk$A*Hmk$int*(Hmk$k-Hmk$m)
  Hmk$S7  <- Hmk$A*(Hmk$m==2)
  Hmk$S8  <- Hmk$int*(Hmk$m==2)
  Hmk$S9  <- Hmk$A*Hmk$int*(Hmk$m==2)
  Hmk$S10 <- Hmk$past_A*Hmk$int*(Hmk$m==2)
  Hmk$S11 <- Hmk$past_int*Hmk$A*(Hmk$m==2)
  Hmk$S12 <- Hmk$int*Hmk$past_int*(Hmk$m==2)
  Hmk$S13 <- Hmk$int*Hmk$A*Hmk$past_int*(Hmk$m==2)

  Hmk <- Hmk[order(Hmk$id, Hmk$m, Hmk$k), ]

  # Nuisance E[phi]
  data2$A1_hat   <- mean(data2$A1)
  data2$int1_hat <- data2$A1_hat * mean(data2$int1[data2$A1==1]) +
    (1 - data2$A1_hat) * mean(data2$int1[data2$A1==0])

  A2_mod <- glm(A2 ~ int1, data = data2[data2$A1==0, ], family = binomial())
  data2$A2_hat <- 0
  if (sum(data2$A1==0) > 0) {
    data2$A2_hat[data2$A1==0] <- predict(A2_mod, newdata = data2[data2$A1==0, ], type = "response")
  }

  int2_mod <- lm(int2 ~ A1 + int1 + A2, data = data2)
  data2$int2_hat <- data2$A2_hat * predict(int2_mod, newdata = transform(data2, A2 = 1)) +
    (1 - data2$A2_hat) * predict(int2_mod, newdata = transform(data2, A2 = 0))

  data2$A1int1_hat <- mean(data2$A1 * data2$int1)

  data2$A2int2     <- data2$A2 * data2$int2
  data2$A2int2_hat <- 0
  idx_ok <- (data2$A1==0 & data2$int1!=2)
  if (sum(idx_ok) > 0) {
    A2int2_mod <- lm(A2int2 ~ int1, data = data2[idx_ok, ])
    data2$A2int2_hat[idx_ok] <- predict(A2int2_mod, newdata = data2[idx_ok, ], type = "response")
  }

  Hmk <- merge(Hmk,
               rbind(data.frame(id = data2$id, m = 1, A_hat    = data2$A1_hat),
                     data.frame(id = data2$id, m = 2, A_hat    = data2$A2_hat)),
               by = c("id","m"), all.x = TRUE)
  Hmk <- merge(Hmk,
               rbind(data.frame(id = data2$id, m = 1, int_hat  = data2$int1_hat),
                     data.frame(id = data2$id, m = 2, int_hat  = data2$int2_hat)),
               by = c("id","m"), all.x = TRUE)
  Hmk <- merge(Hmk,
               rbind(data.frame(id = data2$id, m = 1, Aint_hat = data2$A1int1_hat),
                     data.frame(id = data2$id, m = 2, Aint_hat = data2$A2int2_hat)),
               by = c("id","m"), all.x = TRUE)

  Hmk <- Hmk[order(Hmk$id, Hmk$m, Hmk$k), ]
  Hmk
}

# --------------------------
# Gamma basis and lag-difference helpers (13 cols)
# --------------------------
gamma_basis_matrix <- function(Hmk) {
  F1 <- ifelse((Hmk$k < 1) | (Hmk$m > 1), 0, Hmk$A1)
  F2 <- ifelse((Hmk$k < 1) | (Hmk$m > 1), 0, Hmk$int1)
  F3 <- ifelse((Hmk$k < 1) | (Hmk$m > 1), 0, Hmk$A1  * (Hmk$k - 1))
  F4 <- ifelse((Hmk$k < 1) | (Hmk$m > 1), 0, Hmk$int1* (Hmk$k - 1))
  F5 <- ifelse((Hmk$k < 1) | (Hmk$m > 1), 0, Hmk$int1* Hmk$A1)
  F6 <- ifelse((Hmk$k < 1) | (Hmk$m > 1), 0, Hmk$int1* Hmk$A1 * (Hmk$k - 1))

  F7  <- ifelse((Hmk$k < 2) | (Hmk$m > 2), 0, Hmk$A2)
  F8  <- ifelse((Hmk$k < 2) | (Hmk$m > 2), 0, Hmk$int2)
  F9  <- ifelse((Hmk$k < 2) | (Hmk$m > 2), 0, Hmk$int2 * Hmk$A2)
  F10 <- ifelse((Hmk$k < 2) | (Hmk$m > 2), 0, Hmk$int2 * Hmk$past_A2)
  F11 <- ifelse((Hmk$k < 2) | (Hmk$m > 2), 0, Hmk$A2   * Hmk$past_int2)
  F12 <- ifelse((Hmk$k < 2) | (Hmk$m > 2), 0, Hmk$int2 * Hmk$past_int2)
  F13 <- ifelse((Hmk$k < 2) | (Hmk$m > 2), 0, Hmk$int2 * Hmk$A2 * Hmk$past_int2)

  cbind(F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12,F13)
}

gamma_basis_lagdiff <- function(Hmk, basis_mat) {
  ord <- order(Hmk$id, Hmk$m, Hmk$k)
  Hmk_ord <- Hmk[ord, ]
  B <- basis_mat[ord, , drop = FALSE]
  start <- c(TRUE, diff(Hmk_ord$id) != 0 | diff(Hmk_ord$m) != 0)
  Blag <- B
  Blag[start, ] <- 0
  if (nrow(B) > 1) {
    Blag[!start, ] <- B[which(!start) - 1L, , drop = FALSE]
  }
  D <- B - Blag
  D[order(ord), , drop = FALSE]
}

# --------------------------
# Linear estimator with row masking
# --------------------------
prep_designs <- function(Hmk) {
  # Rows that enter the moment equations (unchanged)
  est_rows <- which(Hmk$k >= Hmk$m & Hmk$m > 0)

  # ---- 1) Compute ydiff on the FULL panel (k = 0,1,2 inside each (id,m)) ----
  ord_full <- order(Hmk$id, Hmk$m, Hmk$k)
  H_full   <- Hmk[ord_full, ]

  # lag within (id,m)
  start <- c(TRUE, diff(H_full$id) != 0 | diff(H_full$m) != 0)
  y     <- H_full$Y
  ylag  <- y
  ylag[start] <- NA_real_
  if (length(y) > 1) {
    ylag[!start] <- y[which(!start) - 1L]
  }
  ydiff_full <- y - ylag
  # back to original row order of Hmk, then subset to estimating rows
  ydiff_full <- ydiff_full[order(ord_full)]
  ydiff <- ydiff_full[est_rows]

  # ---- 2) Basis difference D on estimating rows (unchanged) ----
  basis_all <- gamma_basis_matrix(Hmk)
  D_all     <- gamma_basis_lagdiff(Hmk, basis_all)
  D <- D_all[est_rows, , drop = FALSE]

  # ---- 3)Z = S - E[phi] on estimating rows (unchanged) ----
  H <- Hmk[est_rows, ]  # shorthand

  ZE <- cbind(
    H$A_hat*(H$m==1),
    H$int_hat*(H$m==1),
    H$A_hat*(H$k-H$m)*(H$m==1),
    H$int_hat*(H$k-H$m)*(H$m==1),
    H$Aint_hat*(H$m==1),
    H$Aint_hat*(H$k-H$m)*(H$m==1),

    H$A_hat*(H$m==2),
    H$int_hat*(H$m==2),
    H$Aint_hat*(H$m==2),
    H$int_hat*H$past_A*(H$m==2),
    H$A_hat*H$past_int*(H$m==2),
    H$int_hat*H$past_int*(H$m==2),
    H$Aint_hat*H$past_int*(H$m==2)
  )
  Z <- as.matrix(H[, paste0("S", 1:13)]) - ZE

  # ---- 4) Unified mask: drop rows with NA ydiff or non-finite Z/D ----
  ok <- is.finite(ydiff) &
    apply(Z, 1, function(r) all(is.finite(r))) &
    apply(D, 1, function(r) all(is.finite(r)))
  if (!any(ok)) stop("No estimating rows remain after masking.")

  Hk <- H[ok, , drop = FALSE]
  yk <- ydiff[ok]
  Dk <- D[ok, , drop = FALSE]
  Zk <- Z[ok, , drop = FALSE]

  # ---- 5) Regressor matrices per m  ----
  I1 <- which(Hk$m == 1)
  I2 <- which(Hk$m == 2)
  X1 <- cbind(1, Hk$k[I1])                   # intercept + k   (now k has variation)
  X2 <- model.matrix(~ A1*int1, data = Hk[I2, ])

  list(H = Hk, ydiff = yk, D = Dk, Z = Zk, rows_m1 = I1, rows_m2 = I2,
       X1 = X1, X2 = X2)
}


solve_linear_psi <- function(des) {
  with(des, {
    # (I-P) y by m
    y1 <- ydiff[rows_m1]; XtX1 <- crossprod(X1); Xty1 <- crossprod(X1, y1)
    beta1 <- solve(XtX1, Xty1); r1 <- as.numeric(y1 - X1 %*% beta1)

    y2 <- ydiff[rows_m2]; XtX2 <- crossprod(X2); Xty2 <- crossprod(X2, y2)
    beta2 <- solve(XtX2, Xty2); r2 <- as.numeric(y2 - X2 %*% beta2)

    r <- numeric(nrow(H)); r[rows_m1] <- r1; r[rows_m2] <- r2  # (I-P) y

    # (I-P) D by m
    Dt1 <- D[rows_m1, , drop = FALSE]; XtD1 <- crossprod(X1, Dt1)
    R1  <- Dt1 - X1 %*% solve(XtX1, XtD1)

    Dt2 <- D[rows_m2, , drop = FALSE]; XtD2 <- crossprod(X2, Dt2)
    R2  <- Dt2 - X2 %*% solve(XtX2, XtD2)

    R <- matrix(0, nrow = nrow(D), ncol = ncol(D))
    R[rows_m1, ] <- R1; R[rows_m2, ] <- R2

    # linear solve: (Z'R) psi = Z'r
    M <- crossprod(Z, R)   # 13x13
    b <- crossprod(Z, r)   # 13x1
    psi_hat <- as.numeric(solve(M, b))
    list(psi_hat = psi_hat)
  })
}

estimate_one_linear <- function(data2) {
  Hmk <- build_panel_from_data2(data2)
  des <- prep_designs(Hmk)
  lin <- solve_linear_psi(des)
  lin$psi_hat
}

# --------------------------
# Moving Block Bootstrap over IDs (circular)
# --------------------------
mbb_resample_data2 <- function(data2, L) {
  N <- nrow(data2)
  Bk <- ceiling(N / L)
  starts <- sample.int(N, size = Bk, replace = TRUE)
  idx <- unlist(lapply(starts, function(s) ((s-1) + 0:(L-1)) %% N + 1))
  idx <- idx[1:N]
  boot <- data2[idx, , drop = FALSE]
  rownames(boot) <- NULL
  boot$id <- 1:N
  boot
}

solve_on_data2 <- function(data2) {
  # wrap in try to allow failures without killing the bootstrap
  out <- try(estimate_one_linear(data2), silent = TRUE)
  if (inherits(out, "try-error")) return(rep(NA_real_, 13))
  out
}

mbb_boot <- function(data2_obs, B = 200, L = 5, ci_level = 0.95, verbose = TRUE) {
  alpha <- 1 - ci_level
  p <- 13

  psi0 <- estimate_one_linear(data2_obs)

  PSI <- matrix(NA_real_, B, p)
  for (b in 1:B) {
    d2b <- mbb_resample_data2(data2_obs, L = L)
    PSI[b, ] <- solve_on_data2(d2b)
    if (verbose && (b %% max(1, floor(B/10)) == 0)) {
      message(sprintf("Boot %d / %d", b, B))
    }
  }

  ok <- apply(PSI, 1, function(r) all(is.finite(r)))
  PSIok <- PSI[ok, , drop = FALSE]
  if (nrow(PSIok) == 0) {
    warning("Only 0 successful bootstrap replicates; CIs may be unstable.")
    return(list(psi_hat = psi0,
                boot_se = rep(NA_real_, p),
                ci_percentile = cbind(lo = rep(NA_real_, p), hi = rep(NA_real_, p)),
                ci_basic = cbind(lo = rep(NA_real_, p), hi = rep(NA_real_, p)),
                boot_rep = PSI,
                converged = ok))
  }

  boot_se <- apply(PSIok, 2, sd)

  q_lo <- apply(PSIok, 2, quantile, probs = alpha/2, na.rm = TRUE)
  q_hi <- apply(PSIok, 2, quantile, probs = 1 - alpha/2, na.rm = TRUE)
  ci_pct <- cbind(lo = q_lo, hi = q_hi)

  diffs <- sweep(PSIok, 2, psi0, "-")
  d_lo <- apply(diffs, 2, quantile, probs = 1 - alpha/2, na.rm = TRUE)
  d_hi <- apply(diffs, 2, quantile, probs = alpha/2,     na.rm = TRUE)
  ci_basic <- cbind(lo = psi0 - d_lo, hi = psi0 - d_hi)

  list(
    psi_hat = psi0,
    boot_se = boot_se,
    ci_percentile = ci_pct,
    ci_basic = ci_basic,
    boot_rep = PSI,
    converged = ok
  )
}

# --------------------------
# Monte-Carlo runner
# --------------------------
run_sim <- function(nsims = 100, N = 5000, B = 200, L = 5, ci_level = 0.95, verbose = TRUE) {
  p <- 13
  est_mat <- matrix(NA_real_, nsims, p)
  se_mat  <- matrix(NA_real_, nsims, p)
  cover_pct <- matrix(NA, nsims, p)
  cover_basic <- matrix(NA, nsims, p)

  for (s in 1:nsims) {
    sim  <- simulate_one(N = N)
    d2   <- sim$data2

    # baseline estimate
    psi_hat <- estimate_one_linear(d2)
    est_mat[s, ] <- psi_hat

    # block bootstrap (percentile + basic)
    mbb <- mbb_boot(d2, B = B, L = L, ci_level = ci_level, verbose = FALSE)
    se_mat[s, ] <- mbb$boot_se

    # coverage
    ci_pct   <- mbb$ci_percentile
    ci_basic <- mbb$ci_basic
    cover_pct[s, ]   <- (psi_true >= ci_pct[,1]) & (psi_true <= ci_pct[,2])
    cover_basic[s, ] <- (psi_true >= ci_basic[,1]) & (psi_true <= ci_basic[,2])

    if (verbose && s %% max(1, floor(nsims/10)) == 0) {
      message(sprintf("Sim %d / %d", s, nsims))
    }
  }

  colnames(est_mat) <- colnames(se_mat) <- paste0("psi", 1:p)

  list(
    coverage_percentile = colMeans(cover_pct,   na.rm = TRUE),
    coverage_basic      = colMeans(cover_basic, na.rm = TRUE),
    bias                = colMeans(est_mat, na.rm = TRUE) - psi_true,
    mean_boot_se        = colMeans(se_mat, na.rm = TRUE),
    sd_est              = apply(est_mat, 2, sd, na.rm = TRUE),
    mean_est            = colMeans(est_mat, na.rm = TRUE)
  )
}


nsims <- 500
N     <- 5000
B     <- 200
L     <- 5

res <- run_sim(nsims = nsims, N = N, B = B, L = L, ci_level = 0.95, verbose = TRUE)
print(res$coverage_percentile)
print(res$coverage_basic)
print(res$bias)
print(res$mean_boot_se)
print(res$sd_est)
