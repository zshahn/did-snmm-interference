# install.packages(c("sf","tigris","spdep"))
library(sf)
library(tigris)
library(dplyr)
library(lubridate)
library(tidyr)
library(stringr)
library(spdep)

options(tigris_use_cache = TRUE)   # caches locally after first download

# Lower-res "cartographic boundary" counties
cnty <- counties(year = 2020, cb = TRUE, class = "sf")  # includes AK/HI/PR etc.

# Keep 50 states + DC only
state_keep <- c(state.abb, "DC")
cnty <- cnty[cnty$STUSPS %in% state_keep, ]

# Build rook or queen neighbors
nb_rook  <- poly2nb(cnty, queen = FALSE)  # shared borders only
nb_queen <- poly2nb(cnty, queen = TRUE)   # borders or corners

# Row-normalized weights (W-matrix)
lw_rook  <- nb2listw(nb_rook,  style = "W", zero.policy = TRUE)
W_rook   <- nb2mat(nb_rook,    style = "W", zero.policy = TRUE)

ses <- read_dta("~/snmm/medicaid/master_area_dataset_4MEPS.dta")
sesdic <- labelled::generate_dictionary(ses)

ses = ses %>% mutate(
  lag_mw = lag(st_mw_super)
)
ses$wage_hike = (ses$st_mw_super - ses$lag_mw)>=.25
data = ses[ses$year %in% 2013:2015,c("year","stcofips","stcofips_num","co_tot_pop","st_acaexp")]



library(dplyr)
library(stringr)
library(sf)
library(spdep)
library(tidyr)

pad5 <- function(x) str_pad(as.character(x), 5, pad = "0")

d = data[,c('year','st_acaexp',"stcofips","co_tot_pop")]
names(d)[3] = 'county_code'
names(d)[2] = 'expand'

# Keep state codes in cnty2
cnty2 <- cnty %>%
  mutate(county_code = stringr::str_pad(as.character(GEOID), 5, pad = "0")) %>%
  dplyr::select(county_code, STUSPS, geometry)

# Rebuild rook neighbours with this cnty2 to keep order aligned
nb_full <- spdep::poly2nb(cnty2, queen = FALSE, snap = 1e-6)

# Plain 0/1 adjacency (allows isolates)
A_bin <- spdep::nb2mat(nb_full, style = "B", zero.policy = TRUE)
diag(A_bin) <- 0

# Cross-state mask using states vector
st_vec <- cnty2$STUSPS
cross_mask <- outer(st_vec, st_vec, FUN = function(a, b) a != b)

# Cross-state adjacency
A_cross <- A_bin * cross_mask
diag(A_cross) <- 0

# Keep this order for later joins
cnty_order <- cnty2$county_code



# ---- 3) Full county-year grid + own expansion ----
years_use <- 2013:2015

panel_full <- expand_grid(county_code = cnty_order, year = years_use) %>%
  left_join(d, by = c("year","county_code"))

# If any expand is missing, set to 0 (or drop—your call)
if (anyNA(panel_full$expand)) {
  warning("Missing `expand` values set to 0.")
  panel_full <- panel_full %>% mutate(expand = if_else(is.na(expand), 0L, expand))
}

# ---- 4) Neighbor expand (1 if ANY bordering county expanded that year), using full list ----
neighbor_any <- lapply(years_use, function(yr) {
  a_vec <- panel_full %>%
    filter(year == yr) %>%
    arrange(match(county_code, cnty_order)) %>%
    pull(expand)
  tibble(
    county_code = cnty_order,
    year = yr,
    neighbor_expand = as.integer((A_bin %*% a_vec) > 0)
  )
}) %>% bind_rows()

panel_full <- panel_full %>%
  left_join(neighbor_any, by = c("county_code","year"))

# ---------- CHC (FQHC) service delivery sites → county-year counts (2013–2015) ----------

# Packages
library(httr); library(readr); library(readxl); library(dplyr); library(stringr); library(purrr)

outcome = "ACS_PCT_UNINSURED_BELOW64"
data2013 = read_excel('~/snmm/interference/SDOH_2013_COUNTY_1_0.xlsx',sheet='Data')
data2014 = read_excel('~/snmm/interference/SDOH_2014_COUNTY_1_0.xlsx',sheet='Data')
data2015 = read_excel('~/snmm/interference/SDOH_2015_COUNTY_1_0.xlsx',sheet='Data')
data2016 = read_excel('~/snmm/interference/SDOH_2016_COUNTY_1_0.xlsx',sheet='Data')

data2013 = data2013[,c("COUNTYFIPS","HRSA_MUA_COUNTY","POS_TOT_FQHC",outcome,"POS_FQHC_RATE")]
data2013$year = 2013
data2014 = data2014[,c("COUNTYFIPS","HRSA_MUA_COUNTY","POS_TOT_FQHC",outcome,"POS_FQHC_RATE")]
data2014$year = 2014
data2015 = data2015[,c("COUNTYFIPS","HRSA_MUA_COUNTY","POS_TOT_FQHC",outcome,"POS_FQHC_RATE")]
data2015$year = 2015
data2016 = data2016[,c("COUNTYFIPS","HRSA_MUA_COUNTY","POS_TOT_FQHC",outcome,"POS_FQHC_RATE")]
data2016$year = 2016
data = rbind(data2013,data2014,data2015,data2016)
names(data)[1] = "county_code"
data <- data %>%
  arrange(county_code, year)

# Quick sanity check
data %>% group_by(year) %>% summarise(nonmissing = sum(!is.na(POS_TOT_FQHC)))

# Save
write_csv(data, "~/snmm/interference/county_fqhc_counts_2013_2016.csv")


analysis_df <- panel_full %>%
  left_join(data, by = c("year","county_code")) %>%
  arrange(year, county_code) %>%
  select(year, county_code, expand, neighbor_expand, POS_TOT_FQHC,POS_FQHC_RATE,outcome)

summary(analysis_df)
head(analysis_df)


na_counties = unique(analysis_df$county_code[is.na(analysis_df$POS_FQHC_RATE)])
analysis_df=analysis_df[!analysis_df$county_code%in%na_counties, ]
vars <- c("expand","neighbor_expand",outcome)

library(dplyr); library(tidyr)

# years we’re using
years_use <- c(2013L, 2014L, 2015L)

# panel_full currently has: county_code, year, expand (cumulative 0/1)
# A_bin is the 0/1 adjacency (rows/cols ordered by cnty_order)

# 1) Make incident own exposure: D_{i,t} = expand_{i,t} - expand_{i,t-1} (with expand_{i,2012}=0)
panel_inc <- panel_full %>%
  arrange(match(county_code, cnty_order), year) %>%
  group_by(county_code) %>%
  mutate(expand_lag = dplyr::lag(expand, default = 0L),
         D = pmax(0L, expand - expand_lag)) %>%      # incident adoption (0/1)
  ungroup()

# 2) Build incident neighbor exposure for each year using the FULL graph
neighbor_inc <- lapply(years_use, function(yr) {
  d_vec <- panel_inc %>%
    filter(year == yr) %>%
    arrange(match(county_code, cnty_order)) %>%
    pull(D)

  tibble(
    county_code = cnty_order,
    year = yr,
    H = as.integer((A_cross %*% d_vec) > 0)   # 1 if ANY neighbor newly adopts in year t
  )
}) %>% bind_rows()

# 3) Attach H to panel; keep expand (cumulative) if you still need it elsewhere
panel_inc <- panel_inc %>%
  select(county_code, year, expand, D) %>%
  left_join(neighbor_inc, by = c("county_code","year"))

# 4) Merge outcomes (crude.rate) and go wide for the SNMM design
vars_keep <- c("D","H",outcome)
analysis_inc <- panel_inc %>%
  left_join(analysis_df %>% select(county_code, year, outcome),
            by = c("county_code","year"))

data2 <- analysis_inc %>%
  filter(year %in% years_use) %>%
  arrange(county_code, year) %>%
  distinct(county_code, year, .keep_all = TRUE) %>%
  tidyr::pivot_wider(
    id_cols    = county_code,
    names_from = year,
    values_from = all_of(vars_keep),
    names_glue = "{.value}_{year}"
  ) #%>%

names(data2) = c('county_code','A0','A1','A2','int0','int1','int2','Y0','Y1','Y2')


table(data2$A1)
table(data2$A2)
table(data2$int1)
table(data2$int2)

#fit same model as in simulation
build_panel_from_data2 <- function(ddd) {
  N <- nrow(ddd)
  ddd$id <- 1:N

  Hmk <- data.frame(
    id = rep(1:N, each = 5),
    m  = rep(c(1, 1, 1, 2, 2), times = N),
    k  = rep(c(0, 1, 2, 1, 2), times = N)
  )
  Hmk <- merge(Hmk, ddd[, c("id","A1","int1","A2","int2")], by = "id", all.x = TRUE)
  Hmk <- merge(Hmk, ddd[, c("id","Y0","Y1","Y2")], by = "id", all.x = TRUE)

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
  #Hmk$S9  <- Hmk$A*Hmk$int*(Hmk$m==2)
  Hmk$S10 <- Hmk$past_A*Hmk$int*(Hmk$m==2)
  Hmk$S11 <- Hmk$past_int*Hmk$A*(Hmk$m==2)
  Hmk$S12 <- Hmk$int*Hmk$past_int*(Hmk$m==2)
  #Hmk$S13 <- Hmk$int*Hmk$A*Hmk$past_int*(Hmk$m==2)

  Hmk <- Hmk[order(Hmk$id, Hmk$m, Hmk$k), ]

  # Nuisance E[phi]
  ddd$A1_hat   <- mean(ddd$A1)
  ddd$int1_hat <- ddd$A1_hat * mean(ddd$int1[ddd$A1==1]) +
    (1 - ddd$A1_hat) * mean(ddd$int1[ddd$A1==0])

  A2_mod <- glm(A2 ~ int1, data = ddd[ddd$A1==0, ], family = binomial())
  ddd$A2_hat <- 0
  if (sum(ddd$A1==0) > 0) {
    ddd$A2_hat[ddd$A1==0] <- predict(A2_mod, newdata = ddd[ddd$A1==0, ], type = "response")
  }

  int2_mod <- lm(int2 ~ A1 + int1 + A2, data = ddd)
  ddd$int2_hat <- ddd$A2_hat * predict(int2_mod, newdata = transform(ddd, A2 = 1)) +
    (1 - ddd$A2_hat) * predict(int2_mod, newdata = transform(ddd, A2 = 0))

  ddd$A1int1_hat <- mean(ddd$A1 * ddd$int1)

  ddd$A2int2     <- ddd$A2 * ddd$int2
  ddd$A2int2_hat <- 0
  idx_ok <- (ddd$A1==0 & ddd$int1!=2)
  if (sum(idx_ok) > 0) {
    A2int2_mod <- lm(A2int2 ~ int1, data = ddd[idx_ok, ])
    ddd$A2int2_hat[idx_ok] <- predict(A2int2_mod, newdata = ddd[idx_ok, ], type = "response")
  }

  Hmk <- merge(Hmk,
               rbind(data.frame(id = ddd$id, m = 1, A_hat    = ddd$A1_hat),
                     data.frame(id = ddd$id, m = 2, A_hat    = ddd$A2_hat)),
               by = c("id","m"), all.x = TRUE)
  Hmk <- merge(Hmk,
               rbind(data.frame(id = ddd$id, m = 1, int_hat  = ddd$int1_hat),
                     data.frame(id = ddd$id, m = 2, int_hat  = ddd$int2_hat)),
               by = c("id","m"), all.x = TRUE)
  Hmk <- merge(Hmk,
               rbind(data.frame(id = ddd$id, m = 1, Aint_hat = ddd$A1int1_hat),
                     data.frame(id = ddd$id, m = 2, Aint_hat = ddd$A2int2_hat)),
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
  #F9  <- ifelse((Hmk$k < 2) | (Hmk$m > 2), 0, Hmk$int2 * Hmk$A2)
  F10 <- ifelse((Hmk$k < 2) | (Hmk$m > 2), 0, Hmk$int2 * Hmk$past_A2)
  F11 <- ifelse((Hmk$k < 2) | (Hmk$m > 2), 0, Hmk$A2   * Hmk$past_int2)
  F12 <- ifelse((Hmk$k < 2) | (Hmk$m > 2), 0, Hmk$int2 * Hmk$past_int2)
  #F13 <- ifelse((Hmk$k < 2) | (Hmk$m > 2), 0, Hmk$int2 * Hmk$A2 * Hmk$past_int2)

  cbind(F1,F2,F3,F4,F5,F6,F7,F8,F10,F11,F12)
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

  # ---- 3) Instruments Z = S - E[phi] on estimating rows (unchanged) ----
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
    #H$Aint_hat*(H$m==2),
    H$int_hat*H$past_A*(H$m==2),
    H$A_hat*H$past_int*(H$m==2),
    H$int_hat*H$past_int*(H$m==2)
    #H$Aint_hat*H$past_int*(H$m==2)
  )
  Z <- as.matrix(H[, paste0("S", c(1:8,10:12))]) - ZE

  # ---- 4) Unified mask: drop rows with NA ydiff or non-finite Z/D ----
  ok <- is.finite(ydiff) &
    apply(Z, 1, function(r) all(is.finite(r))) &
    apply(D, 1, function(r) all(is.finite(r)))
  if (!any(ok)) stop("No estimating rows remain after masking.")

  Hk <- H[ok, , drop = FALSE]
  yk <- ydiff[ok]
  Dk <- D[ok, , drop = FALSE]
  Zk <- Z[ok, , drop = FALSE]

  # ---- 5) Regressor matrices per m (same as before) ----
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

estimate_one_linear <- function(d) {
  Hmk <- build_panel_from_data2(d)
  des <- prep_designs(Hmk)
  lin <- solve_linear_psi(des)
  psi_hat = lin$psi_hat
  list(psi_hat,effects_hat = c(blip(m=1,k=1,a=c(1,0),past_a=c(0,0),psi=psi_hat),
                               blip(m=1,k=2,a=c(1,0),past_a=c(0,0),psi=psi_hat),
                               blip(m=1,k=1,a=c(1,1),past_a=c(0,0),psi=psi_hat),
                               blip(m=1,k=2,a=c(1,1),past_a=c(0,0),psi=psi_hat),
                               blip(m=1,k=1,a=c(0,1),past_a=c(0,0),psi=psi_hat),
                               blip(m=1,k=2,a=c(0,1),past_a=c(0,0),psi=psi_hat),

                               blip(m=2,k=2,a=c(1,0),past_a=c(0,0),psi=psi_hat),
                               blip(m=2,k=2,a=c(1,0),past_a=c(0,1),psi=psi_hat),
                               blip(m=2,k=2,a=c(1,1),past_a=c(0,0),psi=psi_hat),
                               blip(m=2,k=2,a=c(1,1),past_a=c(0,1),psi=psi_hat),
                               blip(m=2,k=2,a=c(0,1),past_a=c(0,1),psi=psi_hat),
                               blip(m=2,k=2,a=c(0,1),past_a=c(1,0),psi=psi_hat),
                               blip(m=2,k=2,a=c(0,1),past_a=c(1,1),psi=psi_hat),
                               blip(m=2,k=2,a=c(0,1),past_a=c(0,0),psi=psi_hat)))
}


blip = function(m,k,l=0,a,past_a,psi){
  if(m==1){
    a[1]*psi[1] + a[2]*psi[2] + a[1]*(k-m)*psi[3] + a[2]*(k-m)*psi[4] + a[1]*a[2]*psi[5] +
      a[1]*a[2]*(k-m)*psi[6]
  }else{
    a[1]*psi[7] + a[2]*psi[8] + past_a[1]*a[2]*psi[9] + past_a[2]*a[1]*psi[10] +
      past_a[2]*a[2]*psi[11]
  }
}

estimate_one_linear(data2)
effects_hat = estimate_one_linear(data2)$effects_hat

library(sf)
library(dplyr)
library(purrr)
library(tibble)
library(units)


# -------------------------------------------------------------------
# Spatial block bootstrap for county (MULTIPOLYGON) data
# -------------------------------------------------------------------
# cnty2 : sf with columns `county_code` and `geometry` (MULTIPOLYGON ok)
# analysis_df : data.frame with `county_code` + variables used by stat_fn
# stat_fn     : function(plain_df) -> named numeric vector (the statistic)
# block_km    : hex block spacing in kilometers
# B           : bootstrap replicates
# sample_size_mode: "concat" (classic blocks) or "match" (resample to original N)
# target_crs  : projected CRS in meters (5070 = CONUS Albers)
#
# Returns: tibble with B rows and one column per statistic component
    county_key = "county_code"
    block_km = 75
    B = 500
    sample_size_mode = "concat"#,"match")
    target_crs = 5070
    seed = 12345

  # --- Keep only counties present in analysis; preserve sf class
  counties_use <- cnty2 %>%
    semi_join(data2 %>% distinct(!!sym(county_key)), by = county_key) %>%
    filter(!st_is_empty(geometry)) %>%
    st_transform(target_crs)

  # --- Build hex grid over study extent (in meters)
  cell_m <- set_units(block_km, "km") |> set_units("m") |> drop_units()
  hex <- st_make_grid(st_as_sfc(st_bbox(counties_use)),
                      cellsize = cell_m, square = FALSE) |>
    st_as_sf() |>
    mutate(block_id = row_number())

  # --- One representative point per county (on-surface safer than centroid)
  county_pts <- counties_use %>%
    select(!!sym(county_key)) %>%
    st_point_on_surface()

  # --- Spatial join to assign block_id; nearest fallback for edge cases
  county_pts <- st_join(county_pts, hex["block_id"], join = st_within)
  na_idx <- is.na(county_pts$block_id)
  if (any(na_idx)) {
    nn <- st_nearest_feature(county_pts[na_idx, ], hex)
    county_pts$block_id[na_idx] <- hex$block_id[nn]
  }

  # --- Map county -> block_id and attach to analysis rows
  block_map <- county_pts %>% st_drop_geometry() %>% select(!!sym(county_key), block_id)

  # Note: right_join can drop sf; we want a plain data.frame for stat_fn anyway
  dat <- data2 %>%
    left_join(block_map, by = "county_code")

  if (anyNA(dat$block_id)) {
    warning("Some rows in data2 did not get a block_id (county_code mismatch?).")
  }

  blocks <- sort(unique(dat$block_id))
  n_blocks <- length(blocks)
  N <- nrow(dat)

  # --- One bootstrap replicate: resample blocks with replacement
  one_boot <- function(){
    take <- sample(blocks, size = n_blocks, replace = TRUE)
    boot_dat <- dat[dat$block_id%in%take,]

    if (sample_size_mode == "match") {
      if (nrow(boot_dat) > N) {
        boot_dat <- boot_dat %>% slice_sample(n = N, replace = FALSE)
      } else if (nrow(boot_dat) < N) {
        boot_dat <- boot_dat %>% slice_sample(n = N, replace = TRUE)
      }
    }
    out <- estimate_one_linear(d=boot_dat)
    out
  }

boots = vector(mode='list',length=500)
for(b in 1:500){
  boots[[b]] = one_boot()
}

for(e in 1:13){
  samps = sapply(1:500,function(i)boots[[i]][[2]][e])
  print(c(effects_hat[e]-1.96*sd(samps),effects_hat[e],effects_hat[e]+1.96*sd(samps)))
}


