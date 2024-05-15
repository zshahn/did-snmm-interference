#_________________________________________________________________________________________

# Notes:
# 1 time point
# A_i = (your treatment) depends on U
# Yta format
# Y0 = (baseline Y aka untreated Y) depends on U -> true outcome at t = 0
# Y10 = (Y at time 1 under no treatment) depends on U+constant -> true outcome at t = 1
# A_j = sum of your neighbors treatment (where their treatment depends on U)
# Y11 = Y10 + beta1(A_i) + beta2(A_j) *include interaction? include U?

#_________________________________________________________________________________________

# DATA GENERATION ---------------------------------------------------------

library(tidyverse)
library(nleqslv)
library(furrr) #to parallelize

seed = 6329255 #has to be 7 digits for future_map



nobs=10000
nsims=100
nboot=100

#GENERATE DATA
gen_data = function(self_effect     = 2, # (AR) try varying the relative impact of own vs. neighbor's treatments?
                    neighbor_effect = 0.5,
                    n               = nobs){
  tibble(ind    = 1:n,
         U      = rnorm(n),
         A_i    = rbinom(n, 1, plogis(U)),# (AR) inverse logit is already implemented in R as plogis
         Y0     = rnorm(n, U),
         Y10    = rnorm(n, U+1),
         lag_A  = lag(A_i, default = 0), # (AR) another option is default = rbinom(n, 1, plogis(U))
         lead_A = lead(A_i, default = 0),
         A_j    = lag_A + lead_A, #CREATE A_j: SUM OF NEIGHBORS TREATMENT
         Y11    = Y10 + self_effect*A_i + neighbor_effect*A_j,
         E_Ai   = mean(A_i), #TREATMENT EQUATIONS TO GET ESTIMATES FOR I AND J
         E_Aj   = mean(A_j))
}

# #e.g.
# data = gen_data()



# ESTIMATING EQUATIONS AND SOLVER

#ESTIMATING EQUATION:

estimator = function(data) {
  phi_A <- matrix(c(data$A_i, data$A_j), nrow = nrow(data), ncol = 2, byrow = F)
  colnames(phi_A) <- c("A_i", "A_j")

  E_phi_A <- matrix(c(data$E_Ai, data$E_Aj), nrow = nrow(data), ncol = 2, byrow = F)
  colnames(E_phi_A) <- c("E_Ai", "E_Aj")

  est_eq <- function(psi){

    H1_beta <- data$Y11 - (psi[1]*data$A_i) - (psi[2]*data$A_j)
    H0_beta <-  data$Y0
    colSums((H1_beta - H0_beta)*(phi_A - E_phi_A), na.rm = T)
  }

  #SOLVER
  ss <- nleqslv(x=c(0,0), fn=est_eq)

  tibble(psi1 = ss$x[1], psi2=ss$x[2], termcd = ss$termcd, fvec=ss$fvec[1])
}


# TRUE sampling variance ---------------------------------------------------------------
# REGENERATING THE DATA CREATION PROCESS 500 TIMES
# ______________________________________________________________________________________

true_se = . %>%
  summarise(across(psi1:psi2, sd)) %>%
  pivot_longer(psi1:psi2, names_to='parameter', values_to='true_se')
#
# ntimepoints <- 1
#
# data_inds <- lapply(1:nobs, function(i)((i-ntimepoints)*ntimepoints+1):(i*ntimepoints)) #(AR) I don't understand what's happening here - could you comment?
#
#


# BOOTSTRAP FOR SAMPLE SE -------------------------------------------------

boot_se = function(data, estimator, nboot) {
  tibble(b = 1:nboot,
         data = replicate(n = nboot,
                          expr = data[sample(1:nrow(data), size=nrow(data), replace=TRUE),],
                          simplify = FALSE),
         est = map(data, estimator)) %>%
    unnest(est) %>%
    select(-b, -data) %>%
    summarise(across(psi1:psi2, sd)) %>%
    pivot_longer(everything(), names_to='parameter', values_to='boot_se')
}

# # E.G.
# # get the true SE
# df_sims = tibble(sim = 1:nsims,
#                  data = future_map(1:nsims, ~gen_data()),
#                  est = future_map(data, estimator)) %>%
#   unnest(est)
#
#
#
# df_true_se = df_sims %>% true_se()
#
# #Get the bootstrap se in each dataset
# df_sims_boot = df_sims %>%
#   mutate(boot = future_map(data, boot_se, estimator, nboot)) %>%
#   select(sim, boot) %>%
#   unnest(boot)
#
# #then take the mean of these
# df_boot_se = df_sims_boot %>%
#   group_by(parameter) %>%
#   summarise(boot_se = mean(boot_se))
#
# #compare results
# df_true_se %>%
#   left_join(df_boot_se)



# try multiple scenarios --------------------------------------------------

df_sims_mult = expand_grid(sim = 1:nsims,
                           self_effect = c(0, 1, 2),
                           neighbor_effect = c(0, 1, 2)) %>%
  mutate(data = future_map2(self_effect, neighbor_effect, gen_data, .options = furrr_options(seed=seed)),
         est = future_map(data, estimator)) %>%
  unnest(est)


df_true_se_mult = df_sims_mult %>%
  group_by(self_effect, neighbor_effect) %>%
  true_se()


df_sims_boot_mult = df_sims_mult %>%
  mutate(boot = future_map(data, boot_se, estimator, nboot, .options = furrr_options(seed=seed))) %>%
  unnest(boot)

df_boot_se_mult = df_sims_boot_mult %>%
  group_by(self_effect, neighbor_effect, parameter) %>%
  summarise(boot_se = mean(boot_se)) %>%
  ungroup()

write_csv(
  df_true_se_mult %>%
    left_join(df_boot_se_mult) %>%
    mutate(bias = boot_se - true_se),
  file = 'results.csv')
