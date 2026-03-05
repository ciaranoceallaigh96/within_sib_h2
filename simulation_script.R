#Simulation to verify new within-sibling R² formulas
#Uses Gaussian genetic values 
#Code inspired by Lee et al. -> "A better coefficient of determination for genetic profile analysis", Genetic Epidemiology (2012)
#Uses new within-sib metric on the liability scale and comparable to within-sibling heritability
#oracle PRS = true additive genetic value = perfect/optimal additive PRS predictor

#Note:
#For real-world discordant sampling, P(D) and t_tet should be estimated from an unascertained sibling sample (not discordant pairs only).
#Not modelling relatedness between families, assortative mating, inbreeding, GxG, GxE, indirect effects etc.

library(MASS)
library(fMultivar)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(jsonlite)

R <- 1000   #replicates
K_grid_all   <- c(0.01, 0.05, 0.10, 0.20, 0.50)      #for P(D) vs K
K_grid_sweep <- c(0.01, 0.10, 0.50)                  
h2_grid      <- c(0.10, 0.50, 0.90)                  

d2_grid      <- c(0.0)   #Dominance variance proportion
c2_grid      <- c(0.0)   #Shared environment variance proportion
n_sib_grid   <- c(1000, 5000, 10000, 50000)

K_grid <- K_grid_all

out_dir <- "XXX"
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

cat("Additive-only scenario when dominance/shared env variance (d2_grid/c2_grid) set to 0\n")
cat("t_true <- (0.5*h2 + 0.25*d2 + c2)  # Falconer & Mackay, with Vp=1\n")
cat("Warning: estimate t_tet and P(D) from an unascertained sibling sample (not discordant pairs only).\n\n")

#log10 x-axis
scale_x_n_sib_log10 <- function() {
  scale_x_log10(
    breaks = n_sib_grid,
    labels = n_sib_grid
  )
}

#Convert observed Pearson correlation of 0/1 sib phenotypes to tetrachoric rho
#Adapted from: Yengo et al. (https://github.com/loic-yengo/ConvertFSIR_to_liability)
convPearsonToTetrachoric <- function(r_01, K){
  if (!is.finite(r_01)) return(NA_real_)  #guard against NA cor when margin has zero variance
  tau <- qnorm(1 - K)
  r_01 <- max(min(r_01, 0.999999), -0.999999)
  P_11 <- K*K + r_01 * K*(1-K)
  #clamp to a feasible range for numerical stability 
  P_11 <- min(max(P_11, 1e-12), K - 1e-12)
  mod <- optim(
    par = 0,
    fn = function(rho){
      abs(P_11 - fMultivar::pnorm2d(-tau, -tau, rho = rho)[1])
    },
    method = "Brent",
    lower = -1, upper = 1
  )
  return(mod$par)
}

simulate_cell <- function(rep_id, seed, K, h2, d2, c2, n_sib_pairs){
  set.seed(seed)
  #t from Falconer & Mackay (when Vp=1 on liability scale)
  t_true <- (0.5*h2 + 0.25*d2 + c2)

  #Liability threshold and density at threshold
  tau <- qnorm(1 - K) #thresold 
  z   <- dnorm(tau)

  #Variance components on liability scale
  Va  <- h2
  Vd  <- d2
  Vec <- c2
  Ve  <- 1 - (Va + Vd + Vec) #unique envirornment
  if (Ve <= 0) stop("Invalid variance decomposition")

  #Additive genetic values (A)      Cov(A1,A2)=0.5*Va
  Sigma_A <- matrix(c(Va, 0.5*Va,
                      0.5*Va, Va), 2, 2) #covariance matrix, each sibling has variance Va in additive effects, with covariance 0.5Va between the pair
  A_sibs <- mvrnorm(n_sib_pairs, mu = c(0, 0), Sigma = Sigma_A) #multivariate normal draw

  #Dominance genetic values (D): Cov(D1,D2)=0.25*Vd
  if (Vd > 0) {
    Sigma_D <- matrix(c(Vd, 0.25*Vd,
                        0.25*Vd, Vd), 2, 2) #covariance matrix, each sibling has variance Vd in dominance effects, with covariance 0.25Vd between the pair
    D_sibs <- mvrnorm(n_sib_pairs, mu = c(0, 0), Sigma = Sigma_D) #multivariate normal draw
  } else {
    D_sibs <- matrix(0, n_sib_pairs, 2)
  }

  #Shared environment
  if (Vec > 0) {
    c_shared <- rnorm(n_sib_pairs, mean = 0, sd = sqrt(Vec))
    C_sibs <- cbind(c_shared, c_shared)
  } else {
    C_sibs <- matrix(0, n_sib_pairs, 2)
  }

  #Unique environment
  e_sibs <- matrix(rnorm(n_sib_pairs * 2, mean = 0, sd = sqrt(Ve)), n_sib_pairs, 2)

  #Liability and binary phenotype
  l_sibs <- A_sibs + D_sibs + C_sibs + e_sibs #liability
  y_sibs <- 1L * (l_sibs > tau) #binary phenotype

  #differences
  delta_y   <- y_sibs[,1] - y_sibs[,2]
  delta_prs <- A_sibs[,1] - A_sibs[,2]    #oracle additive PRS differences

  #Pearson on observed 0/1 and tetrachoric estimate
  r_01  <- suppressWarnings(cor(y_sibs[,1], y_sibs[,2]))
  t_tet <- convPearsonToTetrachoric(r_01, K)

  #Discordance
  D   <- (y_sibs[,1] != y_sibs[,2])
  P_D <- mean(D)
  nD  <- sum(D)

  #Within-sib R2
  cov_dy_dprs  <- cov(delta_y, delta_prs)
  var_dprs_all <- var(delta_prs)
  R2_within_true_t <- cov_dy_dprs^2 / (z^2 * var_dprs_all * 2 * (1 - t_true))
  R2_within_tet_t  <- if (is.finite(t_tet)) cov_dy_dprs^2 / (z^2 * var_dprs_all * 2 * (1 - t_tet)) else NA_real_

  #Discordant-sib R2
  cov_D <- NA_real_ #missing value constant 
  if (nD >= 2) cov_D <- cov(delta_y[D], delta_prs[D])

  R2_disc_true_t <- NA_real_
  R2_disc_tet_t  <- NA_real_
  if (!is.na(cov_D)) {
    R2_disc_true_t <- (P_D * cov_D)^2 / (z^2 * var_dprs_all * 2 * (1 - t_true))
    if (is.finite(t_tet)) {
      R2_disc_tet_t  <- (P_D * cov_D)^2 / (z^2 * var_dprs_all * 2 * (1 - t_tet))
    }
  }

  #Within-sib heritability benchmark
  h2w_true_t <- h2 / (2 - 2*t_true)
  h2w_tet_t  <- if (is.finite(t_tet)) h2 / (2 - 2*t_tet) else NA_real_

  ratio_within_true_t <- R2_within_true_t / h2w_true_t
  ratio_disc_true_t   <- R2_disc_true_t   / h2w_true_t

  ratio_within_tet_t <- R2_within_tet_t / h2w_tet_t
  ratio_disc_tet_t   <- R2_disc_tet_t   / h2w_tet_t

  tibble(
    rep = rep_id,
    seed = seed,
    K = K,
    h2 = h2,
    d2 = d2,
    c2 = c2,
    n_sib = n_sib_pairs,
    tau = tau,
    z = z,
    r_01 = r_01,
    t_tet = t_tet,
    t_true = t_true,
    P_D = P_D,
    nD = nD,
    R2_within_true_t = R2_within_true_t,
    R2_within_tet_t  = R2_within_tet_t,
    R2_disc_true_t   = R2_disc_true_t,
    R2_disc_tet_t    = R2_disc_tet_t,
    h2w_true_t = h2w_true_t,
    h2w_tet_t  = h2w_tet_t,
    ratio_within_true_t = ratio_within_true_t,
    ratio_disc_true_t   = ratio_disc_true_t,
    ratio_within_tet_t  = ratio_within_tet_t,
    ratio_disc_tet_t    = ratio_disc_tet_t
  )
}

#simulation grid
grid <- tidyr::crossing(
  rep = 1:R,
  K = K_grid,
  h2 = h2_grid,
  d2 = d2_grid,
  c2 = c2_grid,
  n_sib = n_sib_grid
) %>%
  mutate(seed = 1e6 + rep * 1000 + row_number())

cat(sprintf("\nRunning simulation: %d cells (%d reps x %d K x %d h2 x %d n_sib)\n",
            nrow(grid), R, length(K_grid), length(h2_grid), length(n_sib_grid)))

results_long <- purrr::pmap_dfr(
  list(grid$rep, grid$seed, grid$K, grid$h2, grid$d2, grid$c2, grid$n_sib),
  simulate_cell
)

#full results
write_csv(results_long, file.path(out_dir, "results_long.csv"))
saveRDS(results_long, file.path(out_dir, "results_long.rds"))

#parameters
params <- list(
  R = R,
  K_grid_all = K_grid_all,
  K_grid_sweep = K_grid_sweep,
  h2_grid = h2_grid,
  d2_grid = d2_grid,
  c2_grid = c2_grid,
  n_sib_grid = n_sib_grid,
  out_dir = out_dir,
  fig_dir = fig_dir,
  tab_dir = tab_dir
)
writeLines(toJSON(params, pretty = TRUE, auto_unbox = TRUE), con = file.path(out_dir, "run_parameters.json"))
write_csv(
  tibble::tibble(parameter = names(params), value = vapply(params, function(x) paste(x, collapse = ", "), character(1))),
  file.path(out_dir, "run_parameters.csv")
)

#Summary
summ <- results_long %>%
  group_by(K, h2, d2, c2, n_sib) %>%
  summarise(
    n_reps = n(),
    mean_R2_within_true_t = mean(R2_within_true_t, na.rm = TRUE),
    sd_R2_within_true_t   = sd(R2_within_true_t,   na.rm = TRUE),
    mean_R2_disc_true_t = mean(R2_disc_true_t, na.rm = TRUE),
    sd_R2_disc_true_t   = sd(R2_disc_true_t,   na.rm = TRUE),

    mean_ratio_within_true_t = mean(ratio_within_true_t, na.rm = TRUE),
    sd_ratio_within_true_t   = sd(ratio_within_true_t,   na.rm = TRUE),
    mean_ratio_disc_true_t = mean(ratio_disc_true_t, na.rm = TRUE),
    sd_ratio_disc_true_t   = sd(ratio_disc_true_t,   na.rm = TRUE),

    mean_R2_within_tet_t = mean(R2_within_tet_t, na.rm = TRUE),
    sd_R2_within_tet_t   = sd(R2_within_tet_t,   na.rm = TRUE),
    mean_R2_disc_tet_t = mean(R2_disc_tet_t, na.rm = TRUE),
    sd_R2_disc_tet_t   = sd(R2_disc_tet_t,   na.rm = TRUE),

    mean_ratio_within_tet_t = mean(ratio_within_tet_t, na.rm = TRUE),
    sd_ratio_within_tet_t   = sd(ratio_within_tet_t,   na.rm = TRUE),
    mean_ratio_disc_tet_t = mean(ratio_disc_tet_t, na.rm = TRUE),
    sd_ratio_disc_tet_t   = sd(ratio_disc_tet_t,   na.rm = TRUE),

    mean_P_D = mean(P_D, na.rm = TRUE),
    sd_P_D   = sd(P_D,   na.rm = TRUE),
    mean_nD = mean(nD, na.rm = TRUE),
    sd_nD   = sd(nD,   na.rm = TRUE),

    mean_t_true = mean(t_true, na.rm = TRUE),
    sd_t_true   = sd(t_true,   na.rm = TRUE),
    mean_t_tet = mean(t_tet, na.rm = TRUE),
    sd_t_tet   = sd(t_tet,   na.rm = TRUE),

    mean_h2w_true_t = mean(h2w_true_t, na.rm = TRUE),
    sd_h2w_true_t   = sd(h2w_true_t,   na.rm = TRUE),
    mean_h2w_tet_t = mean(h2w_tet_t, na.rm = TRUE),
    sd_h2w_tet_t   = sd(h2w_tet_t,   na.rm = TRUE),

    .groups = "drop"
  ) %>%
  mutate(
    n_sib_f = factor(n_sib, levels = n_sib_grid),
    K_f_sweep = factor(K, levels = K_grid_sweep),
    h2_f = factor(h2, levels = h2_grid)
  )


write_csv(summ, file.path(out_dir, "summary_full.csv"))
saveRDS(summ, file.path(out_dir, "summary_full.rds"))


                 
write_csv(
  summ %>% select(K, h2, d2, c2, n_sib, n_reps, mean_ratio_within_true_t, sd_ratio_within_true_t),
  file.path(tab_dir, "Table_ratio_R2ws_true_t.csv")
)

write_csv(
  summ %>% select(K, h2, d2, c2, n_sib, n_reps, mean_ratio_disc_true_t, sd_ratio_disc_true_t),
  file.path(tab_dir, "Table_ratio_R2dsp_true_t.csv")
)

tab_nD_by_K_h2_n <- summ %>%
  filter(K %in% K_grid_sweep) %>%
  select(K, h2, n_sib, n_reps, mean_nD, sd_nD, mean_P_D, sd_P_D) %>%
  arrange(K, h2, n_sib)

write_csv(tab_nD_by_K_h2_n, file.path(tab_dir, "Table_nD_by_K_h2_n.csv"))

#Within-sib ratio 
fig_ratio_ws_true <- summ %>%
  filter(K %in% K_grid_sweep) %>%
  mutate(K_f = factor(K, levels = K_grid_sweep))

p_ratio_ws_true <- ggplot(
  fig_ratio_ws_true,
  aes(x = n_sib, y = mean_ratio_within_true_t, group = 1)
) +
  geom_hline(yintercept = 1.0, linetype = "dashed") +
  geom_point(size = 2.4) +
  geom_errorbar(aes(
    ymin = mean_ratio_within_true_t - sd_ratio_within_true_t,
    ymax = mean_ratio_within_true_t + sd_ratio_within_true_t
  ), width = 0) +
  facet_grid(
    rows = vars(h2_f),
    cols = vars(K_f),
    scales = "free_y",
    labeller = labeller(
      K_f  = function(x) paste0("K = ", formatC(as.numeric(as.character(x)), format = "fg", digits = 2)),
      h2_f = function(x) paste0("h² = ", x)
    )
  ) +
  scale_x_n_sib_log10() +
  labs(
    x = "Number of Sibling Pairs",
    y = "Ratio: Within-Sibling R² / h²w"
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "Figure_ratio_R2ws_true_t.png"), p_ratio_ws_true, width = 11, height = 8.5, dpi = 600)
ggsave(file.path(fig_dir, "Figure_ratio_R2ws_true_t.pdf"), p_ratio_ws_true, width = 11, height = 8.5)

#Discordant-sib ratio 
fig_ratio_dsp_true <- summ %>%
  filter(K %in% K_grid_sweep) %>%
  mutate(K_f = factor(K, levels = K_grid_sweep))

p_ratio_dsp_true <- ggplot(
  fig_ratio_dsp_true,
  aes(x = n_sib, y = mean_ratio_disc_true_t, group = 1)
) +
  geom_hline(yintercept = 1.0, linetype = "dashed") +
  geom_point(size = 2.4) +
  geom_errorbar(aes(
    ymin = mean_ratio_disc_true_t - sd_ratio_disc_true_t,
    ymax = mean_ratio_disc_true_t + sd_ratio_disc_true_t
  ), width = 0) +
  facet_grid(
    rows = vars(h2_f),
    cols = vars(K_f),
    scales = "free_y",
    labeller = labeller(
      K_f  = function(x) paste0("K = ", formatC(as.numeric(as.character(x)), format = "fg", digits = 2)),
      h2_f = function(x) paste0("h² = ", x)
    )
  ) +
  scale_x_n_sib_log10() +
  labs(
    x = "Number of Sibling Pairs",
    y = "Ratio: Discordant-Sibling R² / h²w"
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "Figure_ratio_R2dsp_true_t.png"), p_ratio_dsp_true, width = 11, height = 8.5, dpi = 600)
ggsave(file.path(fig_dir, "Figure_ratio_R2dsp_true_t.pdf"), p_ratio_dsp_true, width = 11, height = 8.5)

#t_true vs t_tet (fixed n_sib=10k)
fig_t_vs_ttet <- summ %>%
  filter(n_sib == 10000, K %in% K_grid_sweep)

p_t <- ggplot(fig_t_vs_ttet, aes(x = mean_t_true, y = mean_t_tet, group = h2_f, shape = h2_f, colour = h2_f)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(size = 2.4) +
  geom_errorbar(aes(ymin = mean_t_tet - sd_t_tet, ymax = mean_t_tet + sd_t_tet), width = 0.01) +
  facet_wrap(~ K, scales = "free") +
  labs(
    x = "Intraclass Correlation (t)",
    y = "Tetrachoric t estimate",
    colour = expression(h^2),
    shape  = expression(h^2)
  ) +
  theme_bw()

ggsave(file.path(fig_dir, "Figure_t_true_vs_t_tet.png"), p_t, width = 11, height = 5.5, dpi = 600)
ggsave(file.path(fig_dir, "Figure_t_true_vs_t_tet.pdf"), p_t, width = 11, height = 5.5)

#P(D) vs K 
fig_PD <- summ %>% filter(K %in% K_grid_all)

p_PD <- ggplot(fig_PD, aes(x = K, y = mean_P_D, group = h2_f, shape = h2_f, colour = h2_f)) +
  geom_line() +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = mean_P_D - sd_P_D, ymax = mean_P_D + sd_P_D), width = 0) +
  facet_wrap(~ n_sib, scales = "free_y",
             labeller = labeller(n_sib = function(x) paste0("Number of Sibling Pairs = ", x))) +
  scale_x_log10(breaks = K_grid_all) +
  labs(
    x = "K",
    y = "Probability of Discordance",
    colour = expression(h^2),
    shape  = expression(h^2)
  ) +
  theme_bw()

ggsave(file.path(fig_dir, "Figure_PD_vs_K_log10.png"), p_PD, width = 11, height = 7, dpi = 600)
ggsave(file.path(fig_dir, "Figure_PD_vs_K_log10.pdf"), p_PD, width = 11, height = 7)



#Within-sibling R2 vs n_sib 
fig_metric_r2ws_tet <- summ %>% filter(K %in% K_grid_sweep)

ref_band_ws_tet <- fig_metric_r2ws_tet %>%
  distinct(K, h2, h2_f, n_sib, mean_h2w_tet_t, sd_h2w_tet_t) %>%
  mutate(
    h2w_lo = mean_h2w_tet_t - sd_h2w_tet_t,
    h2w_hi = mean_h2w_tet_t + sd_h2w_tet_t
  )

p_r2ws_tet <- ggplot(
  fig_metric_r2ws_tet,
  aes(x = n_sib, y = mean_R2_within_tet_t, group = h2_f, shape = h2_f, colour = h2_f)
) +
  geom_ribbon(
    data = ref_band_ws_tet,
    aes(x = n_sib, ymin = h2w_lo, ymax = h2w_hi, group = h2_f, fill = h2_f),
    inherit.aes = FALSE,
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(
    data = ref_band_ws_tet,
    aes(x = n_sib, y = mean_h2w_tet_t, group = h2_f, colour = h2_f),
    inherit.aes = FALSE,
    linetype = "dashed",
    size = 0.6
  ) +
  geom_point(size = 2.4) +
  geom_errorbar(aes(
    ymin = mean_R2_within_tet_t - sd_R2_within_tet_t,
    ymax = mean_R2_within_tet_t + sd_R2_within_tet_t
  ), width = 0) +
  facet_wrap(~ K, scales = "free_y",
             labeller = labeller(K = function(x) paste0("K = ", formatC(as.numeric(x), format="fg", digits=2)))) +
  scale_x_n_sib_log10() +
  labs(
    x = "Number of Sibling Pairs",
    y = "Within-Sibling R²",
    colour = expression(h^2),
    shape  = expression(h^2),
    fill   = expression(h^2)
  ) +
  theme_bw()

ggsave(file.path(fig_dir, "Figure_R2_within_vs_n_tet_t.png"), p_r2ws_tet, width = 11, height = 5.5, dpi = 600)
ggsave(file.path(fig_dir, "Figure_R2_within_vs_n_tet_t.pdf"), p_r2ws_tet, width = 11, height = 5.5)

#Discordant-sibling R2 vs n_sib
fig_metric_r2dsp_tet <- summ %>% filter(K %in% K_grid_sweep)

ref_band_dsp_tet <- fig_metric_r2dsp_tet %>%
  distinct(K, h2, h2_f, n_sib, mean_h2w_tet_t, sd_h2w_tet_t) %>%
  mutate(
    h2w_lo = mean_h2w_tet_t - sd_h2w_tet_t,
    h2w_hi = mean_h2w_tet_t + sd_h2w_tet_t
  )

p_r2dsp_tet <- ggplot(
  fig_metric_r2dsp_tet,
  aes(x = n_sib, y = mean_R2_disc_tet_t, group = h2_f, shape = h2_f, colour = h2_f)
) +
  geom_ribbon(
    data = ref_band_dsp_tet,
    aes(x = n_sib, ymin = h2w_lo, ymax = h2w_hi, group = h2_f, fill = h2_f),
    inherit.aes = FALSE,
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(
    data = ref_band_dsp_tet,
    aes(x = n_sib, y = mean_h2w_tet_t, group = h2_f, colour = h2_f),
    inherit.aes = FALSE,
    linetype = "dashed",
    size = 0.6
  ) +
  geom_point(size = 2.4) +
  geom_errorbar(aes(
    ymin = mean_R2_disc_tet_t - sd_R2_disc_tet_t,
    ymax = mean_R2_disc_tet_t + sd_R2_disc_tet_t
  ), width = 0) +
  facet_wrap(~ K, scales = "free_y",
             labeller = labeller(K = function(x) paste0("K = ", formatC(as.numeric(x), format="fg", digits=2)))) +
  scale_x_n_sib_log10() +
  labs(
    x = "Number of Sibling Pairs",
    y = "Discordant-Sibling R²",
    colour = expression(h^2),
    shape  = expression(h^2),
    fill   = expression(h^2)
  ) +
  theme_bw()

ggsave(file.path(fig_dir, "Figure_R2_discordant_vs_n_tet_t.png"), p_r2dsp_tet, width = 11, height = 5.5, dpi = 600)
ggsave(file.path(fig_dir, "Figure_R2_discordant_vs_n_tet_t.pdf"), p_r2dsp_tet, width = 11, height = 5.5)

#Discordant pairs vs cohort size 
fig_nD <- summ %>% filter(K %in% K_grid_sweep)

p_nD <- ggplot(
  fig_nD,
  aes(x = n_sib, y = mean_nD, group = h2_f, shape = h2_f, colour = h2_f)
) +
  geom_line() +
  geom_point(size = 2.4) +
  geom_errorbar(aes(
    ymin = mean_nD - sd_nD,
    ymax = mean_nD + sd_nD
  ), width = 0) +
  facet_wrap(~ K, scales = "free_y",
             labeller = labeller(K = function(x) paste0("K = ", formatC(as.numeric(x), format="fg", digits=2)))) +
  scale_x_n_sib_log10() +
  labs(
    x = "Number of Sibling Pairs",
    y = "Number of Discordant Pairs",
    colour = expression(h^2),
    shape  = expression(h^2)
  ) +
  theme_bw()

ggsave(file.path(fig_dir, "Figure_nD_vs_n.png"), p_nD, width = 11, height = 5.5, dpi = 600)
ggsave(file.path(fig_dir, "Figure_nD_vs_n.pdf"), p_nD, width = 11, height = 5.5)

cat("\nAll outputs saved to: ", normalizePath(out_dir), "\n")
cat(" results_long.csv / results_long.rds (full numerical results)\n")
cat(" summary_full.csv / summary_full.rds (summary)\n")
cat(" run_parameters.json / run_parameters.csv\n")
cat(" figures: PNG (dpi=600) + PDF\n")
cat(" tables: CSV\n")
cat(" Table_nD_by_K_h2_n.csv (discordant counts lookup)\n")
