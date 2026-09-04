# Module C (Paper C): design-statistic ladder + covariance-aware inference layer,
# pinned to the audit pipeline (Tables 4-5). Design statistics are deterministic;
# the covariance-aware pilot and its sizes depend on the wild bootstrap (fixed
# seed) and are checked with wider tolerance.

test_that("block design: Gamma = 0 exactly (Prop. prop-gamma0)", {
  N <- 10; T <- 6; g <- 4
  unit <- rep(1:N, each = T); time <- rep(1:T, times = N)
  ft <- ifelse(unit <= 4, g, NA)
  r <- twfe_design(unit, time, ft)
  expect_lt(r$statistic$Gamma, 1e-8)
  expect_lt(r$statistic$Gamma_cmb, 1e-8)
  expect_lt(r$statistic$Gamma_gt, 1e-8)
  expect_equal(r$statistic$neg_share, 0)
  expect_identical(r$verdict, "CERTIFIED")
  expect_true(any(grepl("block", r$notes)))
})

test_that("three-cohort no-reservoir design (Gamma = 1.54); nested ladder", {
  N <- 48; T <- 12
  unit <- rep(1:N, each = T); time <- rep(1:T, times = N)
  ft <- ifelse(unit <= 16, 3, ifelse(unit <= 32, 7, 11))
  r <- twfe_design(unit, time, ft); st <- r$statistic
  expect_equal(st$Gamma, 1.54, tolerance = 0.02 / 1.54)
  expect_equal(st$neg_share, 0.11, tolerance = 0.01 / 0.11)
  expect_lte(st$Gamma_coh, st$Gamma_cmb + 1e-9)
  expect_lte(st$Gamma_evt, st$Gamma_cmb + 1e-9)
  expect_lte(st$Gamma_cmb, st$Gamma_gt + 1e-9)
  expect_equal(st$Gamma_gt, st$Gamma, tolerance = 1e-10)
  expect_identical(r$verdict, "INCONCLUSIVE")
  expect_equal(r$breakdown, PD$.eta_dagger(0.05, 0.05) / st$Gamma_gt,
               tolerance = 1e-10)
})

test_that("castle design statistics (Table 4)", {
  castle <- read_panel("castle_panel.csv")
  r <- twfe_design(castle$uid, castle$tid, castle$ft); st <- r$statistic
  expect_equal(r$design$n, 550)
  expect_equal(c(r$design$N, r$design$T), c(50, 11))
  expect_equal(st$Gamma, 0.2114, tolerance = 1e-3 / 0.2114)
  expect_equal(st$Gamma_cmb, 0.198, tolerance = 5e-3 / 0.198)
  expect_equal(st$Gamma_gt, st$Gamma, tolerance = 1e-10)
  expect_equal(st$Gamma_evt, 0.168, tolerance = 5e-3 / 0.168)
  expect_equal(st$Gamma_coh, 0.142, tolerance = 5e-3 / 0.142)
  expect_equal(st$neg_share, 0)
  expect_equal(st$N1, 95)
  expect_equal(st$n_w, 34.7127, tolerance = 1e-4)
})

test_that("divorce design statistics (always-treated dropped; Table 4)", {
  divorce <- read_panel("divorce_panel.csv")
  r <- twfe_design(divorce$uid, divorce$tid, divorce$ft); st <- r$statistic
  expect_equal(r$design$N, 49)          # 51 - 2 always-treated
  expect_equal(r$design$n, 1323)
  expect_equal(st$Gamma, 0.6433, tolerance = 2e-3 / 0.6433)
  expect_equal(st$Gamma_cmb, 0.562, tolerance = 5e-3 / 0.562)
  expect_equal(st$Gamma_gt, st$Gamma, tolerance = 1e-10)
  expect_equal(st$Gamma_evt, 0.468, tolerance = 5e-3 / 0.468)
  expect_equal(st$Gamma_coh, 0.381, tolerance = 5e-3 / 0.381)
  expect_equal(st$neg_share, 0.0115, tolerance = 2e-3)
  expect_equal(st$N1, 522)
  expect_true(any(grepl("always-treated", r$notes)))
})

test_that("castle inference: projected-norm upper bound certifies", {
  castle <- read_panel("castle_panel.csv")
  r <- twfe_adequacy(castle$y, castle$uid, castle$tid, castle$ft,
                     bootstrap = 299L, seed = 20260715L)
  st <- r$statistic
  expect_equal(st$beta, 0.081812, tolerance = 1e-3)
  expect_equal(st$sigma, 0.186992, tolerance = 1e-3)
  expect_equal(st$psi_hat, 3.3721, tolerance = 1e-3)
  expect_equal(st$psi_direct, 3.3721, tolerance = 1e-3)
  expect_equal(st$psi_ar1, 1.3358, tolerance = 1e-3)
  expect_equal(st$Gamma_cmb_CR, 0.1079, tolerance = 5e-3 / 0.1079)
  expect_equal(st$pilot_cmb, 0, tolerance = 1e-6)          # noise floor
  expect_equal(st$size_cmb, 0.05, tolerance = 2e-3)
  expect_equal(st$Gamma_gt, st$Gamma, tolerance = 1e-10)
  expect_equal(st$pilot_gt, 0.74787382, tolerance = 1e-7)
  expect_equal(r$eta, 0.08610202, tolerance = 1e-7)
  expect_equal(st$size_gt, 0.05085, tolerance = 2e-3)
  expect_lt(st$K_upper_gt, r$threshold)
  expect_gt(st$K_upper_gt_hc3, r$threshold)
  expect_equal(st$eta_directional, -0.0134423, tolerance = 1e-5)
  expect_equal(st$size_directional, 0.0500207, tolerance = 1e-5)
  expect_equal(st$directional_alignment, -0.164229, tolerance = 1e-5)
  expect_equal(st$sign_reversal_rms, 0.386972, tolerance = 1e-5)
  expect_equal(st$size_realized, st$size_directional, tolerance = 1e-12)
  expect_identical(r$verdict, "CERTIFIED")
  expect_true(!is.null(st$boot) && !is.null(st$boot$norm))
})

test_that("divorce point envelope is large but norm bound is inconclusive", {
  divorce <- read_panel("divorce_panel.csv")
  r <- twfe_adequacy(divorce$y, divorce$uid, divorce$tid, divorce$ft,
                     bootstrap = 299L, seed = 20260715L)
  st <- r$statistic
  expect_equal(st$sigma, 0.1961, tolerance = 1e-3)
  expect_equal(st$psi_hat, 2.91718, tolerance = 2e-3)
  expect_equal(st$psi_direct, 2.91718, tolerance = 2e-3)
  expect_equal(st$psi_ar1, 1.61958, tolerance = 2e-3)
  expect_equal(st$Gamma_cmb_CR, 0.32925, tolerance = 5e-3 / 0.32925)
  expect_equal(st$eta_directional, 1.629204, tolerance = 1e-5)
  expect_equal(st$size_directional, 0.370579, tolerance = 1e-5)
  expect_equal(st$directional_alignment, 0.567753, tolerance = 1e-5)
  expect_equal(st$sign_reversal_rms, 0.0245369, tolerance = 1e-5)
  expect_equal(st$pilot_cmb, 5.29888, tolerance = 1e-3)
  expect_equal(st$size_cmb, 0.414877, tolerance = 2e-3)
  expect_equal(st$Gamma_gt, st$Gamma, tolerance = 1e-10)
  expect_equal(st$pilot_gt, 6.59941945, tolerance = 1e-7)
  expect_equal(r$eta, 2.48546963, tolerance = 1e-7)
  expect_gt(st$size_gt, 0.65)
  expect_lt(st$K_lower_gt, r$threshold)
  expect_gt(st$K_upper_gt, r$threshold)
  expect_gt(st$size_coh, 0.10); expect_gt(st$size_evt, 0.10)
  expect_gte(st$size_cmb, st$size_coh - 1e-9)               # nesting
  expect_identical(r$verdict, "INCONCLUSIVE")
  expect_true(any(grepl("point pilot", r$notes)))
  expect_true(any(grepl("fixed-T", r$notes)))
})

test_that("minimum-wage headline reproduces the JAE point diagnostics", {
  d <- get_dataset("minimum_wage")
  r <- twfe_adequacy(d$y, d$uid, d$tid, d$ft, controls = "never",
                     bootstrap = 19L, seed = 20260828L)
  st <- r$statistic
  expect_equal(st$beta, -0.03660863, tolerance = 1e-7)
  expect_equal(st$Gamma, 0.26964309, tolerance = 1e-7)
  expect_equal(st$Gamma_cmb, 0.25728013, tolerance = 1e-7)
  expect_equal(st$Gamma_gt, st$Gamma, tolerance = 1e-10)
  expect_equal(st$neg_share, 0, tolerance = 1e-12)
  expect_equal(st$psi_direct, 1.38489247, tolerance = 1e-7)
  expect_equal(st$pilot_cmb, 6.46943999, tolerance = 1e-7)
  expect_equal(st$size_cmb, 0.29304460, tolerance = 1e-7)
  expect_equal(st$pilot_gt, 6.46015317, tolerance = 1e-7)
  expect_equal(st$size_gt, 0.31599274, tolerance = 1e-7)
  expect_equal(st$att_target, -0.05158099, tolerance = 1e-7)
  expect_gt(st$K_lower_gt, r$threshold)
  expect_equal(st$eta_directional, 1.37026379, tolerance = 1e-7)
  expect_equal(st$size_directional, 0.27812971, tolerance = 1e-7)
  expect_equal(st$directional_alignment, 0.92476056, tolerance = 1e-7)
  expect_equal(st$sign_reversal_rms, 0.13576699, tolerance = 1e-7)
  expect_identical(r$verdict, "FLAGGED")
})

test_that("exchangeable psi identity (Thm thm-cluster(a))", {
  castle <- read_panel("castle_panel.csv")
  sc <- PD$.staggered_codes(castle$uid, castle$tid, castle$ft)
  D <- PD$.treatment_indicator(castle$tid, castle$ft)
  Dt <- PD$.twoway_demean_codes(D, sc$uid, sc$tid, sc$N, sc$T)
  for (rho_c in c(0.3, 0.5, 0.8)) {
    psi <- PD$.psi_parametric(Dt, sc$uid, sc$tid, rho_c, kind = "exchangeable")
    expect_equal(psi, 1 - rho_c, tolerance = 1e-9)
  }
  expect_equal(PD$.psi_parametric(Dt, sc$uid, sc$tid, 0, kind = "ar1"), 1,
               tolerance = 1e-9)
})

test_that("validation and rendering", {
  unit <- rep(1:6, each = 4); time <- rep(1:4, times = 6)
  expect_error(twfe_design(unit, time, rep(NA, 24)), "treated")
  badft <- ifelse(seq_len(24) %% 3 == 0, 2, 3)
  expect_error(twfe_design(unit, time, badft), "varies within unit")
  expect_error(twfe_design(unit, time, ifelse(unit <= 3, 2.5, NA)),
               "not an observed period")

  castle <- read_panel("castle_panel.csv")
  out <- paste(utils::capture.output(
    print(twfe_adequacy(castle$y, castle$uid, castle$tid, castle$ft, bootstrap = 99L))),
    collapse = "\n")
  expect_match(out, "TWFE Heterogeneity")
  expect_match(out, "restricted ladder")
  expect_match(out, "Point worst-case envelopes")
  expect_match(out, "Boundary-robust group-time K interval")
  expect_match(out, "Directional plug-in")
  expect_false(grepl("Realized-profile", out, fixed = TRUE))
  expect_match(out, "VERDICT: FORMALLY CERTIFIED")
  outd <- paste(utils::capture.output(
    print(twfe_design(castle$uid, castle$tid, castle$ft))), collapse = "\n")
  expect_match(outd, "negative-weight share")
  expect_error(twfe_adequacy(castle$y, castle$uid, castle$tid, castle$ft,
                             cluster = "bad"), "arg")
  expect_error(twfe_adequacy(castle$y, castle$uid, castle$tid, castle$ft,
                             controls = "bad"), "arg")
  expect_error(twfe_adequacy(castle$y, castle$uid, castle$tid, castle$ft,
                             psi = -1), "positive")
  expect_error(twfe_adequacy(castle$y, castle$uid, castle$tid, castle$ft,
                             gamma = 0.6), "gamma")
})
