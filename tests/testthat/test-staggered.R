# Module C (Paper C): design statistic + inference layer, pinned to the
# verified audit pipeline (Paper C Table tab-audit), spec 7.3-7.6.

test_that("block design: Gamma = 0 exactly (Prop. prop-gamma0)", {
  N <- 10; T <- 6; g <- 4
  unit <- rep(1:N, each = T); time <- rep(1:T, times = N)
  ft <- ifelse(unit <= 4, g, NA)
  r <- twfe_design(unit, time, ft)
  expect_lt(r$statistic$Gamma, 1e-8)
  expect_equal(r$statistic$neg_share, 0)
  expect_identical(r$verdict, "CERTIFIED")
  expect_true(any(grepl("block", r$notes)))
})

test_that("three-cohort no-reservoir design (Paper C sims: Gamma = 1.54)", {
  N <- 48; T <- 12
  unit <- rep(1:N, each = T); time <- rep(1:T, times = N)
  ft <- ifelse(unit <= 16, 3, ifelse(unit <= 32, 7, 11))
  r <- twfe_design(unit, time, ft)
  expect_equal(r$statistic$Gamma, 1.54, tolerance = 0.02 / 1.54)
  expect_equal(r$statistic$neg_share, 0.11, tolerance = 0.01 / 0.11)
  expect_identical(r$verdict, "INCONCLUSIVE")
  expect_equal(r$breakdown, PD$.eta_dagger(0.05, 0.05) / r$statistic$Gamma,
               tolerance = 1e-10)
})

test_that("reference case 3: design statistics (spec 7.3)", {
  castle <- read_panel("castle_panel.csv")
  r <- twfe_design(castle$uid, castle$tid, castle$ft)
  expect_equal(r$design$n, 550)
  expect_equal(c(r$design$N, r$design$T), c(50, 11))
  expect_equal(r$statistic$Gamma, 0.211415, tolerance = 1e-4)
  expect_equal(r$statistic$neg_share, 0)
  expect_equal(r$statistic$N1, 95)
  expect_equal(r$statistic$n_w, 34.7127, tolerance = 1e-4)

  divorce <- read_panel("divorce_panel.csv")
  rd <- twfe_design(divorce$uid, divorce$tid, divorce$ft)
  expect_equal(rd$design$n, 1377)
  expect_equal(rd$statistic$Gamma, 0.856801, tolerance = 1e-4)
  expect_equal(rd$statistic$neg_share, 0.072917, tolerance = 1e-4 / 0.072917)
  expect_equal(rd$statistic$N1, 576)
})

test_that("reference cases 4-6: inference layer (audit pipeline values)", {
  castle <- read_panel("castle_panel.csv")
  divorce <- read_panel("divorce_panel.csv")

  r <- twfe_adequacy(castle$y, castle$uid, castle$tid, castle$ft, cluster = "iid")
  st <- r$statistic
  expect_equal(st$beta, 0.081812, tolerance = 1e-3)      # spec 7.6 sanity
  expect_equal(st$sigma, 0.186992, tolerance = 1e-3)
  expect_equal(st$att_bar, 0.109355, tolerance = 1e-3)
  expect_equal(st$cohort_sd_raw, 0.053850, tolerance = 1e-3)
  expect_equal(r$eta, -0.014088, tolerance = 1e-4 / 0.014088)
  expect_equal(r$implied_size, 0.050023, tolerance = 1e-4 / 0.050023)
  expect_equal(st$eta_worst_raw, 0.358708, tolerance = 1e-3)
  expect_equal(st$cohort_sd_shrunk, 0)                    # noise swamps dispersion
  expect_equal(st$eta_worst, 0)
  expect_identical(r$verdict, "CERTIFIED")

  rc <- twfe_adequacy(castle$y, castle$uid, castle$tid, castle$ft)
  stc <- rc$statistic
  expect_equal(stc$rho_ar1, 0.226384, tolerance = 1e-3)
  expect_equal(stc$psi_hat, 1.335767, tolerance = 1e-3)
  expect_equal(stc$Gamma_CR, 0.182924, tolerance = 1e-3)
  expect_equal(rc$eta, -0.012190, tolerance = 1e-4 / 0.012190)
  expect_equal(rc$implied_size, 0.050017, tolerance = 1e-4 / 0.050017)
  expect_identical(rc$verdict, "CERTIFIED")
  expect_equal(stc$psi_driven, 3.372112, tolerance = 1e-2)
  expect_true(any(grepl("cross-check", rc$notes)))

  rd <- twfe_adequacy(divorce$y, divorce$uid, divorce$tid, divorce$ft,
                      cluster = "iid")
  std <- rd$statistic
  expect_equal(std$sigma, 0.198134, tolerance = 1e-3)
  expect_equal(std$att_bar, -0.077839, tolerance = 1e-3)
  expect_equal(std$cohort_sd_raw, 0.305575, tolerance = 1e-3)
  expect_equal(rd$eta, 2.223177, tolerance = 1e-3)
  expect_equal(rd$implied_size, 0.603821, tolerance = 1e-3 / 0.603821)
  expect_equal(std$eta_worst_raw, 12.405, tolerance = 2e-3)
  expect_equal(std$cohort_sd_shrunk, 0.30107, tolerance = 2e-3 / 0.30107)
  expect_equal(std$shrink_factor, 0.985, tolerance = 2e-3 / 0.985)
  expect_identical(rd$verdict, "FLAGGED")

  rdc <- twfe_adequacy(divorce$y, divorce$uid, divorce$tid, divorce$ft)
  stdc <- rdc$statistic
  expect_equal(stdc$rho_ar1, 0.292743, tolerance = 1e-3)
  expect_equal(stdc$psi_hat, 1.615212, tolerance = 1e-3)
  expect_equal(stdc$Gamma_CR, 0.674163, tolerance = 1e-3)
  expect_equal(rdc$eta, 1.749280, tolerance = 1e-3)
  expect_equal(rdc$implied_size, 0.416671, tolerance = 1e-3 / 0.416671)
  expect_identical(rdc$verdict, "FLAGGED")
  expect_true(any(grepl("upper bound", rdc$notes)))
  expect_true(any(grepl("Paper C", rdc$notes)))
})

test_that("exchangeable psi identity (Thm thm-cluster(a))", {
  castle <- read_panel("castle_panel.csv")
  sc <- PD$.staggered_codes(castle$uid, castle$tid, castle$ft)
  D <- PD$.treatment_indicator(castle$tid, castle$ft)
  Dt <- PD$.twoway_demean_codes(D, sc$uid, sc$tid, sc$N, sc$T)
  for (rho_c in c(0.3, 0.5, 0.8)) {
    psi <- PD$.psi_parametric(Dt, sc$uid, sc$tid, rho_c, kind = "exchangeable")
    expect_equal(psi, 1 - rho_c, tolerance = 1e-9)   # d_i'1 = 0 makes this exact
  }
  expect_equal(PD$.psi_parametric(Dt, sc$uid, sc$tid, 0, kind = "ar1"), 1,
               tolerance = 1e-9)
})

test_that("user-supplied cohort effects override the internal pilot", {
  N <- 30; T <- 8
  unit <- rep(1:N, each = T); time <- rep(1:T, times = N)
  ft <- ifelse(unit <= 8, 3, ifelse(unit <= 16, 6, NA))
  k <- seq_along(unit)
  y <- 0.4 * sin(1.1 * k) + 0.05 * unit
  r <- twfe_adequacy(y, unit, time, ft, cluster = "iid",
                     cohort_effects = c(0.5, 0.5), cohort_ses = c(0.1, 0.1))
  expect_equal(r$statistic$cohort_sd_raw, 0)
  expect_equal(r$eta, 0, tolerance = 1e-10)
  expect_identical(r$verdict, "CERTIFIED")
  expect_false(any(grepl("order-of-magnitude", r$notes)))
  r2 <- twfe_adequacy(y, unit, time, ft, cluster = "iid")
  expect_true(any(grepl("order-of-magnitude", r2$notes)))
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
    print(twfe_adequacy(castle$y, castle$uid, castle$tid, castle$ft))),
    collapse = "\n")
  expect_match(out, "TWFE Heterogeneity \\(Paper C\\)")
  expect_match(out, "Gamma")
  expect_match(out, "VERDICT: CERTIFIED")
  outd <- paste(utils::capture.output(
    print(twfe_design(castle$uid, castle$tid, castle$ft))), collapse = "\n")
  expect_match(outd, "negative-weight share")
})
