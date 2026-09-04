# Module B: threshold constants, V-Dem, repeated-report twins, simulation-range
# locks, and honesty machinery. Breakdown expectations are the
# exact-inversion values (the stale linear ones must NOT reappear).

vdem_spec <- function(v, xcol, sdcol) {
  keep <- stats::complete.cases(v$ly, v[[xcol]], v[[sdcol]])
  list(unit = v$iso[keep], time = v$year[keep], y = v$ly[keep],
       x = v[[xcol]][keep], sd = v[[sdcol]][keep])
}

test_that("threshold constants (measurement-error article)", {
  expect_equal(PD$.eta_dagger(0.05, 0.05), 0.652, tolerance = 5e-4 / 0.652)
  expect_equal(PD$.eta_dagger(0.05, 0.01), 0.295, tolerance = 5e-4 / 0.295)
  expect_equal(PD$.eta_quad(0.05, 0.05), 0.661, tolerance = 5e-4 / 0.661)
  expect_equal(PD$.eta_quad(0.05, 0.01), 0.296, tolerance = 1e-3 / 0.296)
  z <- stats::qnorm(0.975)
  expect_equal(z * stats::dnorm(z), 0.11455, tolerance = 2e-5 / 0.11455)
  expect_gt(PD$.eta_quad(0.05, 0.05), PD$.eta_dagger(0.05, 0.05))
  expect_equal(PD$.noncentral_size(0, 0.05), 0.05, tolerance = 1e-12)
  expect_equal(PD$.noncentral_size(-1.3, 0.05), PD$.noncentral_size(1.3, 0.05))
  expect_equal(PD$.noncentral_size(PD$.eta_dagger(0.05, 0.05), 0.05), 0.10,
               tolerance = 1e-9)
  expect_equal(PD$.noncentral_size(3.58, 0.05), 0.9474, tolerance = 5e-5)
  expect_equal(PD$.noncentral_size(3.67, 0.05), 0.9564, tolerance = 5e-5)
  eta_dag <- PD$.eta_dagger(0.05, 0.05)
  z_beta <- stats::qnorm(0.95)
  expect_equal(breakdown_reliability(2, 1, 1), 2 / (2 + eta_dag),
               tolerance = 1e-12)
  expect_equal(certified_breakdown_reliability(2, 1, 1),
               (2 + z_beta) / (2 + z_beta + eta_dag), tolerance = 1e-12)
  expect_gt(certified_breakdown_reliability(2, 1, 1),
            breakdown_reliability(2, 1, 1))
})

test_that("reliability helpers", {
  expect_equal(reliability_from_interval(c(0.1, 0.2), c(0.3, 0.5)), c(0.1, 0.15))
  expect_equal(reliability_from_ratio(0.8, 2), sqrt(0.2) * 2)
  expect_equal(reliability_from_ratio(0.8, 2, scale = "signal"), sqrt(0.25) * 2)
  xrep <- c(-2, -1, 1, 2); zrep <- c(-1.8, -1.2, 0.9, 2.1)
  expect_equal(reliability_from_repeats(xrep, zrep),
               stats::cov(xrep, zrep) / stats::var(xrep))
  expect_equal(reliability_from_repeats(xrep, zrep, method = "equal_variance"),
               1 - stats::var(xrep - zrep) / (2 * stats::var(xrep)))
  expect_error(reliability_from_ratio(1.2, 2))
  expect_error(reliability_from_ratio(0.8, -1))
  expect_error(reliability_from_interval(0.3, 0.2))
  expect_error(reliability_from_repeats(1, 1))
  expect_error(reliability_from_repeats(c(1, 2), 1))
})

test_that("primitive threshold uses exact inversion by default", {
  rho <- 0.2; c2 <- 3; beta0 <- 0.7; sigma <- 1.4
  expected <- beta0^2 * c2^2 * (1 - rho) /
    (sigma^2 * PD$.eta_dagger(0.05, 0.05)^2) - c2
  expect_equal(tau2_crit(rho, c2, beta0, sigma), expected, tolerance = 1e-12)
  z <- stats::qnorm(0.975)
  quad <- beta0^2 * c2^2 * (1 - rho) * z * stats::dnorm(z) /
    (sigma^2 * 0.05) - c2
  expect_equal(tau2_crit(rho, c2, beta0, sigma, method = "quadratic"),
               quad, tolerance = 1e-12)
  expect_gt(tau2_crit(rho, c2, beta0, sigma),
            tau2_crit(rho, c2, beta0, sigma, method = "quadratic"))
})

test_that("reference case 2a: V-Dem two-pole (spec 7.2)", {
  v <- get_dataset("vdem")

  s <- vdem_spec(v, "v2x_polyarchy", "v2x_polyarchy_sd")
  r <- eiv_adequacy(s$y, s$x, s$unit, s$time, sigma_nu = s$sd, pilot = "point")
  expect_equal(r$design$n, 8930)
  expect_equal(r$design$d_K, 221)
  expect_equal(r$statistic$lambda_hat, 0.8984, tolerance = 1e-3 / 0.8984)
  expect_equal(r$statistic$beta_star, 0.06096, tolerance = 2e-3)
  expect_equal(r$statistic$beta_corr, 0.06785, tolerance = 2e-3)
  expect_equal(r$statistic$sigma, 0.32956, tolerance = 2e-3)
  expect_equal(r$design$tau_star2, 127.826, tolerance = 2e-3)
  expect_equal(r$eta, 0.23640, tolerance = 5e-3)
  expect_equal(r$implied_size, 0.05643, tolerance = 5e-4 / 0.05643)
  expect_equal(r$threshold, 0.652, tolerance = 5e-4 / 0.652)
  expect_equal(r$breakdown, 0.762, tolerance = 2e-3 / 0.762)  # fixed point (paper Table 3)
  # point verdict is exactly equivalent to lambda_hat >= breakdown
  expect_identical(r$statistic$lambda_hat >= r$breakdown, r$verdict == "POINT_PASS")
  expect_identical(r$verdict, "POINT_PASS")   # point pilot: a pass, not a certificate
  rc <- eiv_adequacy(s$y, s$x, s$unit, s$time, sigma_nu = s$sd)
  expect_identical(rc$verdict, "CERTIFIED")   # conservative pilot: a certificate
  expect_gt(rc$statistic$eta_upper, rc$eta)
  expect_equal(rc$statistic$breakdown_certified, 0.8513, tolerance = 2e-3)
  expect_equal(rc$statistic$reliability_lower, rc$statistic$lambda_hat)
  expect_equal(rc$statistic$false_certification_bound, 0.05)
  rl <- eiv_adequacy(s$y, s$x, s$unit, s$time, sigma_nu = s$sd,
                     reliability_lower = 0.80, gamma = 0.025,
                     gamma_lambda = 0.025)
  expect_identical(rl$verdict, "FLAGGED")
  expect_equal(rl$statistic$reliability_lower, 0.80)
  expect_equal(rl$statistic$false_certification_bound, 0.05)
  # cluster-robust (country CRVE): paper Table 3 psi_hat = 19.2, still certified
  rcr <- eiv_adequacy(s$y, s$x, s$unit, s$time, sigma_nu = s$sd,
                      pilot = "point", cluster = "crve")
  expect_equal(rcr$statistic$psi_hat, 19.18, tolerance = 1e-2)
  expect_equal(rcr$implied_size, 0.050, tolerance = 1e-3 / 0.050)
  expect_identical(rcr$verdict, "POINT_PASS")
  expect_true("projection_ratio" %in% names(rcr$statistic$cluster))
  expect_equal(rcr$statistic$cluster$projection_ratio, 0.0069,
               tolerance = 2e-4 / 0.0069)

  s <- vdem_spec(v, "v2xlg_legcon", "v2xlg_legcon_sd")
  r <- eiv_adequacy(s$y, s$x, s$unit, s$time, sigma_nu = s$sd, pilot = "point")
  expect_equal(r$design$n, 8529)
  expect_equal(r$statistic$lambda_hat, 0.5472, tolerance = 1e-3 / 0.5472)
  expect_equal(r$eta, 0.8853, tolerance = 5e-3)
  expect_equal(r$implied_size, 0.1435, tolerance = 1e-3 / 0.1435)
  expect_equal(r$breakdown, 0.621, tolerance = 2e-3 / 0.621)  # fixed point (paper Table 3)
  expect_identical(r$verdict, "FLAGGED")
  # the paper's middle case: flagged iid, point pass under country clustering
  rcr <- eiv_adequacy(s$y, s$x, s$unit, s$time, sigma_nu = s$sd,
                      pilot = "point", cluster = "crve")
  expect_equal(rcr$statistic$psi_hat, 24.05, tolerance = 1e-2)
  expect_equal(rcr$implied_size, 0.054, tolerance = 1e-3 / 0.054)
  expect_identical(rcr$verdict, "POINT_PASS")

  # THE naive-pilot danger (Prop. prop-pilot(i)) on real data
  rn <- eiv_adequacy(s$y, s$x, s$unit, s$time, sigma_nu = s$sd, pilot = "naive")
  expect_equal(rn$eta, 0.8853 * 0.5472, tolerance = 1e-2)
  expect_identical(rn$verdict, "POINT_PASS")  # the exact error the diagnostic prevents
  expect_true(any(grepl("ANTI-CONSERVATIVE", rn$notes)))

  s <- vdem_spec(v, "v2x_jucon", "v2x_jucon_sd")
  r <- eiv_adequacy(s$y, s$x, s$unit, s$time, sigma_nu = s$sd, pilot = "point")
  expect_equal(r$design$n, 8889)
  expect_equal(r$statistic$lambda_hat, 0.4125, tolerance = 1e-3 / 0.4125)
  expect_equal(r$eta, 12.10, tolerance = 1e-2)
  expect_equal(r$implied_size, 1, tolerance = 1e-6)
  expect_equal(r$breakdown, 0.929, tolerance = 2e-3 / 0.929)  # fixed point (paper Table 3)
  expect_identical(r$verdict, "FLAGGED")
  expect_true(any(grepl("quadratic", r$notes)))
  # flag SURVIVES clustering: psi_hat = 25.2 but eta_CR = 2.4, size 67%
  rcr <- eiv_adequacy(s$y, s$x, s$unit, s$time, sigma_nu = s$sd,
                      pilot = "point", cluster = "crve")
  expect_equal(rcr$statistic$psi_hat, 25.21, tolerance = 1e-2)
  expect_equal(rcr$eta, 2.41, tolerance = 1e-2)
  expect_equal(rcr$implied_size, 0.674, tolerance = 3e-3 / 0.674)
  expect_identical(rcr$verdict, "FLAGGED")
})

test_that("reference case 2b: gate-1 headline (spec 7.2)", {
  v <- read_gate1()
  r <- eiv_adequacy(v$ly, v$poly, v$iso, v$year, sigma_nu = v$poly_sd,
                    pilot = "point")
  expect_equal(r$statistic$lambda_hat, 0.868, tolerance = 1.5e-3 / 0.868)
  expect_equal(r$statistic$noise_ratio, 0.153, tolerance = 3e-3 / 0.153)
  expect_equal(r$statistic$beta_star, 0.072, tolerance = 1e-3 / 0.072)
  expect_equal(r$statistic$beta_corr, 0.083, tolerance = 1e-3 / 0.083)
})

test_that("repeated-report twins application", {
  tw <- get_dataset("twins")
  needed <- c("DLHRWAGE", "DEDUC1", "DEDUC2", "DTEN", "DMARRIED", "DUNCOV")
  d <- tw[stats::complete.cases(tw[needed]), needed]
  fit <- stats::lm(DLHRWAGE ~ DEDUC1 + DTEN + DMARRIED + DUNCOV, data = d)
  bstar <- unname(stats::coef(fit)[["DEDUC1"]])
  se_star <- unname(stats::coef(summary(fit))["DEDUC1", "Std. Error"])
  W <- stats::model.matrix(~ DTEN + DMARRIED + DUNCOV, data = d)
  x <- qr.resid(qr(W), d$DEDUC1)
  z <- qr.resid(qr(W), d$DEDUC2)
  y <- qr.resid(qr(W), d$DLHRWAGE)
  tau2 <- sum(x^2)
  sigma <- se_star * sqrt(tau2)
  lambda_cov <- reliability_from_repeats(x, z)
  lambda_equal <- reliability_from_repeats(x, z, method = "equal_variance")

  expect_equal(nrow(d), 147L)
  expect_equal(bstar, 0.0908758941551, tolerance = 1e-11)
  expect_equal(se_star, 0.0219814978882, tolerance = 1e-11)
  expect_equal(lambda_cov, 0.5749036782323, tolerance = 1e-11)
  expect_equal(lambda_equal, 0.5514169767299, tolerance = 1e-11)

  rcov <- eiv_adequacy_summary(bstar, sigma, tau2, 147, 4,
                               reliability = lambda_cov, pilot = "point")
  req <- eiv_adequacy_summary(bstar, sigma, tau2, 147, 4,
                              reliability = lambda_equal, pilot = "point")
  expect_equal(rcov$breakdown, 0.863710397039, tolerance = 1e-7)
  expect_equal(req$breakdown, rcov$breakdown, tolerance = 1e-13)
  expect_equal(rcov$statistic$beta_corr, 0.1580715128394, tolerance = 1e-11)
  expect_equal(req$statistic$beta_corr, 0.1648043096058, tolerance = 1e-11)
  # This coefficient agreement is algebraic, not external validation.
  expect_equal(rcov$statistic$beta_corr, sum(x * y) / sum(x * z), tolerance = 1e-11)
  expect_equal(rcov$eta, 3.056917186720, tolerance = 1e-11)
  expect_equal(req$eta, 3.363210998044, tolerance = 1e-11)
  expect_equal(rcov$implied_size, 0.863669337612, tolerance = 1e-11)
  expect_equal(req$implied_size, 0.919728454662, tolerance = 1e-11)
  expect_identical(rcov$verdict, "FLAGGED")
  expect_identical(req$verdict, "FLAGGED")

  rouse <- eiv_adequacy_summary(0.071, 1, 1 / 0.016^2, 445, 4,
                                reliability = 0.748, pilot = "point")
  expect_equal(rouse$breakdown, 0.871831787842, tolerance = 1e-7)
  expect_equal(rouse$statistic$beta_corr, 0.0949197860963, tolerance = 1e-11)
  expect_equal(rouse$eta, 1.494986631016, tolerance = 1e-11)
  expect_equal(rouse$implied_size, 0.321249033980, tolerance = 1e-11)
  expect_identical(rouse$verdict, "FLAGGED")

  out <- paste(utils::capture.output(print(rcov)), collapse = "\n")
  expect_match(out, "n=147")
  expect_match(out, "d_K=4")
  expect_false(grepl("N=0", out))
})

test_that("input validation and edge cases", {
  n0 <- 40
  uid <- rep(1:10, each = 4); tid <- rep(1:4, times = 10)
  k <- seq_len(n0)
  x <- sin(0.8 * k) + 0.1 * uid
  y <- x + 0.2 * cos(1.9 * k)

  expect_error(eiv_adequacy(y, x, uid, tid), "exactly one noise input")
  expect_error(eiv_adequacy(y, x, uid, tid, sigma_nu = 0.1, reliability = 0.9))
  expect_error(eiv_adequacy(y, x, uid, tid, reliability = 1.2))
  expect_error(eiv_adequacy(y, x, uid, tid, codelow = x), "both codelow")
  expect_error(eiv_adequacy(y, x, uid, tid, sigma_nu = c(0.1, 0.2)))
  expect_error(eiv_adequacy(y, x, uid, tid, sigma_nu = -0.1))
  expect_error(eiv_adequacy(y, x, uid, tid, reliability = 0.9,
                            reliability_lower = 0))
  expect_error(eiv_adequacy(y, x, uid, tid, reliability = 0.9,
                            gamma_lambda = 0.96))

  rh <- eiv_adequacy(y, x, uid, tid, sigma_nu = 100)
  expect_identical(rh$verdict, "FLAGGED")
  expect_equal(rh$implied_size, 1)
  expect_true(any(grepl("exceeds", rh$notes)))

  xb <- as.numeric(uid > 5 & tid >= 3)
  rb <- eiv_adequacy(y, xb, uid, tid, reliability = 0.9, pilot = "point")
  expect_true(any(grepl("MISCLASSIFICATION", rb$notes)))

  r1 <- eiv_adequacy(y, x, uid, tid, reliability = 1, pilot = "point")
  expect_equal(r1$eta, 0)
  expect_identical(r1$verdict, "POINT_PASS")
})
