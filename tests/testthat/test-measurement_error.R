# Module B (Paper B): threshold constants, V-Dem two-pole (spec 7.2), gate-1
# headline, PSID application, honesty machinery. Breakdown expectations are the
# exact-inversion values (the stale linear ones must NOT reappear).

vdem_spec <- function(v, xcol, sdcol) {
  keep <- stats::complete.cases(v$ly, v[[xcol]], v[[sdcol]])
  list(unit = v$iso[keep], time = v$year[keep], y = v$ly[keep],
       x = v[[xcol]][keep], sd = v[[sdcol]][keep])
}

test_that("threshold constants (Paper B, Remark rem-exact-cv)", {
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
})

test_that("reliability helpers", {
  expect_equal(reliability_from_interval(c(0.1, 0.2), c(0.3, 0.5)), c(0.1, 0.15))
  expect_equal(reliability_from_ratio(0.8, 2), sqrt(0.25) * 2)
  expect_error(reliability_from_ratio(1.2, 2))
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
  expect_identical(r$statistic$lambda_hat >= r$breakdown, r$verdict == "CERTIFIED")
  expect_identical(r$verdict, "CERTIFIED")
  rc <- eiv_adequacy(s$y, s$x, s$unit, s$time, sigma_nu = s$sd)
  expect_identical(rc$verdict, "CERTIFIED")
  expect_gt(rc$statistic$eta_upper, rc$eta)
  # cluster-robust (country CRVE): paper Table 3 psi_hat = 19.2, still certified
  rcr <- eiv_adequacy(s$y, s$x, s$unit, s$time, sigma_nu = s$sd,
                      pilot = "point", cluster = "crve")
  expect_equal(rcr$statistic$psi_hat, 19.18, tolerance = 1e-2)
  expect_equal(rcr$implied_size, 0.050, tolerance = 1e-3 / 0.050)
  expect_identical(rcr$verdict, "CERTIFIED")

  s <- vdem_spec(v, "v2xlg_legcon", "v2xlg_legcon_sd")
  r <- eiv_adequacy(s$y, s$x, s$unit, s$time, sigma_nu = s$sd, pilot = "point")
  expect_equal(r$design$n, 8529)
  expect_equal(r$statistic$lambda_hat, 0.5472, tolerance = 1e-3 / 0.5472)
  expect_equal(r$eta, 0.8853, tolerance = 5e-3)
  expect_equal(r$implied_size, 0.1435, tolerance = 1e-3 / 0.1435)
  expect_equal(r$breakdown, 0.621, tolerance = 2e-3 / 0.621)  # fixed point (paper Table 3)
  expect_identical(r$verdict, "FLAGGED")
  # the paper's middle case: flagged iid, CERTIFIED under country clustering
  rcr <- eiv_adequacy(s$y, s$x, s$unit, s$time, sigma_nu = s$sd,
                      pilot = "point", cluster = "crve")
  expect_equal(rcr$statistic$psi_hat, 24.05, tolerance = 1e-2)
  expect_equal(rcr$implied_size, 0.054, tolerance = 1e-3 / 0.054)
  expect_identical(rcr$verdict, "CERTIFIED")

  # THE naive-pilot danger (Prop. prop-pilot(i)) on real data
  rn <- eiv_adequacy(s$y, s$x, s$unit, s$time, sigma_nu = s$sd, pilot = "naive")
  expect_equal(rn$eta, 0.8853 * 0.5472, tolerance = 1e-2)
  expect_identical(rn$verdict, "CERTIFIED")   # the exact error Paper B prevents
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

test_that("PSID application via summary form (Paper B app-psid)", {
  bstar <- 0.7332110; sigma <- 4.2465084; tau2 <- 83.1775587
  n <- 4165; d_K <- 595 + 7 - 1

  expect_equal(breakdown_reliability(bstar, sigma, tau2), 0.707,
               tolerance = 2e-3 / 0.707)
  # cluster-robust breakdown at the paper's person-cluster psi = 2.33: 0.61
  expect_equal(breakdown_reliability(bstar, sigma, tau2, psi = 2.330), 0.613,
               tolerance = 2e-3 / 0.613)

  r <- eiv_adequacy_summary(bstar, sigma, tau2, n, d_K, reliability = 0.65,
                            pilot = "point")
  expect_equal(r$eta, 0.848, tolerance = 2e-3 / 0.848)
  expect_equal(r$implied_size, 0.1356, tolerance = 1e-3 / 0.1356)
  expect_identical(r$verdict, "FLAGGED")
  # ... but CERTIFIED under person clustering (paper Table 4: size 8.6%)
  rpsi <- eiv_adequacy_summary(bstar, sigma, tau2, n, d_K, reliability = 0.65,
                               pilot = "point", psi = 2.330)
  expect_equal(rpsi$eta, 0.556, tolerance = 2e-3 / 0.556)
  expect_equal(rpsi$implied_size, 0.086, tolerance = 1e-3 / 0.086)
  expect_identical(rpsi$verdict, "CERTIFIED")
  expect_equal(rpsi$breakdown, 0.613, tolerance = 2e-3 / 0.613)

  rp <- eiv_adequacy_summary(bstar, sigma, tau2, n, d_K, reliability = 0.82,
                             pilot = "point")
  expect_equal(rp$implied_size, 0.0638, tolerance = 1e-3 / 0.0638)
  expect_identical(rp$verdict, "CERTIFIED")
  rf <- eiv_adequacy_summary(bstar, sigma, tau2, n, d_K, reliability = 0.82)
  expect_identical(rf$verdict, "FLAGGED")   # formal certificate fails (paper ddagger)
  expect_equal(rf$statistic$eta_upper, 0.707, tolerance = 2e-3 / 0.707)

  out <- paste(utils::capture.output(print(rf)), collapse = "\n")
  expect_match(out, "n=4165")
  expect_match(out, "d_K=601")
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

  rh <- eiv_adequacy(y, x, uid, tid, sigma_nu = 100)
  expect_identical(rh$verdict, "FLAGGED")
  expect_equal(rh$implied_size, 1)
  expect_true(any(grepl("exceeds", rh$notes)))

  xb <- as.numeric(uid > 5 & tid >= 3)
  rb <- eiv_adequacy(y, xb, uid, tid, reliability = 0.9, pilot = "point")
  expect_true(any(grepl("MISCLASSIFICATION", rb$notes)))

  r1 <- eiv_adequacy(y, x, uid, tid, reliability = 1, pilot = "point")
  expect_equal(r1$eta, 0)
  expect_identical(r1$verdict, "CERTIFIED")
})
