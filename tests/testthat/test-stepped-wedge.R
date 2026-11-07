# Module D (stepped-wedge design; Halkiewicz and Alexeev, arXiv:2610.35853).
# Reference values: tests/testthat/fixtures/sw_reference.csv, a verbatim copy of
# the frozen Python output 2026-11-07_sw_package_reference.csv (swdesign.py full
# GLS/TWFE, Theorem 2 and Theorem 3 closed forms, Corollary 1, Haines cells).

ref <- utils::read.csv(testthat::test_path("fixtures", "sw_reference.csv"),
                       stringsAsFactors = FALSE)
num <- function(x) if (is.na(x) || identical(x, "")) NULL else as.numeric(x)
des <- function(s) as.numeric(strsplit(s, ";", fixed = TRUE)[[1]])
cov_args <- function(r) {
  a <- list(estimator = r$estimator)
  if (r$estimator == "GLS") {
    if (!is.null(num(r$phi))) a$phi <- num(r$phi)
    if (!is.null(num(r$icc))) { a$icc <- num(r$icc); a$cac <- num(r$cac); a$m <- num(r$m) }
  }
  a
}
TOL <- 1e-10

test_that("Theorem 2 closed form: Gamma^2 along t(q, phi), T = 3..10", {
  t_of <- function(T, q, phi) if (is.infinite(phi)) (2 * T - 1 - q * T) / ((T - 1) * (T + 1)) else
    ((1 - q) + (2 * T - 1 - q * T) * phi) / ((T - 1) * ((1 + q) + (T + 1) * phi))
  G2 <- function(T, t) (T - 2) * ((T - 3) * (T + 1) * t^2 + 15) / (5 * T * (3 - (T + 1) * t)^2)
  for (T in 3:10) for (q in c(0, 0.1, 0.25, 0.5, 0.8)) for (phi in c(0, 0.5, 2, 10, Inf)) {
    d <- c(rep((1 - q) / (T - 1), T - 1), q)
    est <- if (is.infinite(phi)) "TWFE" else "GLS"
    g <- if (est == "TWFE") sw_gamma(design = d, estimator = "TWFE") else sw_gamma(design = d, phi = phi)
    expect_equal(g^2, G2(T, t_of(T, q, phi)), tolerance = TOL)
    # weights lie on the one-parameter curve W_e ~ (T-e)(1-te)
    W <- if (est == "TWFE") sw_exposure_weights(d, estimator = "TWFE") else sw_exposure_weights(d, phi = phi)
    t <- t_of(T, q, phi); e <- seq_len(T - 1)
    expect_equal(unname(W), 6 * (T - e) * (1 - t * e) / (T * (T - 1) * (3 - (T + 1) * t)), tolerance = TOL)
  }
})

test_that("Theorem 1 (TWFE complete wedge) closed form, T = 3..10", {
  for (T in 3:10) {
    G2 <- (19 * T^3 - 31 * T^2 - 2 * T + 12) / (5 * T * (T - 2) * (T + 1))
    expect_equal(sw_gamma(design = c(rep(1, T - 1), 0), estimator = "TWFE")^2, G2, tolerance = TOL)
  }
})

test_that("Theorem 3: Var_ETI(q) = u'(B - sC)^{-1}u / s; ETI cost is unimodal in q", {
  for (T in c(4, 5, 7)) for (phi in c(0.5, 2, Inf)) {
    est <- if (is.infinite(phi)) "TWFE" else "GLS"
    f <- function(q) {
      d <- c(rep((1 - q) / (T - 1), T - 1), q)
      if (est == "TWFE") sw_eti_cost(d, estimator = "TWFE")$eti_cost else sw_eti_cost(d, phi = phi)$eti_cost
    }
    v <- vapply(seq(0, 0.9, by = 0.01), f, 0)
    s <- sign(diff(v)); s <- s[s != 0]
    expect_true(s[1] < 0)                       # falls at q = 0
    expect_lte(sum(s[-1] != s[-length(s)]), 1)  # at most one turn
  }
})

test_that("reference values: balanced, integer and nested designs (Python, 1e-10)", {
  rows <- ref[ref$kind %in% c("balanced", "integer", "nested") & ref$quantity != "phi_prime", ]
  cases <- split(rows, rows$case)
  for (cs in cases) {
    r <- cs[1, ]; d <- des(r$design); a <- cov_args(r)
    if (r$kind == "nested") a <- list(estimator = "GLS", icc = num(r$icc), cac = num(r$cac), m = num(r$m))
    W <- do.call(sw_exposure_weights, c(list(design = d), a))
    ec <- do.call(sw_eti_cost, c(list(design = d), a))
    for (i in seq_len(nrow(cs))) {
      q <- cs$quantity[i]; val <- cs$value[i]
      got <- switch(q,
        Gamma = sw_gamma(W), kappa = sw_radius(W), it_ratio = ec$it_ratio,
        eti_cost = ec$eti_cost, var_it = ec$var_it, var_eti = ec$var_eti,
        unname(W[as.integer(sub("W_", "", q))]))
      expect_equal(got, val, tolerance = TOL, label = paste(r$case, q))
    }
  }
  pp <- ref[ref$quantity == "phi_prime", ]
  for (i in seq_len(nrow(pp)))
    expect_equal(sw_phi_nested(num(pp$icc[i]), num(pp$cac[i]), num(pp$m[i])), pp$value[i], tolerance = TOL)
})

test_that("paper numbers: kappa* 0.133 and ETI cost 2.27 (T = 5, phi = 2); capped search frontier", {
  expect_equal(round(sw_radius(design = c(6, 6, 6, 6, 0), phi = 2), 3), 0.133)
  expect_equal(round(sw_eti_cost(c(6, 6, 6, 6, 0), phi = 2)$eti_cost, 2), 2.27)
  s <- sw_design_search(G = 24, T = 5, cap = 6, max_wait = 12, phi = 2)
  fr <- s$frontier[match(c(4, 6, 8, 12), s$frontier$max_wait), ]
  expect_equal(round(fr$eti_cost, 2), c(1.52, 1.38, 1.29, 1.23))
  for (k in c(0, 4, 6, 8, 12)) {
    r <- ref[ref$case == sprintf("search_eti_phi2.0_cap6_wl%d", k), ]
    f <- s$frontier[s$frontier$max_wait == k, ]
    expect_equal(f$eti_cost, r$value, tolerance = TOL)
    expect_equal(as.numeric(f[1, c("start2", "start3", "start4", "start5", "wait")]), des(r$design))
  }
  rr <- ref[ref$case == "search_radius_phi2_complete", ]
  b <- sw_design_search(24, 5, max_wait = 0, phi = 2, objective = "radius", max_it_ratio = 1)$best
  expect_equal(b$radius, rr$value, tolerance = TOL)
  expect_equal(as.numeric(b[1, 1:5]), des(rr$design))
})

test_that("Haines actual cells reproduce kappa* = 0.1244 from the bundled data", {
  h <- get_dataset("haines2017")
  n <- with(h, tapply(rep(1, nrow(h)), list(ward, period), sum))
  starts <- tapply(h$start, h$ward, `[`, 1)
  ord <- order(starts, as.numeric(names(starts)))
  r <- ref[ref$case == "haines_actual", ]
  icc <- num(r$icc[1])
  W <- sw_exposure_weights(c(2, 2, 2, 2, 2, 2, 0), icc = icc, cac = 1, sizes = n[ord, ])
  for (e in 1:6) expect_equal(unname(W[e]), r$value[r$quantity == paste0("W_", e)], tolerance = TOL)
  k <- sw_radius(W)
  expect_equal(k, r$value[r$quantity == "kappa"], tolerance = TOL)
  expect_equal(round(k, 4), 0.1244)
})

test_that("unequal sizes reduce to equal sizes; nested reduces to phi'", {
  sz <- matrix(20, 24, 5)
  d <- c(5, 5, 5, 5, 4)
  expect_equal(sw_exposure_weights(d, icc = 0.05, cac = 0.5, sizes = sz),
               sw_exposure_weights(d, icc = 0.05, cac = 0.5, m = 20), tolerance = TOL)
  expect_equal(sw_exposure_weights(d, icc = 0.05, cac = 0.5, m = 20),
               sw_exposure_weights(d, phi = sw_phi_nested(0.05, 0.5, 20)), tolerance = TOL)
  expect_equal(sw_eti_cost(d, estimator = "TWFE", sizes = sz)$eti_cost,
               sw_eti_cost(d, estimator = "TWFE")$eti_cost, tolerance = TOL)
  expect_equal(sw_phi_nested(0.1, 1, 30), 30 * 0.1 / 0.9, tolerance = TOL)
})

test_that("input validation", {
  expect_error(sw_exposure_weights(c(0, 6, 6, 6, 6), phi = 2), "start period 2")
  expect_error(sw_exposure_weights(c(6, 6), phi = 2), "length T >= 3")
  expect_error(sw_exposure_weights(c(6, 6, 6, 6, 0), phi = 2, icc = 0.1, m = 10), "not both")
  expect_error(sw_exposure_weights(c(6, 6, 6, 6, 0)), "GLS needs")
  expect_error(sw_eti_cost(c(24, 0, 0, 0, 0), phi = 2), "does not identify")
  expect_error(sw_gamma(c(0.5, 0.5), target = c(1, 1)), "sum to one")
})
