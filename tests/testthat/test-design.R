# Shared infrastructure (spec section 2) + reference case 7.1.

test_that("union-find d_K", {
  d <- design_summary(rep(1:5, each = 4), rep(1:4, times = 5))
  expect_equal(c(d$N, d$T, d$n), c(5, 4, 20))
  expect_equal(d$ncomponents, 1)
  expect_equal(d$d_K, 5 + 4 - 1)
  expect_equal(d$rho, d$d_K / 20)

  # two disconnected blocks
  uid2 <- c(rep(1:3, each = 2), rep(4:6, each = 2))
  tid2 <- c(rep(1:2, times = 3), rep(3:4, times = 3))
  d2 <- design_summary(uid2, tid2)
  expect_equal(d2$ncomponents, 2)
  expect_equal(d2$d_K, 6 + 4 - 2)

  d3 <- design_summary(rep(1, 7), 1:7)
  expect_equal(d3$d_K, 1 + 7 - 1)

  d4 <- design_summary(c("a", "a", "b", "b"), c(2001, 2002, 2001, 2002))
  expect_equal(c(d4$N, d4$T, d4$d_K), c(2, 2, 3))
})

test_that("two-way demeaning matches the dense projection (unbalanced)", {
  uid <- integer(0); tid <- integer(0)
  for (i in 1:8) for (t in 1:6) {
    if ((i + 2 * t) %% 7 == 0) next
    uid <- c(uid, i); tid <- c(tid, t)
  }
  n <- length(uid)
  expect_lt(n, 48)
  x <- sin(0.7 * seq_len(n) + 20260710 %% 13) + 0.3 * seq_len(n)^2 / n

  xt <- twoway_demean(x, uid, tid)
  D <- dense_dummies(uid, tid, 8, 6)
  M <- diag(n) - D %*% pinv_ref(crossprod(D)) %*% t(D)
  expect_lt(max(abs(xt - as.numeric(M %*% x))), 1e-8)

  # idempotence and cluster-centering (Paper C, Lemma lem-center)
  expect_lt(max(abs(twoway_demean(xt, uid, tid) - xt)), 1e-8)
  for (i in 1:8) expect_lt(abs(sum(xt[uid == i])), 1e-8)
})

test_that("block-design within treatment is uniform on treated cells", {
  N <- 10; T <- 6; g <- 4
  uid <- rep(1:N, each = T); tid <- rep(1:T, times = N)
  D <- as.numeric(uid <= 4 & tid >= g)
  Dt <- twoway_demean(D, uid, tid)
  w <- Dt[D == 1]
  expect_lt(max(abs(w - w[1])), 1e-8)
})

test_that("reference case 1: V-Dem design primitive (spec 7.1)", {
  v <- read_gate1()
  d <- design_summary(v$iso, v$year, x = v$poly)
  expect_equal(d$n, 8704)
  expect_equal(d$N, 174)
  expect_equal(d$T, 60)
  expect_equal(d$d_K, 233)
  expect_equal(d$ncomponents, 1)
  expect_equal(d$rho, 0.0268, tolerance = 5e-4 / 0.0268)

  # gate-1 headline lambda_hat = 0.868 from pure infrastructure arithmetic
  a_hat <- mean(v$poly_sd^2) * (d$n - d$d_K)
  lambda_hat <- 1 - a_hat / d$tau_star2
  expect_equal(lambda_hat, 0.868, tolerance = 1e-3 / 0.868)
})

test_that("AdequacyReport construction and printing", {
  d <- design_summary(rep(1:5, each = 4), rep(1:4, times = 5))
  r <- PD$.new_AdequacyReport("measurement_error", d, list(lambda_hat = 0.868),
                              0.28, 0.652, 0.82, 0.053, "CERTIFIED", 0.05, 0.05,
                              "corrected, conservative pilot used")
  out <- paste(utils::capture.output(print(r)), collapse = "\n")
  expect_match(out, "Panel Adequacy Report — Measurement Error \\(Paper B\\)")
  expect_match(out, "lambda_hat = 0.868")
  expect_match(out, "\\|eta\\| = 0.280")
  expect_match(out, "Threshold \\(delta=0.05\\) = 0.652")
  expect_match(out, "Implied size of nominal 5% test: 5.3%")
  expect_match(out, "VERDICT: CERTIFIED at delta=0.05")
  expect_match(out, "Note: corrected, conservative pilot")

  expect_error(PD$.new_AdequacyReport("conformal", d, list(), NULL, NULL, NULL,
                                      NULL, "CERTIFIED", 0.05, 0.05, character(0)),
               "unknown pathology")
  expect_error(PD$.new_AdequacyReport("leverage", d, list(), NULL, NULL, NULL,
                                      NULL, "MAYBE", 0.05, 0.05, character(0)),
               "unknown verdict")
})
