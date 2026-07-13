# Module A (Paper A): exact identities vs dense computation + Table 4 sizes.

test_that("implied sizes (Thm 3.1 / Cor 3.1; Table 4 predictions)", {
  expect_equal(PD$.size_hc0(0.052, 0.05), 0.056, tolerance = 1e-3 / 0.056)
  expect_equal(PD$.size_hc0(0.104, 0.05), 0.064, tolerance = 1e-3 / 0.064)
  expect_equal(PD$.size_hc0(0.25, 0.05), 0.090, tolerance = 1e-3 / 0.090)
  expect_equal(PD$.size_hc0(0.50, 0.05), 0.166, tolerance = 1e-3 / 0.166)
  expect_lt(PD$.size_hc3(0.50, 0.05), 0.01)
  expect_equal(PD$.size_hc3(0, 0.05), 0.05, tolerance = 1e-10)
  expect_equal(PD$.rho_dagger(0.05, 0.05), 0.2957, tolerance = 1e-3 / 0.2957)
  expect_equal(PD$.size_hc0(PD$.rho_dagger(0.05, 0.05), 0.05), 0.10,
               tolerance = 1e-9)
})

test_that("FE leverage diagonal", {
  N <- 6; T <- 5
  uid <- rep(1:N, each = T); tid <- rep(1:T, times = N)
  p <- fe_leverage(uid, tid)
  expect_true(all(abs(p - (1 / T + 1 / N - 1 / (N * T))) < 1e-10))
  expect_equal(sum(p), N + T - 1, tolerance = 1e-8)

  uidu <- integer(0); tidu <- integer(0)
  for (i in 1:7) for (t in 1:5) {
    if ((2 * i + t) %% 6 == 0) next
    uidu <- c(uidu, i); tidu <- c(tidu, t)
  }
  D <- dense_dummies(uidu, tidu, 7, 5)
  P <- D %*% pinv_ref(crossprod(D)) %*% t(D)
  pu <- fe_leverage(uidu, tidu)
  expect_lt(max(abs(pu - diag(P))), 1e-8)
})

test_that("exact HC identities under uniform full leverage", {
  N <- 6; T <- 4
  s <- c(1, 1, 1, -1, -1, -1); r <- c(1, 1, -1, -1); cc <- 0.7
  uid <- rep(1:N, each = T); tid <- rep(1:T, times = N)
  n <- N * T
  k <- seq_len(n)
  x <- cc * s[uid] * r[tid] + 0.3 * uid + 0.5 * tid
  y <- 1.5 * x + sin(3.7 * k)

  rep_ <- leverage_report(y, x, uid, tid)
  st <- rep_$statistic
  Hbar <- 1 / T + 1 / N
  expect_equal(st$max_leverage, Hbar, tolerance = 1e-8)
  expect_equal(st$leverage_spread, 1, tolerance = 1e-8)
  expect_equal(st$se_hc1, st$se_hc2, tolerance = 1e-8)
  expect_equal(st$se_hc0 / st$se_hc2, sqrt(1 - Hbar), tolerance = 1e-8)
  expect_equal(st$se_hc3 / st$se_hc2, 1 / sqrt(1 - Hbar), tolerance = 1e-8)
  expect_equal(st$se_cjn, st$se_naive * sqrt(n / (n - (N + T - 1) - 1)),
               tolerance = 1e-12)
})

test_that("matches dense full-regression computation (unbalanced)", {
  uid <- integer(0); tid <- integer(0)
  for (i in 1:9) for (t in 1:6) {
    if ((i + 3 * t) %% 8 == 0) next
    uid <- c(uid, i); tid <- c(tid, t)
  }
  n <- length(uid); N <- 9; T <- 6
  k <- seq_len(n)
  x <- sin(1.3 * k) + 0.2 * uid - 0.1 * tid
  y <- 0.8 * x + cos(2.1 * k)

  D <- dense_dummies(uid, tid, N, T)
  A <- cbind(D, x)
  Ainv <- pinv_ref(crossprod(A))
  coefs <- Ainv %*% crossprod(A, y)
  H_dense <- diag(A %*% Ainv %*% t(A))
  u_dense <- as.numeric(y - A %*% coefs)
  d_K <- N + T - 1
  dof <- n - d_K - 1
  xt_dense <- as.numeric((diag(n) - D %*% pinv_ref(crossprod(D)) %*% t(D)) %*% x)
  tau2 <- sum(xt_dense^2)

  rep_ <- leverage_report(y, x, uid, tid)
  st <- rep_$statistic
  expect_equal(st$beta, coefs[length(coefs)], tolerance = 1e-8)
  expect_equal(st$max_leverage, max(H_dense), tolerance = 1e-8)
  expect_equal(st$se_cjn, sqrt(sum(u_dense^2) / dof / tau2), tolerance = 1e-8)
  expect_equal(st$se_hc0, sqrt(sum(xt_dense^2 * u_dense^2)) / tau2, tolerance = 1e-8)
  expect_equal(st$se_hc2, sqrt(sum(xt_dense^2 * u_dense^2 / (1 - H_dense))) / tau2,
               tolerance = 1e-8)
  expect_equal(st$se_hc3, sqrt(sum(xt_dense^2 * u_dense^2 / (1 - H_dense)^2)) / tau2,
               tolerance = 1e-8)
  expect_equal(rep_$design$d_K, d_K)
})

test_that("verdict logic and rendering", {
  N <- 50; T <- 20
  uid <- rep(1:N, each = T); tid <- rep(1:T, times = N)
  n <- N * T; k <- seq_len(n)
  x <- sin(0.9 * k) + 0.1 * uid
  y <- 2 * x + 0.3 * cos(1.7 * k)
  rep_ <- leverage_report(y, x, uid, tid)
  expect_lt(rep_$design$rho, 0.1)
  expect_identical(rep_$verdict, "CERTIFIED")
  expect_false(any(grepl("HC3 over-correcting", rep_$notes)))
  expect_equal(rep_$implied_size, PD$.size_hc0(rep_$design$rho, 0.05),
               tolerance = 1e-12)

  N2 <- 8; T2 <- 4
  uid2 <- rep(1:N2, each = T2); tid2 <- rep(1:T2, times = N2)
  n2 <- N2 * T2; k2 <- seq_len(n2)
  x2 <- cos(1.1 * k2) + 0.2 * tid2
  y2 <- x2 + 0.5 * sin(2.3 * k2)
  rep2 <- leverage_report(y2, x2, uid2, tid2)
  expect_gt(rep2$design$rho, PD$.rho_dagger(0.05, 0.05))
  expect_identical(rep2$verdict, "FLAGGED")
  expect_gt(rep2$implied_size, 0.10)
  expect_true(any(grepl("HC3 over-correcting", rep2$notes)))
  expect_equal(rep2$breakdown, PD$.rho_dagger(0.05, 0.05), tolerance = 1e-12)

  out <- paste(utils::capture.output(print(rep2)), collapse = "\n")
  expect_match(out, "Leverage / Variance \\(Paper A\\)")
  expect_match(out, "SE\\(beta\\)")
  expect_match(out, "VERDICT: FLAGGED")

  expect_error(leverage_report(y2[1:5], x2[1:5], uid2[1:5], tid2[1:5]))
  expect_error(leverage_report(y2, as.numeric(uid2), uid2, tid2),
               "no within variation")
})
