# Paper A: exact inference under concentrated identifying variation.
# Mirrors julia/test/cycle_tests.jl. Monte Carlo replication counts are smaller
# than the Julia suite's because R is slower; tolerances are widened to match,
# and every randomized block is seeded.

test_that("kappa = 1 on designs that are a single cycle", {
  # one four-cycle: the contrast spans the whole cycle space
  x <- c(1, 3, 0, 2); u <- c(1, 2, 2, 1); tt <- c(1, 1, 2, 2)
  for (m in c("greedy", "structured", "sparse"))
    expect_equal(cycle_capture(x, u, tt, method = m), 1, tolerance = 1e-10)
  # one digon (a repeated cell)
  expect_equal(cycle_capture(c(2, 5), c(1, 1), c(1, 1), method = "greedy"), 1,
               tolerance = 1e-10)
})

test_that("contrasts annihilate BOTH fixed effects exactly", {
  set.seed(11)
  N <- 8; T <- 6
  u <- rep(1:N, each = T); tt <- rep(1:T, times = N)
  x <- stats::rnorm(N * T)
  for (m in c("greedy", "structured", "sparse")) {
    cs <- cycle_contrasts(x, u, tt, method = m)
    expect_gt(length(cs$rows), 0)
    for (c in seq_along(cs$rows)) {
      du <- tapply(cs$signs[[c]], u[cs$rows[[c]]], sum)
      dt <- tapply(cs$signs[[c]], tt[cs$rows[[c]]], sum)
      # exact, not approximate: this is what makes the test exact
      expect_true(all(du == 0))
      expect_true(all(dt == 0))
    }
    # supports are edge-disjoint (needed for sign-flip exchangeability)
    allr <- unlist(cs$rows)
    expect_identical(length(unique(allr)), length(allr))
    expect_gt(cs$kappa, 0)
    expect_lte(cs$kappa, 1 + 1e-12)
  }
})

test_that("loadings are exact: v'x equals v'xtilde", {
  set.seed(3)
  N <- 6; T <- 5
  u <- rep(1:N, each = T); tt <- rep(1:T, times = N)
  x <- stats::rnorm(N * T)
  cs <- cycle_contrasts(x, u, tt)
  xt <- twoway_demean(x, u, tt)
  for (c in seq_along(cs$rows)) {
    L <- length(cs$rows[[c]])
    b_raw <- sum(cs$signs[[c]] * x[cs$rows[[c]]]) / sqrt(L)
    b_res <- sum(cs$signs[[c]] * xt[cs$rows[[c]]]) / sqrt(L)
    expect_equal(b_raw, b_res, tolerance = 1e-9)
    expect_equal(cs$loadings[c], b_raw, tolerance = 1e-12)
  }
})

test_that("packing is deterministic (the cross-language parity requirement)", {
  set.seed(5)
  N <- 7; T <- 7
  u <- rep(1:N, each = T); tt <- rep(1:T, times = N)
  x <- as.numeric(sample(1:9, N * T, replace = TRUE))   # discrete: many exact ties
  expect_identical(cycle_capture(x, u, tt), cycle_capture(x, u, tt))
  expect_identical(cycle_capture(x, u, tt, method = "greedy"),
                   cycle_capture(x, u, tt, method = "greedy"))
})

test_that("structured packing dominates naive greedy", {
  set.seed(9)
  N <- 10; T <- 8
  u <- rep(1:N, each = T); tt <- rep(1:T, times = N)
  x <- stats::rnorm(N * T)
  expect_gt(cycle_capture(x, u, tt, method = "structured"),
            cycle_capture(x, u, tt, method = "greedy"))
})

test_that("sparse packing handles the mobility-network geometry", {
  # two-edge movers between firm pairs: the AKM geometry, where the dense
  # routine's O(N^2 T^2) enumeration is infeasible
  set.seed(21)
  nw <- 400; nf <- 20
  a <- rep(seq_len(nw), each = 2)
  b <- as.integer(sample(seq_len(nf), 2 * nw, replace = TRUE))
  x <- stats::rnorm(length(a))
  cs <- cycle_contrasts(x, a, b, method = "sparse")
  expect_gt(length(cs$rows), 0)
  expect_lte(cs$kappa, 1 + 1e-12)
  for (c in seq_along(cs$rows)) {
    expect_equal(sum(tapply(cs$signs[[c]], a[cs$rows[[c]]], sum)^2), 0)
    expect_equal(sum(tapply(cs$signs[[c]], b[cs$rows[[c]]], sum)^2), 0)
  }
  allr <- unlist(cs$rows)
  expect_identical(length(unique(allr)), length(allr))
  # extreme (nested) pairing beats naive greedy on this geometry
  expect_gt(cs$kappa, cycle_capture(x, a, b, method = "greedy"))
  expect_identical(cs$kappa, cycle_contrasts(x, a, b, method = "sparse")$kappa)
  expect_true(support_compatibility(cs, a)$compatible)
  expect_gt(signflip_test(stats::rnorm(length(x)), cs, nflips = 49,
                          blocks = a)$effective_C, 0)
})

test_that("custom multiway contrasts enforce FE and block compatibility", {
  worker <- c(1, 1, 2, 2)
  firm <- c(10, 10, 20, 20)
  year <- c(2000, 2001, 2000, 2001)
  x <- c(0, 0, 0, 1)
  cs <- contrast_system(x, list(worker, firm, year), list(1:4),
                        list(c(-1, 1, 1, -1)), blocks = worker)
  expect_equal(cs$kappa, 1, tolerance = 1e-12)
  expect_true(support_compatibility(cs, worker)$compatible)
  expect_lt(max(abs(multiway_demean(x, list(worker, firm, year)) -
                    twoway_demean(x, worker, year))), 1e-10)

  u <- rep(1:5, each = 4); tt <- rep(1:4, times = 5)
  xx <- sin(0.7 * seq_along(u)); yy <- cos(1.1 * xx)
  cs2 <- cycle_contrasts(xx, u, tt)
  expect_false(support_compatibility(cs2, u)$compatible)
  expect_error(signflip_test(yy, cs2, nflips = 49, blocks = u),
               "not unions")
})

test_that("cluster diagnostics aggregate at the dependence-block level", {
  set.seed(2026)
  N <- 12; worker <- rep(seq_len(N), each = 2); year <- rep(1:2, times = N)
  x <- stats::rnorm(2 * N); y <- 0.4 * x + stats::rnorm(2 * N)
  rep_ <- cycle_report(y, x, worker, year, blocks = worker,
                       interval = FALSE, nflips = 49)

  xt <- twoway_demean(x, worker, year)
  block_mass <- as.numeric(rowsum(xt^2, worker, reorder = FALSE))
  expect_identical(rep_$statistic$concentration_level, "block")
  expect_equal(rep_$statistic$concentration_blocks, N)
  expect_equal(rep_$statistic$lambda_n, max(block_mass) / sum(block_mass))

  yt <- twoway_demean(y, worker, year)
  u <- yt - rep_$statistic$beta_ols * xt
  cluster_score <- as.numeric(rowsum(xt * u, worker, reorder = FALSE))
  score_share <- cluster_score^2 / sum(cluster_score^2)
  expect_equal(rep_$statistic$score_lambda_n, max(score_share))
  expect_equal(rep_$statistic$score_H_n, sum(score_share^2))
})

test_that("only treatment-loaded supports determine orbit granularity", {
  cs <- cycle_contrasts(c(1, 1, 0, 2), c(1, 1, 2, 2), c(1, 1, 2, 2),
                        method = "greedy")
  out <- signflip_test(c(0.2, -0.1, 0.4, -0.3), cs, nflips = 99)
  expect_identical(out$C, 2L)
  expect_identical(out$effective_C, 1L)
  expect_equal(out$full_enumeration_floor, 1)
})

test_that("real Spec-G panel: V_n and lambda_n match the article", {
  d <- get_dataset("fscore")
  cy <- paste(d$country, d$year)
  cs <- cycle_contrasts(d$fscore, d$uid, cy)
  # V_n is a property of the design, not of the packing: exact match
  expect_equal(cs$V_n, 461.693976, tolerance = 1e-5 / 461.693976)
  xt <- twoway_demean(d$fscore, d$uid, cy)
  expect_equal(max(xt^2) / cs$V_n, 0.036, tolerance = 5e-4 / 0.036)
  # kappa is a heuristic lower bound and tie-break dependent; the article's run
  # reports 0.795. Ours is deterministic and must be at least as good.
  expect_gte(cs$kappa, 0.79)
  expect_lt(cycle_capture(d$fscore, d$uid, cy, method = "greedy"), cs$kappa)
})

test_that("sign-flip test is EXACT under a non-Gaussian concentrated design", {
  # The property the module exists for: the conventional t-test fails badly
  # here, the exact test does not.
  set.seed(7)
  N <- 10; T <- 4
  u <- rep(1:N, each = T); tt <- rep(1:T, times = N)
  n <- N * T
  x <- stats::rnorm(n); x[1] <- x[1] * 12          # concentrate the variation
  reps <- 300; nflips <- 199
  rej_c <- 0; rej_t <- 0
  for (i in seq_len(reps)) {
    eps <- sample(c(-1, 1), n, TRUE) * stats::rnorm(n)^2   # symmetric, heavy-tailed
    eps[1] <- eps[1] * 8
    if (signflip_test(eps, x, u, tt, beta0 = 0, nflips = nflips)$p <= 0.05)
      rej_c <- rej_c + 1
    if (abs(leverage_report(eps, x, u, tt)$statistic$t_hc2) > 1.96)
      rej_t <- rej_t + 1
  }
  expect_lt(rej_c / reps, 0.11)     # exact: holds its nominal level
  expect_gt(rej_t / reps, 0.20)     # conventional: badly over-rejects
})

test_that("confidence set inverts the test, and is affine-exact", {
  set.seed(2)
  N <- 8; T <- 5
  u <- rep(1:N, each = T); tt <- rep(1:T, times = N)
  n <- N * T
  x <- stats::rnorm(n)
  y <- 0.7 * x + stats::rnorm(n)
  set.seed(4)
  itv <- signflip_interval(y, x, u, tt, alpha = 0.05, ngrid = 2001, nflips = 999)
  expect_lt(itv$lo, itv$beta_tilde)          # centred on the contrast estimate
  expect_gt(itv$hi, itv$beta_tilde)
  # the p-value at beta_tilde is maximal (the statistic vanishes there)
  set.seed(4)
  expect_gt(signflip_test(y, x, u, tt, beta0 = itv$beta_tilde, nflips = 999)$p, 0.05)
  # a value far outside the set is rejected
  set.seed(4)
  far <- itv$hi + 10 * (itv$hi - itv$lo)
  expect_lte(signflip_test(y, x, u, tt, beta0 = far, nflips = 999)$p, 0.05)
})

test_that("report, verdicts and the 2^(1-C) floor", {
  set.seed(13)
  N <- 25; T <- 10                            # big enough that lambda_n is small
  u <- rep(1:N, each = T); tt <- rep(1:T, times = N)
  n <- N * T
  x <- stats::rnorm(n); y <- stats::rnorm(n)
  set.seed(1)
  rep_ <- cycle_report(y, x, u, tt, nflips = 499)
  expect_lt(rep_$statistic$lambda_n, 0.10)    # genuinely diffuse
  expect_identical(rep_$verdict, "POINT_PASS")
  expect_equal(rep_$statistic$min_pvalue,
               2^(1 - rep_$statistic$effective_C))
  expect_equal(rep_$statistic$se_price, 1 / sqrt(rep_$statistic$kappa))
  expect_true(is.finite(rep_$statistic$score_lambda_n))
  expect_gt(rep_$statistic$score_n_eff, 0)
  expect_true(grepl("Concentrated Identifying Variation",
                    paste(utils::capture.output(print(rep_)), collapse = " ")))

  # concentrated design -> FLAGGED, and it says why
  set.seed(13)
  xc <- stats::rnorm(n); xc[1] <- xc[1] * 30
  set.seed(1)
  repc <- cycle_report(y, xc, u, tt, nflips = 499, interval = FALSE)
  expect_identical(repc$verdict, "FLAGGED")
  expect_true(any(grepl("CONCENTRATION WARNING", repc$notes)))
  concise <- paste(utils::capture.output(print(repc)), collapse = "\n")
  expect_match(concise, "Diagnostic notes hidden", fixed = TRUE)
  expect_false(grepl("CONCENTRATION WARNING", concise, fixed = TRUE))
  expect_false(grepl("\nNote:", concise, fixed = TRUE))
  detailed <- paste(utils::capture.output(print(repc, notes = TRUE)), collapse = "\n")
  expect_match(detailed, "Diagnostic notes for", fixed = TRUE)
  expect_match(detailed, "CONCENTRATION WARNING", fixed = TRUE)
  standalone <- paste(utils::capture.output(show_notes(repc)), collapse = "\n")
  expect_match(standalone, "1. ", fixed = TRUE)

  # too few supports for any level-alpha test to exist
  rep2 <- cycle_report(c(1, 3, 0, 2), c(1, 3, 0, 2), c(1, 2, 2, 1), c(1, 1, 2, 2),
                       nflips = 99, interval = FALSE)
  expect_identical(rep2$statistic$C, 1L)
  expect_identical(rep2$statistic$effective_C, 1L)
  expect_equal(rep2$statistic$min_pvalue, 1)
  expect_identical(rep2$verdict, "INCONCLUSIVE")   # 2^(1-1) = 1 > alpha
  expect_match(rep2$statistic$reason, "may repair", fixed = TRUE)
  concise2 <- paste(utils::capture.output(print(rep2)), collapse = "\n")
  expect_match(concise2, "Reason:", fixed = TRUE)
  expect_match(concise2, "may repair", fixed = TRUE)

  # Binary-treatment granularity is structural, holding the panel fixed.
  unit_b <- rep(1:6, each = 2)
  time_b <- rep(1L, 12)
  y_b <- sin(1:12) + 0.1 * cos(2 * (1:12))
  x_low <- rep(c(1, 0), 6)
  x_low[7:12] <- 0                         # n1 = 3
  x_high <- rep(c(1, 0), 6)                # n1 = 6

  pre_low <- applicable(x_low, unit_b, time_b)
  pre_high <- applicable(x_high, unit_b, time_b)
  expect_false(pre_low$ok)
  expect_match(pre_low$reason, "n1 = 3 treated", fixed = TRUE)
  expect_match(pre_low$reason, "Not repairable by repacking", fixed = TRUE)
  expect_true(pre_high$ok)

  low <- cycle_report(y_b, x_low, unit_b, time_b,
                      method = "greedy", interval = FALSE)
  high <- cycle_report(y_b, x_high, unit_b, time_b,
                       method = "greedy", interval = FALSE)
  expect_identical(low$verdict, "INCONCLUSIVE")
  expect_identical(low$statistic$n_treated, 3L)
  expect_false(low$statistic$binary_floor_ok)
  expect_match(low$statistic$reason, "floor 2^(1-3) = 0.25", fixed = TRUE)
  expect_false(identical(high$verdict, "INCONCLUSIVE"))
  expect_identical(high$statistic$effective_C, 6L)
  expect_true(high$statistic$binary_floor_ok)

  x_mostly <- 1 - x_low
  pre_mostly <- applicable(x_mostly, unit_b, time_b)
  expect_false(pre_mostly$ok)
  expect_match(pre_mostly$reason, "n1 = 9 treated and n0 = 3", fixed = TRUE)
  mostly <- cycle_report(y_b, x_mostly, unit_b, time_b,
                         method = "greedy", interval = FALSE)
  expect_identical(mostly$verdict, "INCONCLUSIVE")
  expect_identical(mostly$statistic$effective_C, 3L)

  row <- adequacy_row(x_high, unit_b, time_b, method = "greedy")
  expect_identical(row$n_treated, 6L)
  expect_true(row$binary_floor_ok)
  continuous <- adequacy_row(as.numeric(1:12), unit_b, time_b,
                             method = "greedy")
  expect_true(is.na(continuous$n_treated))
  expect_true(continuous$binary_floor_ok)
})
