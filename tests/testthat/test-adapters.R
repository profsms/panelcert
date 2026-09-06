# Model-ingestion adapters: the diagnostic run on a fitted model must equal
# the diagnostic run on the raw data the model was fitted to (spec section 0:
# consume the already-estimated model, never re-specify).

test_that("fixest adapter reproduces the raw-data diagnostic (Module B)", {
  skip_if_not_installed("fixest")
  v <- get_dataset("vdem")
  d <- v[stats::complete.cases(v$ly, v$v2x_polyarchy, v$v2x_polyarchy_sd), ]

  m <- fixest::feols(ly ~ v2x_polyarchy | iso + year, data = d)
  r_model <- eiv_adequacy(m, sigma_nu = d$v2x_polyarchy_sd, pilot = "point")
  r_raw <- eiv_adequacy(d$ly, d$v2x_polyarchy, d$iso, d$year,
                        sigma_nu = d$v2x_polyarchy_sd, pilot = "point")
  expect_equal(r_model$statistic$lambda_hat, r_raw$statistic$lambda_hat,
               tolerance = 1e-10)
  expect_equal(r_model$eta, r_raw$eta, tolerance = 1e-10)
  expect_identical(r_model$verdict, r_raw$verdict)
  expect_equal(r_model$design$d_K, r_raw$design$d_K)
})

test_that("fixest adapter reproduces the raw-data diagnostic (Module C)", {
  skip_if_not_installed("fixest")
  castle <- read_panel("castle_panel.csv")
  castle$D <- as.numeric(!is.na(castle$ft) & castle$tid >= castle$ft)

  m <- fixest::feols(y ~ D | uid + tid, data = castle)
  r_model <- twfe_adequacy(m, first_treat = castle$ft)
  r_raw <- twfe_adequacy(castle$y, castle$uid, castle$tid, castle$ft)
  expect_equal(r_model$statistic$Gamma, r_raw$statistic$Gamma, tolerance = 1e-10)
  expect_equal(r_model$statistic$psi_hat, r_raw$statistic$psi_hat,
               tolerance = 1e-10)
  expect_equal(r_model$eta, r_raw$eta, tolerance = 1e-10)
  expect_identical(r_model$verdict, r_raw$verdict)
})

test_that("plm adapter reproduces the raw-data diagnostic", {
  skip_if_not_installed("plm")
  v <- get_dataset("vdem")
  d <- v[stats::complete.cases(v$ly, v$v2x_polyarchy, v$v2x_polyarchy_sd), ]

  m <- plm::plm(ly ~ v2x_polyarchy, data = d, index = c("iso", "year"),
                model = "within", effect = "twoways")
  # plm sorts its pdata.frame by panel index.  Per-observation pilot inputs to
  # a fitted-model adapter must follow that estimation-sample order.
  mi <- plm::index(m)
  model_key <- paste(as.character(mi[[1]]), as.character(mi[[2]]), sep = "\r")
  data_key <- paste(as.character(d$iso), as.character(d$year), sep = "\r")
  sd_model <- d$v2x_polyarchy_sd[match(model_key, data_key)]
  r_model <- eiv_adequacy(m, sigma_nu = sd_model, pilot = "point")
  r_raw <- eiv_adequacy(d$ly, d$v2x_polyarchy, d$iso, d$year,
                        sigma_nu = d$v2x_polyarchy_sd, pilot = "point")
  expect_equal(r_model$statistic$lambda_hat, r_raw$statistic$lambda_hat,
               tolerance = 1e-8)
  expect_equal(r_model$eta, r_raw$eta, tolerance = 1e-8)
  expect_identical(r_model$verdict, r_raw$verdict)

  # a one-way plm model must trigger BOTH guards: the specification warning
  # and the coefficient-mismatch warning (its one-way coefficient differs
  # from the two-way reproduction the diagnostic is defined on)
  m1 <- plm::plm(ly ~ v2x_polyarchy, data = d, index = c("iso", "year"),
                 model = "within", effect = "individual")
  w <- capture_warnings(eiv_adequacy(m1, sigma_nu = d$v2x_polyarchy_sd,
                                     pilot = "point"))
  expect_true(any(grepl("not a two-way within model", w)))
  expect_true(any(grepl("does not match the fitted plm coefficient", w)))
})

test_that("lm adapter reproduces the raw-data diagnostic", {
  n <- 120
  df <- data.frame(u = factor(rep(1:12, each = 10)),
                   t = factor(rep(1:10, times = 12)))
  k <- seq_len(n)
  df$x <- sin(0.9 * k) + 0.1 * as.integer(df$u)
  df$y <- 1.5 * df$x + cos(1.3 * k)

  m <- stats::lm(y ~ x + u + t, data = df)
  r_model <- leverage_report(m, x = "x", unit = "u", time = "t")
  r_raw <- leverage_report(df$y, df$x, df$u, df$t)
  expect_equal(r_model$statistic$beta, r_raw$statistic$beta, tolerance = 1e-10)
  expect_equal(r_model$statistic$se_hc2, r_raw$statistic$se_hc2,
               tolerance = 1e-10)
  expect_identical(r_model$verdict, r_raw$verdict)
  # and both match what lm itself estimated
  expect_equal(r_model$statistic$beta, unname(stats::coef(m)["x"]),
               tolerance = 1e-8)

  expect_error(leverage_report(m, x = "nope", unit = "u", time = "t"),
               "not a variable")

  # The measurement-error adapter preserves an observed nuisance control and
  # uses its leverage in the observation-specific error trace.
  df$z <- sin(0.37 * k) + 0.02 * as.integer(df$u)
  df$y2 <- 1.1 * df$x - 0.6 * df$z + cos(0.71 * k)
  sd <- 0.03 + 0.0002 * k
  mc <- stats::lm(y2 ~ x + z + u + t, data = df)
  e_model <- eiv_adequacy(mc, x = "x", unit = "u", time = "t",
                          sigma_nu = sd, pilot = "point")
  e_raw <- eiv_adequacy(df$y2, df$x, df$u, df$t, controls = df$z,
                        sigma_nu = sd, pilot = "point")
  expect_equal(e_model$statistic$beta_star, unname(stats::coef(mc)["x"]),
               tolerance = 1e-8)
  expect_equal(e_model$statistic$lambda_hat,
               e_raw$statistic$lambda_hat, tolerance = 1e-10)
})

test_that("adapter guards: multiple regressors and mismatched first_treat", {
  skip_if_not_installed("fixest")
  n <- 120
  df <- data.frame(u = rep(1:12, each = 10), t = rep(1:10, times = 12))
  k <- seq_len(n)
  df$x1 <- sin(0.9 * k); df$x2 <- cos(0.4 * k)
  df$y <- df$x1 + 0.5 * df$x2 + sin(2.2 * k)
  m2 <- fixest::feols(y ~ x1 + x2 | u + t, data = df)
  expect_error(eiv_adequacy(m2, reliability = 0.9), "exactly one")

  df$D <- as.numeric(df$u <= 4 & df$t >= 5)
  m <- fixest::feols(y ~ D | u + t, data = df)
  expect_error(twfe_adequacy(m, first_treat = rep(5, 10)),
               "estimation sample")
})
