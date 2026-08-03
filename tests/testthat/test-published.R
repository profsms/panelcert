test_that("published Paper A capture values are regression-locked", {
  fixture <- function(name) testthat::test_path("fixtures", name)

  match <- utils::read.csv(fixture("kss_match.csv"), check.names = FALSE)
  match_greedy <- cycle_contrasts(match$x, match$unit, match$time,
                                  method = "greedy")
  match_designed <- cycle_contrasts(match$x, match$unit, match$time,
                                    method = "sparse")
  expect_equal(match_greedy$kappa, 0.2638214123104785, tolerance = 1e-12)
  expect_equal(match_designed$kappa, 0.5110281257551292, tolerance = 1e-12)
  expect_gte(match_greedy$kappa, 0.262 - 5e-4)
  expect_gte(match_designed$kappa, 0.509 - 5e-4)

  wage <- utils::read.csv(fixture("kss_wage.csv"), check.names = FALSE)
  supports <- utils::read.csv(fixture("kss_wage_supports.csv"),
                              check.names = FALSE)
  wage_rows <- lapply(seq_len(nrow(supports)), function(i)
    as.integer(unlist(supports[i, ], use.names = FALSE)))
  wage_weights <- replicate(nrow(supports), c(1, -1, -1, 1),
                            simplify = FALSE)
  wage_system <- contrast_system(wage$x,
    list(wage$worker, wage$firm, wage$year), wage_rows, wage_weights)
  expect_equal(wage_system$kappa, 0.6112378951948423, tolerance = 1e-12)

  data("fscore", package = "panelcert")
  fscore_time <- interaction(fscore$country, fscore$year, drop = TRUE)
  fscore_greedy <- cycle_contrasts(fscore$fscore, fscore$uid, fscore_time,
                                   method = "greedy")
  fscore_designed <- cycle_contrasts(fscore$fscore, fscore$uid, fscore_time,
                                     method = "structured")
  expect_equal(fscore_greedy$kappa, 0.3741655924479094, tolerance = 1e-12)
  expect_equal(fscore_designed$kappa, 0.827387880261079, tolerance = 1e-12)
  expect_gte(fscore_greedy$kappa, 0.365 - 5e-4)
  expect_gte(fscore_designed$kappa, 0.795 - 5e-4)

  dense <- utils::read.csv(fixture("calibrated_dense.csv"),
                           check.names = FALSE)
  dense_greedy <- cycle_contrasts(dense$x, dense$unit, dense$time,
                                  method = "greedy")
  dense_designed <- cycle_contrasts(dense$x, dense$unit, dense$time,
                                    method = "structured")
  expect_equal(dense_greedy$kappa, 0.4388570425010178, tolerance = 1e-12)
  expect_equal(dense_designed$kappa, 0.834126282995861, tolerance = 1e-12)
  expect_gte(dense_greedy$kappa, 0.251 - 5e-4)
  expect_gte(dense_designed$kappa, 0.828 - 5e-4)

  data("grunfeld", package = "panelcert")
  grunfeld_design <- design_summary(grunfeld$firm, grunfeld$year,
    x = grunfeld$capital, controls = grunfeld$value)
  grunfeld_score <- score_concentration(grunfeld$invest, grunfeld$capital,
    grunfeld$firm, grunfeld$year, controls = grunfeld$value)
  grunfeld_system <- cycle_contrasts(grunfeld$capital, grunfeld$firm,
    grunfeld$year, method = "structured", controls = grunfeld$value)
  expect_equal(grunfeld_design$lambda_n, 0.20649852072450153,
               tolerance = 1e-12)
  expect_equal(grunfeld_design$n_eff, 16.559732274829294,
               tolerance = 1e-11)
  expect_equal(grunfeld_score$lambda_score, 0.738792811202239,
               tolerance = 1e-12)
  expect_equal(grunfeld_score$n_eff_score, 1.8035450833677689,
               tolerance = 1e-12)
  expect_length(grunfeld_system$rows, 32L)
  expect_equal(grunfeld_system$kappa, 0.626959467924497,
               tolerance = 1e-12)
  expect_equal(grunfeld_system$max_share, 0.3522385481988268,
               tolerance = 1e-12)
  for (i in seq_along(grunfeld_system$rows)) {
    r <- grunfeld_system$rows[[i]]
    q <- grunfeld_system$signs[[i]] / sqrt(length(r))
    expect_lte(abs(sum(q * grunfeld$value[r])), 2e-12)
    expect_lte(max(abs(rowsum(q, grunfeld$firm[r])[, 1L])), 2e-12)
    expect_lte(max(abs(rowsum(q, grunfeld$year[r])[, 1L])), 2e-12)
  }

  designs <- list(
    list(match$x, match$unit, match$time, match_greedy, match_designed),
    list(fscore$fscore, fscore$uid, fscore_time,
         fscore_greedy, fscore_designed),
    list(dense$x, dense$unit, dense$time, dense_greedy, dense_designed)
  )
  for (z in designs) {
    d <- design_summary(z[[2]], z[[3]], x = z[[1]])
    expect_true(d$lambda_n > 0 && d$lambda_n <= 1)
    expect_true(d$n_eff >= 1 && d$n_eff <= d$n)
    expect_true(z[[4]]$kappa >= 0 && z[[4]]$kappa <= z[[5]]$kappa)
    expect_lte(z[[5]]$kappa, 1 + 1e-12)
  }

  label_supports <- function(cs, unit, time) {
    sort(vapply(cs$rows, function(rows) {
      paste(sort(paste(unit[rows], time[rows], sep = ":")), collapse = "|")
    }, character(1)))
  }
  perm <- rev(seq_len(nrow(dense)))
  permuted <- cycle_contrasts(dense$x[perm], dense$unit[perm],
                              dense$time[perm], method = "structured")
  expect_equal(permuted$kappa, dense_designed$kappa, tolerance = 1e-12)
  expect_identical(label_supports(permuted, dense$unit[perm], dense$time[perm]),
                   label_supports(dense_designed, dense$unit, dense$time))
})
