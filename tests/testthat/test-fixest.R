# External consistency: the package's Frisch-Waugh reproduction of the user's
# FE regression must agree with fixest::feols on the same data — the package
# CONSUMES an already-estimated model and must never disagree with it (spec 0).

test_that("FWL core agrees with fixest::feols on the V-Dem panel", {
  skip_if_not_installed("fixest")
  v <- get_dataset("vdem")
  keep <- stats::complete.cases(v$ly, v$v2x_polyarchy, v$v2x_polyarchy_sd)
  d <- v[keep, ]

  m <- fixest::feols(ly ~ v2x_polyarchy | iso + year, data = d)
  r <- eiv_adequacy(d$ly, d$v2x_polyarchy, d$iso, d$year,
                    sigma_nu = d$v2x_polyarchy_sd, pilot = "point")

  # same point estimate the user already has
  expect_equal(unname(stats::coef(m)["v2x_polyarchy"]), r$statistic$beta_star,
               tolerance = 1e-8)

  # same within transformation: fixest::demean vs alternating projections
  xt_fixest <- as.numeric(fixest::demean(d["v2x_polyarchy"],
                                         d[c("iso", "year")])[, 1])
  xt_ours <- twoway_demean(d$v2x_polyarchy, d$iso, d$year)
  expect_lt(max(abs(xt_fixest - xt_ours)), 1e-7)
  expect_equal(sum(xt_fixest^2), r$design$tau_star2, tolerance = 1e-8)
})
