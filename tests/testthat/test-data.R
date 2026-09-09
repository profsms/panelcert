# The bundled datasets are the article-reproducibility surface: after
# library(panelcert) they must be present, correctly shaped, and drive the
# diagnostics directly with no external download.

test_that("all bundled datasets ship and are correctly shaped", {
  vdem    <- get_dataset("vdem")
  psid    <- get_dataset("psid")
  twins   <- get_dataset("twins")
  castle  <- get_dataset("castle")
  divorce <- get_dataset("divorce")
  minimum_wage <- get_dataset("minimum_wage")
  brazil <- get_dataset("brazil")
  fscore  <- get_dataset("fscore")

  expect_s3_class(vdem, "data.frame")
  expect_equal(nrow(vdem), 8931)
  expect_true(all(c("iso", "year", "ly", "v2x_polyarchy", "v2x_polyarchy_sd") %in%
                    names(vdem)))
  expect_equal(dim(psid), c(4165L, 4L))
  expect_equal(dim(twins), c(183L, 16L))
  expect_equal(sum(stats::complete.cases(twins[c("DLHRWAGE", "DEDUC1", "DEDUC2",
                                                 "DTEN", "DMARRIED", "DUNCOV")])),
               147L)
  expect_equal(dim(castle), c(550L, 4L))
  expect_equal(dim(divorce), c(1377L, 4L))
  expect_equal(dim(minimum_wage), c(15988L, 4L))
  expect_equal(length(unique(minimum_wage$uid)), 2284L)
  expect_equal(length(unique(minimum_wage$tid)), 7L)
  expect_equal(dim(brazil), c(34080L, 4L))
  expect_equal(length(unique(brazil$uid)), 2840L)
  expect_equal(length(unique(brazil$tid)), 12L)
  # ft convention preserved: NA = never-treated
  expect_equal(sum(is.na(castle$ft)), 319)   # 29 never-treated states x 11 years
  expect_equal(sum(is.na(divorce$ft)), 540)  # 20 never-treated states x 27 years
  expect_true(all(is.na(brazil$ft) | brazil$ft >= 2L))

  # Paper A firm panel
  expect_equal(dim(fscore), c(217L, 6L))
  expect_equal(length(unique(fscore$uid)), 19L)
  expect_true(all(c("uid", "year", "country", "status", "fscore", "ret") %in%
                    names(fscore)))
})

test_that("the paper applications run straight off the bundled data", {
  # Paper B lead application, verbatim from the dataset
  d <- vdem[stats::complete.cases(vdem$ly, vdem$v2x_polyarchy, vdem$v2x_polyarchy_sd), ]
  rb <- eiv_adequacy(d$ly, d$v2x_polyarchy, d$iso, d$year, sigma_nu = d$v2x_polyarchy_sd)
  expect_identical(rb$verdict, "CERTIFIED")

  # Paper C long-panel application, verbatim from the dataset. Its lower bound
  # does not clear the threshold and the saturated upper procedure is
  # unavailable because the projected score covariance is rank deficient.
  rc <- twfe_adequacy(divorce$y, divorce$uid, divorce$tid, divorce$ft)
  expect_identical(rc$verdict, "INCONCLUSIVE")

  # Diffuse-regime companion application, verbatim from the dataset.
  ra <- leverage_report(log1p(fscore$ret), fscore$fscore, fscore$uid, fscore$year)
  expect_equal(ra$design$n, 217L)
  expect_equal(round(ra$statistic$beta, 5), -0.00147)   # paper -0.00147
  expect_equal(round(ra$statistic$se_hc2, 4), 0.0153)   # HC2 SE, paper 0.0153
  expect_equal(round(ra$design$tau_star2), 559)         # pooled tau^2 ~ 559
  mE <- fscore$country == "Poland"                       # Spec E headline |t|
  re <- leverage_report(log1p(fscore$ret[mE]), fscore$fscore[mE],
                        fscore$uid[mE], fscore$year[mE])
  expect_equal(round(abs(re$statistic$beta / re$statistic$se_hc2), 2), 0.61)
})
