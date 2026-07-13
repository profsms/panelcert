# The bundled datasets are the article-reproducibility surface: after
# library(panelcert) they must be present, correctly shaped, and drive the
# diagnostics directly with no external download.

test_that("all four datasets ship and are correctly shaped", {
  vdem    <- get_dataset("vdem")
  psid    <- get_dataset("psid")
  castle  <- get_dataset("castle")
  divorce <- get_dataset("divorce")

  expect_s3_class(vdem, "data.frame")
  expect_equal(nrow(vdem), 8931)
  expect_true(all(c("iso", "year", "ly", "v2x_polyarchy", "v2x_polyarchy_sd") %in%
                    names(vdem)))
  expect_equal(dim(psid), c(4165L, 4L))
  expect_equal(dim(castle), c(550L, 4L))
  expect_equal(dim(divorce), c(1377L, 4L))
  # ft convention preserved: NA = never-treated
  expect_equal(sum(is.na(castle$ft)), 319)   # 29 never-treated states x 11 years
  expect_equal(sum(is.na(divorce$ft)), 540)  # 20 never-treated states x 27 years
})

test_that("the paper applications run straight off the bundled data", {
  # Paper B lead application, verbatim from the dataset
  d <- vdem[stats::complete.cases(vdem$ly, vdem$v2x_polyarchy, vdem$v2x_polyarchy_sd), ]
  rb <- eiv_adequacy(d$ly, d$v2x_polyarchy, d$iso, d$year, sigma_nu = d$v2x_polyarchy_sd)
  expect_identical(rb$verdict, "CERTIFIED")

  # Paper C flagged pole, verbatim from the dataset
  rc <- twfe_adequacy(divorce$y, divorce$uid, divorce$tid, divorce$ft)
  expect_identical(rc$verdict, "FLAGGED")
})
