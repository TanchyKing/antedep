test_that("new packaged datasets have expected structure", {
  data("cattle_growth", package = "antedep")
  data("cochlear_implant", package = "antedep")
  data("cochlear_implant_cat", package = "antedep")

  expect_true(is.list(cattle_growth))
  expect_true(is.matrix(cattle_growth$y))
  expect_equal(ncol(cattle_growth$y), 11)
  expect_equal(length(cattle_growth$blocks), nrow(cattle_growth$y))

  expect_true(is.list(cochlear_implant))
  expect_true(is.matrix(cochlear_implant$y))
  expect_equal(ncol(cochlear_implant$y), 4)
  expect_equal(length(cochlear_implant$blocks), nrow(cochlear_implant$y))
  expect_true(anyNA(cochlear_implant$y))

  expect_true(is.list(cochlear_implant_cat))
  expect_true(is.matrix(cochlear_implant_cat$y))
  expect_type(cochlear_implant_cat$y, "integer")
  expect_equal(cochlear_implant_cat$n_categories, 11L)
  expect_true(anyNA(cochlear_implant_cat$y))
})
