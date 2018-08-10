context("misc")

#------------------------------------------------
# test assert_malecot_project
test_that("assert_malecot_project working correctly", {
  p <- malecot_project()
  expect_true(assert_malecot_project(p))

  expect_error(assert_malecot_project(1))
})
