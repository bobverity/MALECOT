context("assertions")

#------------------------------------------------
# test assert_numeric
test_that("assert_numeric works", {
  expect_true(assert_numeric(5))
  expect_true(assert_numeric(-5:5))

  expect_error(assert_numeric("foo"))
})

#------------------------------------------------
# test assert_pos
test_that("assert_pos works", {
  expect_true(assert_pos(seq(1, 5, 0.5)))
  expect_true(assert_pos(seq(0, 5, 0.5), zero_allowed = TRUE))

  expect_error(assert_pos(seq(-1, -5, -0.5)))
  expect_error(assert_pos(seq(0, 5, 0.5), zero_allowed = FALSE))
  expect_error(assert_pos(seq(-5, 5, 0.5), zero_allowed = TRUE))
  expect_error(assert_pos("foo"))
})

#------------------------------------------------
# test assert_int
test_that("assert_int works", {
  expect_true(assert_int(-5:5))

  expect_error(assert_int(0.5))
  expect_error(assert_int("foo"))
})

#------------------------------------------------
# test assert_pos_int
test_that("assert_pos_int works", {
  expect_true(assert_pos_int(1:5))
  expect_true(assert_pos_int(0:5, zero_allowed = TRUE))

  expect_error(assert_pos_int(-1:-5))
  expect_error(assert_pos_int(0:5, zero_allowed = FALSE))
  expect_error(assert_pos_int(-5:5, zero_allowed = TRUE))
  expect_error(assert_pos_int("foo"))
})

#------------------------------------------------
# test assert_bounded
test_that("assert_bounded works", {
  expect_true(assert_bounded(seq(0, 1, 0.1)))
  expect_true(assert_bounded(seq(-5, 5, 0.1), left = -5, right = 5))

  expect_error(assert_bounded(-1, left = 0))
  expect_error(assert_bounded(2, right = 0))
  expect_error(assert_bounded(0, left = 0, inclusive_left = FALSE))
  expect_error(assert_bounded(1, right = 1, inclusive_right = FALSE))
  expect_error(assert_bounded("foo"))
})

#------------------------------------------------
# test assert_same_length
test_that("assert_same_length works", {
  expect_true(assert_same_length(1:3, c("a", "b", "c"), list("foo", 1, 0.5)))

  expect_error(assert_same_length(1:3, c("a", "b")))
})

#------------------------------------------------
# test assert_matrix
test_that("assert_matrix works", {
  m <- matrix(1)
  expect_true(assert_matrix(m))

  expect_error(assert_matrix(1))
  expect_error(assert_matrix(1:5))
  expect_error(assert_matrix("foo"))
})

#------------------------------------------------
# test assert_square_matrix
test_that("assert_square_matrix works", {
  expect_true(assert_square_matrix(matrix(0, 3, 3)))

  expect_error(assert_square_matrix(matrix(0, 2, 3)))
  expect_error(assert_square_matrix(1))
  expect_error(assert_square_matrix(1:5))
  expect_error(assert_square_matrix("foo"))
})

#------------------------------------------------
# test assert_symmetric_matrix
test_that("assert_symmetric_matrix works", {
  m0 <- matrix(1:16, 4, 4)
  m1 <- m0 + t(m0)
  expect_true(assert_symmetric_matrix(m1))

  expect_error(assert_symmetric_matrix(m0))
  expect_error(assert_symmetric_matrix(1))
  expect_error(assert_symmetric_matrix(1:5))
  expect_error(assert_symmetric_matrix("foo"))
})

#------------------------------------------------
# test assert_malecot_project
test_that("assert_malecot_project works", {
  p <- malecot_project()
  expect_true(assert_malecot_project(p))

  expect_error(assert_malecot_project(1))
})
