
#------------------------------------------------
# is numeric
assert_numeric <- function(x, name = deparse(substitute(x))) {
  if (!is.numeric(x)) {
    stop(sprintf("'%s' must be numeric", name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is positive (with or without zero allowed)
assert_pos <- function(x, zero_allowed = TRUE, name = deparse(substitute(x))) {
  assert_numeric(x, name)
  if (zero_allowed) {
    if (!all(x>=0)) {
      stop(sprintf("'%s' must be greater than or equal to zero", name), call. = FALSE)
    }
  } else {
    if (!all(x>0)) {
      stop(sprintf("'%s' must be greater than zero", name), call. = FALSE)
    }
  }
  return(TRUE)
}

#------------------------------------------------
# is integer
assert_int <- function(x, name = deparse(substitute(x))) {
  assert_numeric(x, name)
  if (!isTRUE(all.equal(x, as.integer(x)))) {
    stop(sprintf("'%s' must be integer valued", name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is positive integer (with or without zero allowed)
assert_pos_int <- function(x, zero_allowed = TRUE, name = deparse(substitute(x))) {
  assert_int(x, name)
  assert_pos(x, zero_allowed, name)
  return(TRUE)
}

#------------------------------------------------
# is between bounds (inclusive or exclusive)
assert_bounded <- function(x, left = 0, right = 1, inclusive_left = TRUE, inclusive_right = TRUE, name = deparse(substitute(x))) {
  assert_numeric(x, name)
  if (inclusive_left) {
    if (!all(x>=left)) {
      stop(sprintf("'%s' must be greater than or equal to %s", name, left), call. = FALSE)
    }
  } else {
    if (!all(x>left)) {
      stop(sprintf("'%s' must be greater than %s", name, left), call. = FALSE)
    }
  }
  if (inclusive_right) {
    if (!all(x<=right)) {
      stop(sprintf("'%s' must be less than or equal to %s", name, right), call. = FALSE)
    }
  } else {
    if (!all(x<right)) {
      stop(sprintf("'%s' must be less than %s", name, right), call. = FALSE)
    }
  }
  return(TRUE)
}

#------------------------------------------------
# objects all same length
assert_same_length <- function(...) {
  l <- mapply(length, list(...))
  if (!length(unique(l)) == 1) {
    dots <- match.call(expand.dots = FALSE)$...
    dot_names <- paste(sapply(dots, deparse), collapse=", ")
    stop(sprintf("variables %s must be the same length", dot_names), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is matrix
assert_matrix <- function(x, name = deparse(substitute(x))) {
  if (!is.matrix(x)) {
    stop(sprintf("'%s' must be a matrix", name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is square matrix
assert_square_matrix <- function(x, name = deparse(substitute(x))) {
  assert_matrix(x, name)
  if (nrow(x) != ncol(x)) {
    stop(sprintf("'%s' must be a square matrix", name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is symmetric matrix
assert_symmetric_matrix <- function(x, name = deparse(substitute(x))) {
  assert_square_matrix(x, name)
  if (!isSymmetric(x)) {
    stop(sprintf("'%s' must be a symmetric matrix", name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is class covfefe_project
#assert_covfefe_project <- function(x, name = deparse(substitute(x))) {
#  if (!is.covfefe_project(x)) {
#    stop(sprintf("'%s' must be of class 'covfefe_project'", name), call. = FALSE)
#  }
#  return(TRUE)
#}
