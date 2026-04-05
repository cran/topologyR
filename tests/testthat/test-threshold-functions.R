# ---- Threshold function tests ----

test_that("calculate_thresholds returns all five methods", {
  th <- calculate_thresholds(rnorm(50))
  expect_named(th, c("mean_diff", "median_diff", "sd", "iqr", "dbscan"))
  expect_true(all(vapply(th, is.numeric, logical(1))))
  expect_true(all(vapply(th, function(x) x > 0, logical(1))))
})

test_that("calculate_thresholds validates input", {
  expect_error(calculate_thresholds("text"), "numeric")
  expect_error(calculate_thresholds(1), "at least 2")
  expect_error(calculate_thresholds(c(1, NA)), "NA")
})

test_that("calculate_topology returns positive integer", {
  bs <- calculate_topology(rnorm(20), threshold = 0.5)
  expect_true(is.numeric(bs))
  expect_true(bs >= 2L)  # at least {} and V
})

test_that("calculate_topology validates input", {
  expect_error(calculate_topology(1:5, -1), "non-negative")
  expect_error(calculate_topology("text", 1), "numeric")
})

test_that("analyze_topology_factors returns data.frame with plot attr", {
  set.seed(42)
  result <- analyze_topology_factors(rnorm(30), plot = TRUE)
  expect_s3_class(result, "data.frame")
  expect_true("factor" %in% names(result))
  expect_true("base_size" %in% names(result))
  expect_false(is.null(attr(result, "plot")))
})

test_that("analyze_topology_factors validates input", {
  expect_error(analyze_topology_factors("text"), "numeric")
  expect_error(analyze_topology_factors(c(1, NA)), "NA")
})

test_that("visualize_topology_thresholds returns data.frame", {
  set.seed(42)
  result <- visualize_topology_thresholds(rnorm(20), plot = TRUE)
  expect_s3_class(result, "data.frame")
  expect_true("method" %in% names(result))
  plots <- attr(result, "plots")
  expect_true(is.list(plots))
  expect_equal(length(plots), 3L)
})

# ---- Legacy connectivity function tests ----

test_that("is_topology_connected works on basic examples", {
  expect_true(is_topology_connected(list(c(1, 2, 3))))
  expect_false(is_topology_connected(list()))
  expect_true(is_topology_connected(list(c(1, 2), c(2, 3))))
})

test_that("is_topology_connected2 works on basic examples", {
  expect_true(is_topology_connected2(list(c(1, 2, 3))))
  expect_false(is_topology_connected2(list()))
})

test_that("is_topology_connected_manual checks coverage", {
  expect_true(is_topology_connected_manual(list(c(1, 2, 3), c(3, 4, 5))))
  expect_false(is_topology_connected_manual(list(c(1, 2), c(4, 5))))
})
