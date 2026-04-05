# ---- HVG tests ----

test_that("HVG handles empty and single-element series", {
  g0 <- horizontal_visibility_graph(numeric(0))
  expect_equal(g0$n, 0L)
  expect_equal(g0$n_edges, 0L)

  g1 <- horizontal_visibility_graph(42)
  expect_equal(g1$n, 1L)
  expect_equal(g1$n_edges, 0L)
})

test_that("HVG adjacent points are always connected", {
  g <- horizontal_visibility_graph(c(5, 3, 1, 4, 2))
  # All consecutive pairs must be edges
  for (i in 1:4) {
    has_edge <- any(
      (g$edges$from == i & g$edges$to == i + 1L) |
      (g$edges$from == i + 1L & g$edges$to == i)
    )
    expect_true(has_edge, info = paste("Adjacent edge", i, i + 1))
  }
})

test_that("HVG blocks correctly when intermediate value >= min", {
  # y = [5, 3, 3]: edge (1,3) should NOT exist because y[2]=3 >= min(5,3)=3
  g <- horizontal_visibility_graph(c(5, 3, 3))
  edge_1_3 <- any(
    (g$edges$from == 1 & g$edges$to == 3) |
    (g$edges$from == 3 & g$edges$to == 1)
  )
  expect_false(edge_1_3)
})

test_that("HVG allows edges when intermediates are strictly below min", {
  # y = [3, 1, 2, 5]: edge (1,3) exists because y[2]=1 < min(3,2)=2
  g <- horizontal_visibility_graph(c(3, 1, 2, 5))
  edge_1_3 <- any(g$edges$from == 1 & g$edges$to == 3)
  expect_true(edge_1_3)

  # y = [3, 1, 2, 5]: edge (1,4) exists because all intermediates < min(3,5)=3
  edge_1_4 <- any(
    (g$edges$from == 1 & g$edges$to == 4) |
    (g$edges$from == 4 & g$edges$to == 1)
  )
  expect_true(edge_1_4)
})

test_that("HVG handles equal values correctly", {
  # All equal: only adjacent edges
  g <- horizontal_visibility_graph(c(3, 3, 3, 3))
  expect_equal(g$n_edges, 3L)  # only 3 adjacent pairs
})

test_that("HVG known example: exact edge set", {
  # y = [2, 1, 2, 1, 3]
  g <- horizontal_visibility_graph(c(2, 1, 2, 1, 3))
  expected_edges <- matrix(c(
    1, 2,
    1, 3,
    2, 3,
    3, 4,
    3, 5,
    4, 5
  ), ncol = 2, byrow = TRUE)

  for (k in seq_len(nrow(expected_edges))) {
    a <- expected_edges[k, 1]
    b <- expected_edges[k, 2]
    has <- any(
      (g$edges$from == a & g$edges$to == b) |
      (g$edges$from == b & g$edges$to == a)
    )
    expect_true(has, info = paste("Expected edge", a, b))
  }
  expect_equal(g$n_edges, 6L)
})


# ---- NVG tests ----

test_that("NVG handles empty and single-element series", {
  g0 <- natural_visibility_graph(numeric(0))
  expect_equal(g0$n, 0L)
  expect_equal(g0$n_edges, 0L)

  g1 <- natural_visibility_graph(42)
  expect_equal(g1$n, 1L)
  expect_equal(g1$n_edges, 0L)
})

test_that("NVG has at least as many edges as HVG", {
  set.seed(123)
  series <- rnorm(50)
  g_hvg <- horizontal_visibility_graph(series)
  g_nvg <- natural_visibility_graph(series)
  expect_gte(g_nvg$n_edges, g_hvg$n_edges)
})

test_that("NVG adjacent points are always connected", {
  g <- natural_visibility_graph(c(5, 3, 1, 4, 2))
  for (i in 1:4) {
    has_edge <- any(
      (g$edges$from == i & g$edges$to == i + 1L) |
      (g$edges$from == i + 1L & g$edges$to == i)
    )
    expect_true(has_edge, info = paste("Adjacent edge", i, i + 1))
  }
})

test_that("NVG blocks when intermediate is above line", {
  # y = [3, 8, 1, 2, 10]: point 0 is NOT visible from point 4
  # because the line from (0,3) to (4,10) at t=1 is 4.75, y[1]=8 >= 4.75
  g <- natural_visibility_graph(c(3, 8, 1, 2, 10))
  edge_1_5 <- any(
    (g$edges$from == 1 & g$edges$to == 5) |
    (g$edges$from == 5 & g$edges$to == 1)
  )
  expect_false(edge_1_5)
})

test_that("NVG known small example", {
  g <- natural_visibility_graph(c(3, 1, 2, 5))
  # All pairs visible: (1,2), (1,3), (1,4), (2,3), (2,4), (3,4)
  expect_equal(g$n_edges, 6L)
})


# ---- Input validation ----

test_that("visibility graph functions validate input", {
  expect_error(horizontal_visibility_graph("text"), "numeric")
  expect_error(horizontal_visibility_graph(c(1, NA)), "NA")
  expect_error(horizontal_visibility_graph(c(1, Inf)), "Inf")
  expect_error(natural_visibility_graph("text"), "numeric")
  expect_error(natural_visibility_graph(c(1, NA)), "NA")
})

test_that("adjacency list is symmetric", {
  g <- horizontal_visibility_graph(c(3, 1, 4, 1, 5))
  for (i in seq_along(g$adjacency)) {
    for (j in g$adjacency[[i]]) {
      expect_true(i %in% g$adjacency[[j]],
                  info = paste(i, "in adj of", j, "but not vice versa"))
    }
  }
})
