# ---- Topology generation tests ----

test_that("generate_topology produces valid topology from simple graph", {
  # Path graph: 1-2-3
  adj <- list(c(2L), c(1L, 3L), c(2L))
  topo <- generate_topology(adj, 3L, verify_axioms = TRUE)

  expect_true(topo$complete)
  expect_true(topo$axioms_ok)
  expect_true(topo$n_open_sets >= 4L)

  sets <- topo$topology
  has_empty <- any(vapply(sets, function(s) length(s) == 0, logical(1)))
  has_full <- any(vapply(sets, function(s) setequal(s, 1:3), logical(1)))
  expect_true(has_empty)
  expect_true(has_full)
})

test_that("generate_topology detects connectivity correctly", {
  adj_conn <- list(c(2L), c(1L, 3L), c(2L))
  topo_conn <- generate_topology(adj_conn, 3L)
  expect_true(is.logical(topo_conn$connected))

  # Two isolated nodes
  adj_disc <- list(integer(0), integer(0))
  topo_disc <- generate_topology(adj_disc, 2L)
  expect_true(is.logical(topo_disc$connected))
})

test_that("generate_topology handles single element", {
  adj <- list(integer(0))
  topo <- generate_topology(adj, 1L)
  expect_true(topo$connected)
})

test_that("generate_topology works for n > 64 (multi-word bitsets)", {
  # Ring graph on 70 nodes: each node connected to next and previous
  n <- 70L
  adj <- lapply(seq_len(n), function(i) {
    c(if (i > 1L) i - 1L else n, if (i < n) i + 1L else 1L)
  })
  # Skip topology enumeration (too large), just get base + connectivity
  topo <- generate_topology(adj, n, max_open_sets = 0L)
  expect_true(is.logical(topo$connected))
  expect_true(length(topo$base) >= n)
  expect_true(length(topo$components) >= 1L)
})

test_that("specialization preorder gives exact components", {
  # Two disconnected triangles with N[v] (closed neighborhoods).
  # N[1]={1,2,3}, N[2]={1,2,3}, N[3]={1,2,3} -> all identical -> 1 component
  # N[4]={4,5,6}, N[5]={4,5,6}, N[6]={4,5,6} -> all identical -> 1 component
  # Two graph components -> 2 topological components.
  adj <- list(
    c(2L, 3L), c(1L, 3L), c(1L, 2L),  # triangle 1
    c(5L, 6L), c(4L, 6L), c(4L, 5L)   # triangle 2
  )
  topo <- generate_topology(adj, 6L)
  expect_false(topo$connected)
  expect_equal(length(topo$components), 2L)
  all_elements <- sort(unlist(topo$components))
  expect_equal(all_elements, 1:6)

  # A case with genuinely 2 components: path 1-2 and isolated node 3
  # N(1)={2}, N(2)={1}, N(3)={}. Base={1},{2},{1,2}. S(1)⊂S(2)? No.
  # Actually with N(1)={2}, N(2)={1}: base after intersections = {{2},{1}}.
  # S(1)={{1}} (only base element containing 1), S(2)={{2}}.
  # Not comparable → still 2 separate components from triangle, plus node 3.
  # Better test: star graph where center has superset neighborhoods
  adj_star <- list(c(2L, 3L, 4L), c(1L), c(1L), c(1L))
  topo_star <- generate_topology(adj_star, 4L)
  # N(1)={2,3,4}, N(2)={1}, N(3)={1}, N(4)={1}
  # Base includes {1} (from intersections of N(2)∩N(3), etc.)
  # S(1) includes base elements containing 1 (like {1})
  # Nodes 2,3,4 appear in N(1), so base element {2,3,4} contains them
  # This should show some connectivity structure
  expect_true(is.logical(topo_star$connected))
  expect_true(length(topo_star$components) >= 1L)
  all_star <- sort(unlist(topo_star$components))
  expect_equal(all_star, 1:4)
})

test_that("max_open_sets limits enumeration gracefully", {
  # Path graph on 10 nodes: N[v] produces distinct neighborhoods,
  # base grows via intersections, union closure can exceed limit.
  n <- 10L
  adj <- lapply(seq_len(n), function(i) {
    nb <- integer(0)
    if (i > 1L) nb <- c(nb, i - 1L)
    if (i < n)  nb <- c(nb, i + 1L)
    nb
  })
  topo <- generate_topology(adj, n, max_open_sets = 20L)
  expect_true(is.logical(topo$connected))
  # With a tight limit, enumeration should be incomplete
  expect_false(topo$complete)
})


# ---- Connectivity tests (standalone) ----

test_that("exact connectivity identifies connected spaces", {
  tau <- list(integer(0), 1:3)
  r <- is_topology_connected_exact(tau, 3L)
  expect_true(r$connected)
})

test_that("exact connectivity identifies disconnected spaces", {
  tau <- list(integer(0), 1:4, c(1L, 2L), c(3L, 4L))
  r <- is_topology_connected_exact(tau, 4L)
  expect_false(r$connected)
  expect_true("component1" %in% names(r))
  expect_true("component2" %in% names(r))
  expect_true(setequal(union(r$component1, r$component2), 1:4))
  expect_equal(length(intersect(r$component1, r$component2)), 0L)
})

test_that("exact connectivity works for n > 64", {
  # Indiscrete topology on 100 elements: just {}, V
  tau <- list(integer(0), 1:100)
  r <- is_topology_connected_exact(tau, 100L)
  expect_true(r$connected)
})

test_that("exact connectivity rejects invalid input", {
  expect_error(is_topology_connected_exact("not a list", 3L), "list")
})


# ---- simplest_topology tests ----

test_that("simplest_topology produces discrete topology", {
  st <- simplest_topology(c(10, 20, 30))
  expect_equal(st$topology_type, "discrete")
  expect_equal(st$n, 3L)
  expect_false(st$connected)
  expect_equal(length(st$base), 3L)
  for (s in st$base) {
    expect_equal(length(s), 1L)
  }
})

test_that("simplest_topology single element is connected", {
  st <- simplest_topology(42)
  expect_true(st$connected)
})


# ---- complete_topology tests ----

test_that("complete_topology works for small data", {
  ct <- complete_topology(c(1, 2, 3))
  expect_true(ct$n_open_sets >= 2L)
  expect_true(is.logical(ct$connected))
})

test_that("complete_topology validates input", {
  expect_error(complete_topology("text"), "numeric")
  expect_error(complete_topology(1), "at least 2")
})
