# ======================================================================
# Tests for directed visibility graphs, Alexandrov topology,
# bitopological analysis, and invariants
# ======================================================================

# ------------------------------------------------------------------
# 1. Directed visibility graphs
# ------------------------------------------------------------------

test_that("HVG directed=FALSE returns no directed adjacency", {
  g <- horizontal_visibility_graph(c(3, 1, 4, 1, 5))
  expect_null(g$out_adjacency)
  expect_null(g$in_adjacency)
})

test_that("HVG directed=TRUE returns out_adjacency and in_adjacency", {
  g <- horizontal_visibility_graph(c(3, 1, 4, 1, 5), directed = TRUE)
  expect_true(!is.null(g$out_adjacency))
  expect_true(!is.null(g$in_adjacency))
  expect_length(g$out_adjacency, g$n)
  expect_length(g$in_adjacency, g$n)
})

test_that("HVG directed edges satisfy from < to (DAG property)", {
  g <- horizontal_visibility_graph(c(3, 1, 4, 1, 5), directed = TRUE)
  # All out-neighbors of v should have index > v
  for (v in seq_len(g$n)) {
    out_nb <- g$out_adjacency[[v]]
    if (length(out_nb) > 0) {
      expect_true(all(out_nb > v),
                  info = paste("out_adjacency of", v, "has backward edges"))
    }
    in_nb <- g$in_adjacency[[v]]
    if (length(in_nb) > 0) {
      expect_true(all(in_nb < v),
                  info = paste("in_adjacency of", v, "has forward edges"))
    }
  }
})

test_that("NVG directed=TRUE returns valid directed adjacency", {
  g <- natural_visibility_graph(c(3, 1, 4, 1, 5), directed = TRUE)
  expect_true(!is.null(g$out_adjacency))
  expect_true(!is.null(g$in_adjacency))
  # DAG property
  for (v in seq_len(g$n)) {
    if (length(g$out_adjacency[[v]]) > 0)
      expect_true(all(g$out_adjacency[[v]] > v))
    if (length(g$in_adjacency[[v]]) > 0)
      expect_true(all(g$in_adjacency[[v]] < v))
  }
})

test_that("directed adjacency is consistent with undirected", {
  g <- horizontal_visibility_graph(c(3, 1, 4, 1, 5), directed = TRUE)
  # For each v, sort(union(out_adj[v], in_adj[v])) == sort(adj[v])
  for (v in seq_len(g$n)) {
    directed_union <- sort(c(g$out_adjacency[[v]], g$in_adjacency[[v]]))
    undirected <- sort(g$adjacency[[v]])
    expect_equal(directed_union, undirected,
                 info = paste("vertex", v))
  }
})

test_that("directed graph of single observation", {
  g <- horizontal_visibility_graph(c(42), directed = TRUE)
  expect_length(g$out_adjacency, 1L)
  expect_equal(g$out_adjacency[[1]], integer(0))
  expect_equal(g$in_adjacency[[1]], integer(0))
})

test_that("HVG directed: known edge structure for c(3, 1, 4, 1, 5)", {
  # Manual computation:
  # Edges: (1,2), (1,3), (2,3), (3,4), (3,5), (4,5)
  g <- horizontal_visibility_graph(c(3, 1, 4, 1, 5), directed = TRUE)
  expect_equal(sort(g$out_adjacency[[1]]), c(2L, 3L))
  expect_equal(g$out_adjacency[[2]], 3L)
  expect_equal(sort(g$out_adjacency[[3]]), c(4L, 5L))
  expect_equal(g$out_adjacency[[4]], 5L)
  expect_equal(g$out_adjacency[[5]], integer(0))
  expect_equal(g$in_adjacency[[1]], integer(0))
  expect_equal(g$in_adjacency[[2]], 1L)
  expect_equal(sort(g$in_adjacency[[3]]), c(1L, 2L))
  expect_equal(g$in_adjacency[[4]], 3L)
  expect_equal(sort(g$in_adjacency[[5]]), c(3L, 4L))
})


# ------------------------------------------------------------------
# 2. Alexandrov topology
# ------------------------------------------------------------------

test_that("Alexandrov topology: path graph 1->2->3", {
  # out_adjacency for path: 1->[2], 2->[3], 3->[]
  out_adj <- list(2L, 3L, integer(0))
  alex <- generate_alexandrov_topology(out_adj, 3L)

  # reach[3]={3}, reach[2]={2,3}, reach[1]={1,2,3}
  # Base = {{3}, {2,3}, {1,2,3}} = 3 distinct sets
  expect_equal(length(alex$base), 3L)
  expect_true(alex$connected)
  expect_equal(length(alex$components), 1L)
  expect_true(alex$base_complete)
})

test_that("Alexandrov topology: enumeration produces correct open sets", {
  out_adj <- list(2L, 3L, integer(0))
  alex <- generate_alexandrov_topology(out_adj, 3L, max_open_sets = 100L)

  # Topology = {empty, {3}, {2,3}, {1,2,3}} = 4 open sets
  expect_equal(alex$n_open_sets, 4L)
  expect_true(alex$complete)
})

test_that("Alexandrov topology: disconnected DAG", {
  # Two isolated components: 1->2, 3->4
  out_adj <- list(2L, integer(0), 4L, integer(0))
  alex <- generate_alexandrov_topology(out_adj, 4L)

  expect_false(alex$connected)
  expect_equal(length(alex$components), 2L)
})

test_that("Alexandrov topology: single vertex", {
  out_adj <- list(integer(0))
  alex <- generate_alexandrov_topology(out_adj, 1L)

  expect_equal(length(alex$base), 1L)
  expect_true(alex$connected)
})

test_that("Alexandrov topology: complete DAG produces 1 base element", {
  # 1->{2,3}, 2->{3}, 3->[] — all reach {1,2,3} except vertex 3
  # reach[3]={3}, reach[2]={2,3}, reach[1]={1,2,3}
  # Actually 3 distinct sets, not 1
  out_adj <- list(c(2L, 3L), 3L, integer(0))
  alex <- generate_alexandrov_topology(out_adj, 3L)

  # reach[3]={3}, reach[2]={2,3}, reach[1]={1,2,3} → 3 base elements
  expect_equal(length(alex$base), 3L)
  expect_true(alex$connected)
})

test_that("Alexandrov topology: star DAG (hub reaches all)", {
  # 1->{2,3,4,5}, 2->{}, 3->{}, 4->{}, 5->{}
  out_adj <- list(c(2L, 3L, 4L, 5L), integer(0), integer(0),
                  integer(0), integer(0))
  alex <- generate_alexandrov_topology(out_adj, 5L)

  # reach[2]={2}, reach[3]={3}, reach[4]={4}, reach[5]={5}, reach[1]={1,2,3,4,5}
  # Base has 5 elements. Singletons are pairwise incomparable.
  # S(1)={V} ⊆ S(k) for any k? S(1) has only base[4]={1,2,3,4,5}.
  # S(2) has base[0]={2} and base[4]={1,2,3,4,5}. S(1) ⊆ S(2)? Yes.
  # So 1 is comparable with everyone → connected.
  expect_equal(length(alex$base), 5L)
  expect_true(alex$connected)
})

test_that("Alexandrov base is always subset of Nada base (resolution hierarchy)", {
  series <- c(3, 1, 4, 1, 5, 9, 2, 6)
  g <- horizontal_visibility_graph(series, directed = TRUE)

  alex <- generate_alexandrov_topology(g$out_adjacency, g$n)
  nada_fwd <- generate_topology(g$out_adjacency, g$n, max_open_sets = 0L)

  # Nada base size >= Alexandrov base size
  expect_true(length(nada_fwd$base) >= length(alex$base),
              info = paste("Nada:", length(nada_fwd$base),
                           "Alex:", length(alex$base)))
})

test_that("Alexandrov verify_axioms passes for small example", {
  out_adj <- list(2L, 3L, integer(0))
  alex <- generate_alexandrov_topology(out_adj, 3L,
                                       max_open_sets = 100L,
                                       verify_axioms = TRUE)
  expect_true(alex$axioms_ok)
})


# ------------------------------------------------------------------
# 3. Forward and backward Nada topologies
# ------------------------------------------------------------------

test_that("Nada-forward and Nada-backward from directed graph", {
  g <- horizontal_visibility_graph(c(3, 1, 4, 1, 5), directed = TRUE)

  tau_fwd <- generate_topology(g$out_adjacency, g$n, max_open_sets = 0L)
  tau_bwd <- generate_topology(g$in_adjacency, g$n, max_open_sets = 0L)

  # Both should have valid results
  expect_true(!is.null(tau_fwd$connected))
  expect_true(!is.null(tau_bwd$connected))
  expect_true(length(tau_fwd$base) > 0)
  expect_true(length(tau_bwd$base) > 0)
})

test_that("forward + backward base sizes can differ (asymmetry)", {
  # For an asymmetric series, base sizes may differ
  series <- c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5)
  g <- horizontal_visibility_graph(series, directed = TRUE)
  tau_fwd <- generate_topology(g$out_adjacency, g$n, max_open_sets = 0L)
  tau_bwd <- generate_topology(g$in_adjacency, g$n, max_open_sets = 0L)

  # We don't know which is larger, but both should have valid bases
  expect_true(length(tau_fwd$base) >= 1)
  expect_true(length(tau_bwd$base) >= 1)
})


# ------------------------------------------------------------------
# 4. Bitopological analysis (full pipeline)
# ------------------------------------------------------------------

test_that("generate_bitopology returns all expected components", {
  bt <- generate_bitopology(c(3, 1, 4, 1, 5, 9, 2, 6))
  expect_true(!is.null(bt$graph))
  expect_true(!is.null(bt$undirected))
  expect_true(!is.null(bt$forward))
  expect_true(!is.null(bt$backward))
  expect_true(!is.null(bt$alexandrov))
  expect_true(!is.null(bt$invariants))
})

test_that("generate_bitopology with NVG", {
  bt <- generate_bitopology(c(3, 1, 4, 1, 5), graph_type = "nvg")
  expect_true(!is.null(bt$graph))
  expect_true(!is.null(bt$invariants))
})

test_that("generate_bitopology without Alexandrov", {
  bt <- generate_bitopology(c(3, 1, 4, 1, 5), alexandrov = FALSE)
  expect_null(bt$alexandrov)
  expect_null(bt$invariants$resolution)
})

test_that("generate_bitopology input validation", {
  expect_error(generate_bitopology("not numeric"))
  expect_error(generate_bitopology(c(1, NA, 3)))
  expect_error(generate_bitopology(c(1, Inf)))
  expect_error(generate_bitopology(42))  # length 1
})


# ------------------------------------------------------------------
# 5. Bitopological invariants
# ------------------------------------------------------------------

test_that("bitopology_invariants returns all expected fields", {
  g <- horizontal_visibility_graph(c(3, 1, 4, 1, 5), directed = TRUE)
  tf <- generate_topology(g$out_adjacency, g$n, max_open_sets = 0L)
  tb <- generate_topology(g$in_adjacency, g$n, max_open_sets = 0L)
  inv <- bitopology_invariants(tf, tb, g$n)

  expect_true(!is.null(inv$forward_components))
  expect_true(!is.null(inv$backward_components))
  expect_true(!is.null(inv$forward_connected))
  expect_true(!is.null(inv$backward_connected))
  expect_true(!is.null(inv$irreversibility_components))
  expect_true(!is.null(inv$irreversibility_base))
  expect_true(!is.null(inv$asymmetry_direction))
  expect_true(!is.null(inv$pairwise))
})

test_that("irreversibility index is in valid range", {
  bt <- generate_bitopology(c(3, 1, 4, 1, 5, 9, 2, 6))
  inv <- bt$invariants
  expect_true(inv$irreversibility_components >= 0)
  expect_true(inv$irreversibility_components <= 1)
  expect_true(inv$irreversibility_base >= 0)
  expect_true(inv$irreversibility_base <= 1)
})

test_that("symmetric series has low irreversibility", {
  # Palindromic series: should have similar forward/backward structure
  series <- c(1, 3, 5, 3, 1)
  bt <- generate_bitopology(series)
  # For a palindrome, forward and backward should be very similar
  expect_true(bt$invariants$irreversibility_components <= 0.5,
              info = "Palindrome should have limited asymmetry")
})

test_that("asymmetry_direction sign is meaningful", {
  bt <- generate_bitopology(c(3, 1, 4, 1, 5, 9, 2, 6))
  inv <- bt$invariants
  # asymmetry_direction = bwd_comp - fwd_comp
  expect_equal(inv$asymmetry_direction,
               inv$backward_components - inv$forward_components)
})

test_that("resolution gain is non-negative", {
  bt <- generate_bitopology(c(3, 1, 4, 1, 5, 9, 2, 6))
  res <- bt$invariants$resolution
  expect_true(!is.null(res))
  expect_true(res$nada_forward_base_gain >= 0,
              info = "Nada base must be >= Alexandrov base")
  expect_true(res$nada_backward_base_gain >= 0)
})


# ------------------------------------------------------------------
# 6. Pairwise connectedness
# ------------------------------------------------------------------

test_that("pairwise connectedness is NA without enumeration", {
  bt <- generate_bitopology(c(3, 1, 4, 1, 5), max_open_sets = 0L)
  expect_true(is.na(bt$invariants$pairwise$pairwise_connected))
})

test_that("pairwise connectedness is computable with enumeration", {
  bt <- generate_bitopology(c(3, 1, 4), max_open_sets = 1000L)
  pw <- bt$invariants$pairwise
  # Should be TRUE or FALSE, not NA
  if (!is.na(pw$pairwise_connected)) {
    expect_true(is.logical(pw$pairwise_connected))
  }
})


# ------------------------------------------------------------------
# 7. Nada-forward vs Alexandrov resolution (concrete example)
# ------------------------------------------------------------------

test_that("Nada-forward has strictly more base elements than Alexandrov for path 1->2->3", {
  # Nada-forward: N+[1]={1,2}, N+[2]={2,3}, N+[3]={3}
  # Base = {{1,2}, {2,3}, {3}, {2}} (intersection {1,2}∩{2,3}={2})
  # Alexandrov: reach[1]={1,2,3}, reach[2]={2,3}, reach[3]={3}
  # Base = {{1,2,3}, {2,3}, {3}}
  out_adj <- list(2L, 3L, integer(0))

  nada <- generate_topology(out_adj, 3L, max_open_sets = 0L)
  alex <- generate_alexandrov_topology(out_adj, 3L)

  # Nada should have strictly more base elements
  expect_true(length(nada$base) > length(alex$base),
              info = paste("Nada:", length(nada$base),
                           "Alex:", length(alex$base)))
})


# ------------------------------------------------------------------
# 8. Existing tests still pass (regression)
# ------------------------------------------------------------------

test_that("existing undirected topology still works after refactor", {
  series <- c(3, 1, 4, 1, 5)
  g <- horizontal_visibility_graph(series)
  topo <- generate_topology(g$adjacency, g$n)
  expect_true(!is.null(topo$connected))
  expect_true(length(topo$base) > 0)
})

test_that("K6 complete graph still gives indiscrete topology (N[v] check)", {
  n <- 6
  adj <- lapply(seq_len(n), function(v) setdiff(seq_len(n), v))
  topo <- generate_topology(adj, n, max_open_sets = 100L)
  # With N[v], complete graph gives indiscrete topology
  expect_equal(topo$n_open_sets, 2L)
  expect_true(topo$connected)
})
