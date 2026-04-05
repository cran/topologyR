// ======================================================================
// alexandrov_engine.cpp — Alexandrov topology from directed visibility graph
//
// Given a directed acyclic graph (DAG) where edges go from earlier to
// later time (as produced by directed visibility graphs), computes the
// Alexandrov topology via bitset-based reachability propagation.
//
// The Alexandrov topology is defined by: open sets = upsets of the
// reachability relation. The minimal open set containing vertex v is
// reach(v) = {v} ∪ {w : v can reach w via directed paths}.
//
// For a visibility graph DAG with natural topological ordering (time
// series index order), reachability is computed in O(n * m / 64) by
// reverse-order propagation with bitset OR.
//
// The hierarchy τ_A ⊆ τ_Nada holds in general: every upset is a union
// of Nada base elements, but intersections of Nada neighborhoods can
// produce sets that are not upsets, giving the Nada topology strictly
// finer resolution.
//
// References:
//   Alexandrov, P. (1937). Diskrete Räume.
//   McCord, M.C. (1966). Singular homology groups and homotopy groups
//     of finite topological spaces.
//   Stong, R.E. (1966). Finite topological spaces.
// ======================================================================

#include "bitset_ops.h"
#include <queue>
#include <string>


// ======================================================================
// Section 1: Templatized Alexandrov pipeline
// ======================================================================

template<int W>
static Rcpp::List alexandrov_impl(
    Rcpp::List out_adjacency,
    int n_elements,
    int max_open_sets,
    bool enumerate_topology,
    bool check_connected,
    bool verify_axioms
) {
  using Ops = BitOps<W>;

  uint64_t full_set[W];
  compute_full_set(full_set, W, n_elements);
  uint64_t empty_set[W];
  Ops::zero(empty_set);

  // ------------------------------------------------------------------
  // Step 1: Compute reachability via reverse-order propagation
  //
  // For a directed visibility graph, vertex indices respect topological
  // order (edges only go from lower to higher index). Processing in
  // reverse order guarantees that reach[w] is complete when we process v.
  //
  // reach[v] = {v} ∪ ⋃_{w ∈ out(v)} reach[w]
  //
  // Cost: O(n * avg_degree * W) = O(n * m / 64)
  // ------------------------------------------------------------------
  std::vector<uint64_t> reach(static_cast<size_t>(n_elements) * W, 0);

  // Initialize self-bits
  for (int v = 0; v < n_elements; v++) {
    bs_set_bit(&reach[static_cast<size_t>(v) * W], v);
  }

  // Propagate in reverse topological order
  for (int v = n_elements - 1; v >= 0; v--) {
    Rcpp::IntegerVector out_nb = out_adjacency[v];
    uint64_t* rv = &reach[static_cast<size_t>(v) * W];
    for (int i = 0; i < out_nb.size(); i++) {
      int w = out_nb[i] - 1;  // 1-based to 0-based
      if (w >= 0 && w < n_elements) {
        const uint64_t* rw = &reach[static_cast<size_t>(w) * W];
        Ops::bor(rv, rv, rw);
      }
    }
  }

  // ------------------------------------------------------------------
  // Step 2: Deduplicate upsets to form Alexandrov base
  //
  // The minimal open set containing v is reach[v].
  // The collection {reach[v] : v ∈ V} forms the minimal basis.
  // At most n distinct sets (often fewer if vertices share reachability).
  // ------------------------------------------------------------------
  FlatVec<W> base;
  FlatHash<W> base_set(n_elements * 2);

  // Also store the reach sets as a FlatVec for easy access
  // (subbase = base for Alexandrov, since upsets are already a base)
  for (int v = 0; v < n_elements; v++) {
    const uint64_t* rv = &reach[static_cast<size_t>(v) * W];
    if (base_set.insert(rv)) {
      base.push(rv);
    }
  }

  // ------------------------------------------------------------------
  // Step 3: Connectivity via specialization preorder on the base
  //
  // Same algorithm as the Nada engine. For each element x, compute
  // S(x) = {B ∈ base : x ∈ B} as a bitset over the base indices.
  // Two elements are comparable iff S(x) ⊆ S(y) or S(y) ⊆ S(x).
  // Connected components of the comparability graph = topological
  // connected components (McCord 1966, Stong 1966).
  // ------------------------------------------------------------------
  bool connected = true;
  std::vector<std::vector<int>> components_vec;

  if (check_connected && n_elements >= 2) {
    int B = base.size();

    if (B == 0) {
      connected = true;
      components_vec.push_back(std::vector<int>());
      for (int x = 0; x < n_elements; x++)
        components_vec[0].push_back(x + 1);
    } else {
      int W_base = (B + 63) / 64;
      std::vector<uint64_t> S(static_cast<size_t>(n_elements) * W_base, 0);

      for (int b = 0; b < B; b++) {
        for (int x = 0; x < n_elements; x++) {
          if (bs_test_bit(base.get(b), x)) {
            S[static_cast<size_t>(x) * W_base + b / 64] |=
              (1ULL << (b % 64));
          }
        }
      }

      std::vector<int> comp_id(n_elements, -1);
      int n_components = 0;

      for (int start = 0; start < n_elements; start++) {
        if (comp_id[start] != -1) continue;

        std::queue<int> q;
        q.push(start);
        comp_id[start] = n_components;

        while (!q.empty()) {
          int x = q.front();
          q.pop();
          const uint64_t* Sx = &S[static_cast<size_t>(x) * W_base];

          for (int y = 0; y < n_elements; y++) {
            if (comp_id[y] != -1) continue;
            const uint64_t* Sy = &S[static_cast<size_t>(y) * W_base];

            bool xy = true, yx = true;
            for (int w = 0; w < W_base; w++) {
              if ((Sx[w] & Sy[w]) != Sx[w]) xy = false;
              if ((Sy[w] & Sx[w]) != Sy[w]) yx = false;
              if (!xy && !yx) break;
            }

            if (xy || yx) {
              comp_id[y] = n_components;
              q.push(y);
            }
          }
        }

        n_components++;
      }

      connected = (n_components == 1);

      components_vec.resize(n_components);
      for (int x = 0; x < n_elements; x++) {
        components_vec[comp_id[x]].push_back(x + 1);
      }
    }

  } else if (n_elements == 1) {
    connected = true;
    components_vec.push_back({1});
  }

  // ------------------------------------------------------------------
  // Step 4: Optionally enumerate Alexandrov topology (closure under unions)
  //
  // Same wavefront pattern as the Nada engine.
  // ------------------------------------------------------------------
  FlatVec<W> topology_vec;
  bool topo_complete = false;
  int iter_topo = 0;
  bool axioms_ok = true;
  std::string axiom_failure_msg;
  bool axioms_checked = false;

  if (enumerate_topology) {
    FlatHash<W> tau_set(base.size() * 4 + 16);

    tau_set.insert(empty_set);
    topology_vec.push(empty_set);
    tau_set.insert(full_set);
    topology_vec.push(full_set);

    bool exceeded = false;
    for (int i = 0; i < base.size() && !exceeded; i++) {
      if (tau_set.insert(base.get(i))) {
        topology_vec.push(base.get(i));
        if (tau_set.size() > max_open_sets) exceeded = true;
      }
    }

    FlatVec<W> frontier_union;
    uint64_t union_scratch[W];

    if (!exceeded) {
      for (int i = 0; i < topology_vec.size(); i++) {
        frontier_union.push(topology_vec.get(i));
      }
    }

    while (frontier_union.size() > 0 && !exceeded) {
      iter_topo++;
      FlatVec<W> new_frontier;
      int tau_size = topology_vec.size();

      for (int fi = 0; fi < frontier_union.size() && !exceeded; fi++) {
        const uint64_t* f = frontier_union.get(fi);
        for (int ti = 0; ti < tau_size && !exceeded; ti++) {
          Ops::bor(union_scratch, f, topology_vec.get(ti));
          if (Ops::is_zero(union_scratch)) continue;
          if (Ops::equal(union_scratch, full_set)) continue;

          if (tau_set.insert(union_scratch)) {
            topology_vec.push(union_scratch);
            new_frontier.push(union_scratch);

            if (tau_set.size() > max_open_sets) {
              exceeded = true;
            }
          }
        }
      }

      frontier_union = std::move(new_frontier);
    }

    topo_complete = !exceeded && (frontier_union.size() == 0);

    if (verify_axioms && topo_complete) {
      axioms_checked = true;
      int tau_n = topology_vec.size();
      uint64_t ax_scratch[W];
      for (int i = 0; i < tau_n && axioms_ok; i++) {
        for (int j = i + 1; j < tau_n && axioms_ok; j++) {
          Ops::band(ax_scratch, topology_vec.get(i), topology_vec.get(j));
          if (!tau_set.contains(ax_scratch)) {
            axioms_ok = false;
            axiom_failure_msg =
              "Intersection of two open sets not in topology. "
              "Bug in Alexandrov construction.";
          }
        }
      }
    }
  }

  // ------------------------------------------------------------------
  // Step 5: Pack results for R
  // ------------------------------------------------------------------
  // Subbase = base for Alexandrov (upsets are already intersection-closed)
  Rcpp::List subbase_list(base.size());
  for (int i = 0; i < base.size(); i++) {
    subbase_list[i] = decode_bitset(base.get(i), W, n_elements);
  }

  Rcpp::List base_list(base.size());
  for (int i = 0; i < base.size(); i++) {
    base_list[i] = decode_bitset(base.get(i), W, n_elements);
  }

  Rcpp::List comp_list(components_vec.size());
  for (size_t i = 0; i < components_vec.size(); i++) {
    comp_list[i] = Rcpp::wrap(components_vec[i]);
  }

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("subbase") = subbase_list,
    Rcpp::Named("base") = base_list,
    Rcpp::Named("base_complete") = true,  // always complete for Alexandrov
    Rcpp::Named("connected") = check_connected
      ? Rcpp::wrap(connected)
      : Rcpp::wrap(NA_LOGICAL),
    Rcpp::Named("components") = comp_list,
    Rcpp::Named("iterations_base") = 0  // no iteration needed
  );

  if (enumerate_topology) {
    Rcpp::List topo_list(topology_vec.size());
    for (int i = 0; i < topology_vec.size(); i++) {
      topo_list[i] = decode_bitset(topology_vec.get(i), W, n_elements);
    }
    result["topology"] = topo_list;
    result["n_open_sets"] = topology_vec.size();
    result["complete"] = topo_complete;
    result["iterations_topo"] = iter_topo;

    if (axioms_checked) {
      result["axioms_ok"] = axioms_ok;
      if (!axioms_ok) {
        result["axiom_failure"] = axiom_failure_msg;
      }
    }
  } else {
    result["topology"] = R_NilValue;
    result["n_open_sets"] = NA_INTEGER;
    result["complete"] = false;
    result["iterations_topo"] = 0;
  }

  return result;
}


// ======================================================================
// Section 2: Runtime-W fallback for W > 3
// ======================================================================

static Rcpp::List alexandrov_runtime(
    Rcpp::List out_adjacency, int n_elements, int W,
    int max_open_sets, bool enumerate_topology,
    bool check_connected, bool verify_axioms
) {
  std::vector<uint64_t> full_set(W), empty_set(W, 0);
  compute_full_set(full_set.data(), W, n_elements);

  // Step 1: Reachability propagation
  std::vector<uint64_t> reach(static_cast<size_t>(n_elements) * W, 0);
  for (int v = 0; v < n_elements; v++) {
    bs_set_bit(&reach[static_cast<size_t>(v) * W], v);
  }
  for (int v = n_elements - 1; v >= 0; v--) {
    Rcpp::IntegerVector out_nb = out_adjacency[v];
    uint64_t* rv = &reach[static_cast<size_t>(v) * W];
    for (int i = 0; i < out_nb.size(); i++) {
      int w = out_nb[i] - 1;
      if (w >= 0 && w < n_elements) {
        const uint64_t* rw = &reach[static_cast<size_t>(w) * W];
        rtops::bor(rv, rv, rw, W);
      }
    }
  }

  // Step 2: Deduplicate upsets
  RtVec base(W); RtHash base_set(W, n_elements * 2);
  for (int v = 0; v < n_elements; v++) {
    const uint64_t* rv = &reach[static_cast<size_t>(v) * W];
    if (base_set.insert(rv)) base.push(rv);
  }

  // Step 3: Connectivity
  bool connected = true; std::vector<std::vector<int>> components_vec;
  if (check_connected && n_elements >= 2) {
    int B = base.size();
    if (B == 0) {
      connected = true;
      components_vec.push_back(std::vector<int>());
      for (int x = 0; x < n_elements; x++) components_vec[0].push_back(x + 1);
    } else {
      int Wb = (B + 63) / 64;
      std::vector<uint64_t> S(static_cast<size_t>(n_elements) * Wb, 0);
      for (int b = 0; b < B; b++)
        for (int x = 0; x < n_elements; x++)
          if (bs_test_bit(base.get(b), x))
            S[static_cast<size_t>(x) * Wb + b / 64] |= (1ULL << (b % 64));

      std::vector<int> comp_id(n_elements, -1); int nc = 0;
      for (int s = 0; s < n_elements; s++) {
        if (comp_id[s] != -1) continue;
        std::queue<int> q; q.push(s); comp_id[s] = nc;
        while (!q.empty()) {
          int x = q.front(); q.pop();
          const uint64_t* Sx = &S[static_cast<size_t>(x) * Wb];
          for (int y = 0; y < n_elements; y++) {
            if (comp_id[y] != -1) continue;
            const uint64_t* Sy = &S[static_cast<size_t>(y) * Wb];
            bool xy = true, yx = true;
            for (int w = 0; w < Wb; w++) {
              if ((Sx[w] & Sy[w]) != Sx[w]) xy = false;
              if ((Sy[w] & Sx[w]) != Sy[w]) yx = false;
              if (!xy && !yx) break;
            }
            if (xy || yx) { comp_id[y] = nc; q.push(y); }
          }
        }
        nc++;
      }
      connected = (nc == 1);
      components_vec.resize(nc);
      for (int x = 0; x < n_elements; x++) components_vec[comp_id[x]].push_back(x + 1);
    }
  } else if (n_elements == 1) {
    connected = true; components_vec.push_back({1});
  }

  // Step 4: Optional enumeration
  RtVec topology_vec(W); bool topo_complete = false; int iter_topo = 0;
  bool axioms_ok = true; std::string axiom_fail; bool axioms_checked = false;
  if (enumerate_topology) {
    RtHash tau_set(W, base.size() * 4 + 16);
    tau_set.insert(empty_set.data()); topology_vec.push(empty_set.data());
    tau_set.insert(full_set.data()); topology_vec.push(full_set.data());
    bool exceeded = false;
    for (int i = 0; i < base.size() && !exceeded; i++) {
      if (tau_set.insert(base.get(i))) {
        topology_vec.push(base.get(i));
        if (tau_set.size() > max_open_sets) exceeded = true;
      }
    }
    RtVec fu(W); std::vector<uint64_t> uscratch(W);
    if (!exceeded) for (int i = 0; i < topology_vec.size(); i++) fu.push(topology_vec.get(i));
    while (fu.size() > 0 && !exceeded) {
      iter_topo++; RtVec nf(W); int ts = topology_vec.size();
      for (int fi = 0; fi < fu.size() && !exceeded; fi++) {
        const uint64_t* f = fu.get(fi);
        for (int ti = 0; ti < ts && !exceeded; ti++) {
          rtops::bor(uscratch.data(), f, topology_vec.get(ti), W);
          if (rtops::is_zero(uscratch.data(), W)) continue;
          if (rtops::equal(uscratch.data(), full_set.data(), W)) continue;
          if (tau_set.insert(uscratch.data())) {
            topology_vec.push(uscratch.data()); nf.push(uscratch.data());
            if (tau_set.size() > max_open_sets) exceeded = true;
          }
        }
      }
      fu = std::move(nf);
    }
    topo_complete = !exceeded && (fu.size() == 0);
    if (verify_axioms && topo_complete) {
      axioms_checked = true; int tn = topology_vec.size();
      std::vector<uint64_t> ax(W);
      for (int i = 0; i < tn && axioms_ok; i++)
        for (int j = i + 1; j < tn && axioms_ok; j++) {
          rtops::band(ax.data(), topology_vec.get(i), topology_vec.get(j), W);
          if (!tau_set.contains(ax.data())) {
            axioms_ok = false; axiom_fail = "Intersection not in topology.";
          }
        }
    }
  }

  // Step 5: Pack
  Rcpp::List sl(base.size()); for (int i = 0; i < base.size(); i++) sl[i] = decode_bitset(base.get(i), W, n_elements);
  Rcpp::List bl(base.size()); for (int i = 0; i < base.size(); i++) bl[i] = decode_bitset(base.get(i), W, n_elements);
  Rcpp::List cl(components_vec.size()); for (size_t i = 0; i < components_vec.size(); i++) cl[i] = Rcpp::wrap(components_vec[i]);
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("subbase") = sl, Rcpp::Named("base") = bl,
    Rcpp::Named("base_complete") = true,
    Rcpp::Named("connected") = check_connected ? Rcpp::wrap(connected) : Rcpp::wrap(NA_LOGICAL),
    Rcpp::Named("components") = cl, Rcpp::Named("iterations_base") = 0);
  if (enumerate_topology) {
    Rcpp::List tl(topology_vec.size()); for (int i = 0; i < topology_vec.size(); i++) tl[i] = decode_bitset(topology_vec.get(i), W, n_elements);
    result["topology"] = tl; result["n_open_sets"] = topology_vec.size();
    result["complete"] = topo_complete; result["iterations_topo"] = iter_topo;
    if (axioms_checked) { result["axioms_ok"] = axioms_ok; if (!axioms_ok) result["axiom_failure"] = axiom_fail; }
  } else {
    result["topology"] = R_NilValue; result["n_open_sets"] = NA_INTEGER;
    result["complete"] = false; result["iterations_topo"] = 0;
  }
  return result;
}


// ======================================================================
// Section 3: Runtime dispatch -> compile-time template
// ======================================================================

// [[Rcpp::export]]
Rcpp::List alexandrov_topology_engine(
    Rcpp::List out_adjacency,
    int n_elements,
    int max_open_sets = 0,
    bool enumerate_topology = false,
    bool check_connected = true,
    bool verify_axioms = false
) {
  if (n_elements < 1) Rcpp::stop("n_elements must be >= 1.");

  int W = (n_elements + 63) / 64;

  switch (W) {
    case 1: return alexandrov_impl<1>(out_adjacency, n_elements, max_open_sets,
                                      enumerate_topology, check_connected,
                                      verify_axioms);
    case 2: return alexandrov_impl<2>(out_adjacency, n_elements, max_open_sets,
                                      enumerate_topology, check_connected,
                                      verify_axioms);
    case 3: return alexandrov_impl<3>(out_adjacency, n_elements, max_open_sets,
                                      enumerate_topology, check_connected,
                                      verify_axioms);
    default:
      return alexandrov_runtime(out_adjacency, n_elements, W,
                                max_open_sets, enumerate_topology,
                                check_connected, verify_axioms);
  }
}
