// ======================================================================
// topology_engine.cpp — templatized multi-word bitset topology engine
//
// Template parameter W (number of uint64_t words per bitset) is resolved
// at compile time.  W=1 (n <= 64) produces single-instruction bitwise
// ops with zero loop overhead.  W>1 uses word-by-word loops.
//
// Architecture (Strategy 3 — base as primary object):
//   1. Subbase (neighborhoods) -> encoded as bitsets
//   2. Base (closure under finite intersections, wavefront)
//      with max_base_sets safety limit
//   3. Connectivity via specialization preorder on the base (exact,
//      polynomial, no topology enumeration needed)
//   4. Topology enumeration (optional, closure under unions, with
//      max_open_sets safety limit)
// ======================================================================

#include "bitset_ops.h"
#include <queue>
#include <string>


// ======================================================================
// Section 5: Templatized pipeline implementation
// ======================================================================

template<int W>
static Rcpp::List engine_impl(
    Rcpp::List adjacency,
    int n_elements,
    int max_open_sets,
    int max_base_sets,
    bool enumerate_topology,
    bool check_connected,
    bool verify_axioms
) {
  using Ops = BitOps<W>;

  uint64_t full_set[W];
  compute_full_set(full_set, W, n_elements);
  uint64_t empty_set[W];
  Ops::zero(empty_set);
  uint64_t scratch[W];

  // ------------------------------------------------------------------
  // Step 1: Encode subbase (neighborhoods)
  // ------------------------------------------------------------------
  FlatVec<W> subbase;
  FlatHash<W> subbase_set(n_elements * 2);

  for (int v = 0; v < n_elements; v++) {
    Rcpp::IntegerVector neighbors = adjacency[v];
    Ops::zero(scratch);
    bs_set_bit(scratch, v);  // N[v]: closed neighborhood (Nada et al. 2018)
    for (int i = 0; i < neighbors.size(); i++) {
      int idx = neighbors[i] - 1;
      if (idx >= 0 && idx < n_elements) bs_set_bit(scratch, idx);
    }
    if (subbase_set.insert(scratch)) {
      subbase.push(scratch);
    }
  }

  // ------------------------------------------------------------------
  // Step 2: Base = closure of subbase under finite intersections
  //         with max_base_sets safety limit
  // ------------------------------------------------------------------
  FlatVec<W> base;
  FlatHash<W> base_set(subbase.size() * 4 + 16);
  bool base_complete = true;

  for (int i = 0; i < subbase.size(); i++) {
    base_set.insert(subbase.get(i));
    base.push(subbase.get(i));
  }

  FlatVec<W> frontier_inter;
  for (int i = 0; i < subbase.size(); i++) {
    frontier_inter.push(subbase.get(i));
  }

  int iter_base = 0;
  int max_iter_base = n_elements * 2 + 100;
  uint64_t inter_scratch[W];
  bool base_exceeded = false;

  while (frontier_inter.size() > 0 && iter_base < max_iter_base && !base_exceeded) {
    iter_base++;
    FlatVec<W> new_frontier;
    int bsz = base.size();

    for (int fi = 0; fi < frontier_inter.size() && !base_exceeded; fi++) {
      const uint64_t* f = frontier_inter.get(fi);
      for (int bi = 0; bi < bsz && !base_exceeded; bi++) {
        Ops::band(inter_scratch, f, base.get(bi));
        if (!Ops::is_zero(inter_scratch) &&
            base_set.insert(inter_scratch)) {
          base.push(inter_scratch);
          new_frontier.push(inter_scratch);
          if (base_set.size() > max_base_sets) {
            base_exceeded = true;
            base_complete = false;
          }
        }
      }
    }

    frontier_inter = std::move(new_frontier);
  }

  if (frontier_inter.size() > 0 && !base_exceeded) {
    // max_iter_base reached without convergence
    base_complete = false;
  }

  // ------------------------------------------------------------------
  // Step 3: Connectivity via specialization preorder on the base
  // ------------------------------------------------------------------
  bool connected = true;
  std::vector<std::vector<int>> components_vec;

  if (check_connected && n_elements >= 2) {
    int B = base.size();

    if (B == 0) {
      // Empty base: topology is indiscrete ({empty, V}), which is connected.
      connected = true;
      components_vec.push_back(std::vector<int>());
      for (int x = 0; x < n_elements; x++)
        components_vec[0].push_back(x + 1);
    } else {
      // S(x) as bitsets over {0..B-1}
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

      // BFS on comparability graph using runtime-W subset check
      // (W_base is data-dependent, cannot be templatized)
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

            // Runtime subset check (W_base varies with |base|)
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
  // Step 4: Optionally enumerate topology (closure under unions)
  // ------------------------------------------------------------------
  FlatVec<W> topology_vec;
  bool topo_complete = false;
  int iter_topo = 0;
  bool axioms_ok = true;
  std::string axiom_failure_msg;
  bool axioms_checked = false;

  if (enumerate_topology) {
    // Warn if enumerating a known-disconnected space
    if (check_connected && !connected) {
      Rcpp::warning(
        "Enumerating topology of a space already determined to be "
        "disconnected (%d components). Consider using max_open_sets = 0 "
        "to skip enumeration when only connectivity is needed.",
        static_cast<int>(components_vec.size()));
    }

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
              "Bug in base construction.";
          }
        }
      }
    }
  }

  // ------------------------------------------------------------------
  // Step 5: Pack results for R
  // ------------------------------------------------------------------
  Rcpp::List subbase_list(subbase.size());
  for (int i = 0; i < subbase.size(); i++) {
    subbase_list[i] = decode_bitset(subbase.get(i), W, n_elements);
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
    Rcpp::Named("base_complete") = base_complete,
    Rcpp::Named("connected") = check_connected
      ? Rcpp::wrap(connected)
      : Rcpp::wrap(NA_LOGICAL),
    Rcpp::Named("components") = comp_list,
    Rcpp::Named("iterations_base") = iter_base
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
// Section 6: Runtime dispatch -> compile-time template
//
// W=1 (n<=64): all BitOps compile to single instructions
// W=2 (n<=128): unrolled to pairs
// W=3 (n<=192): unrolled to triples
// W>3: general runtime loop (still fast, but with loop overhead)
// ======================================================================

// Runtime-W fallback for W > 3 (uses the original non-templatized path)
static Rcpp::List engine_runtime(
    Rcpp::List adjacency, int n_elements, int W,
    int max_open_sets, int max_base_sets,
    bool enumerate_topology, bool check_connected, bool verify_axioms);

// [[Rcpp::export]]
Rcpp::List generate_topology_engine(
    Rcpp::List adjacency,
    int n_elements,
    int max_open_sets = 1000000,
    int max_base_sets = 100000,
    bool enumerate_topology = true,
    bool check_connected = true,
    bool verify_axioms = false
) {
  if (n_elements < 1) Rcpp::stop("n_elements must be >= 1.");

  int W = (n_elements + 63) / 64;

  switch (W) {
    case 1: return engine_impl<1>(adjacency, n_elements, max_open_sets,
                                  max_base_sets, enumerate_topology,
                                  check_connected, verify_axioms);
    case 2: return engine_impl<2>(adjacency, n_elements, max_open_sets,
                                  max_base_sets, enumerate_topology,
                                  check_connected, verify_axioms);
    case 3: return engine_impl<3>(adjacency, n_elements, max_open_sets,
                                  max_base_sets, enumerate_topology,
                                  check_connected, verify_axioms);
    default:
      return engine_runtime(adjacency, n_elements, W,
                            max_open_sets, max_base_sets,
                            enumerate_topology, check_connected,
                            verify_axioms);
  }
}


// ======================================================================
// Section 7: Runtime-W fallback for W > 3
//
// Uses the non-templatized operations with runtime W parameter.
// This avoids instantiating templates for every possible W.
// ======================================================================

static Rcpp::List engine_runtime(
    Rcpp::List adjacency, int n_elements, int W,
    int max_open_sets, int max_base_sets,
    bool enumerate_topology, bool check_connected, bool verify_axioms
) {
  std::vector<uint64_t> full_set(W), empty_set(W, 0), scratch(W);
  compute_full_set(full_set.data(), W, n_elements);

  // Step 1: Subbase (N[v] — closed neighborhood)
  RtVec subbase(W); RtHash subbase_set(W, n_elements * 2);
  for (int v = 0; v < n_elements; v++) {
    Rcpp::IntegerVector nb = adjacency[v];
    std::memset(scratch.data(), 0, W*8);
    bs_set_bit(scratch.data(), v);  // N[v]: closed neighborhood (Nada et al. 2018)
    for (int i = 0; i < nb.size(); i++) { int idx = nb[i]-1; if (idx>=0 && idx<n_elements) bs_set_bit(scratch.data(), idx); }
    if (subbase_set.insert(scratch.data())) subbase.push(scratch.data());
  }

  // Step 2: Base
  RtVec base(W); RtHash base_set(W, subbase.size()*4+16);
  bool base_complete = true;
  for (int i = 0; i < subbase.size(); i++) { base_set.insert(subbase.get(i)); base.push(subbase.get(i)); }
  RtVec frontier_inter(W);
  for (int i = 0; i < subbase.size(); i++) frontier_inter.push(subbase.get(i));
  int iter_base = 0, max_iter_base = n_elements*2+100;
  std::vector<uint64_t> iscratch(W);
  bool base_exceeded = false;
  while (frontier_inter.size() > 0 && iter_base < max_iter_base && !base_exceeded) {
    iter_base++; RtVec nf(W); int bsz = base.size();
    for (int fi = 0; fi < frontier_inter.size() && !base_exceeded; fi++) {
      const uint64_t* f = frontier_inter.get(fi);
      for (int bi = 0; bi < bsz && !base_exceeded; bi++) {
        rtops::band(iscratch.data(), f, base.get(bi), W);
        if (!rtops::is_zero(iscratch.data(), W) && base_set.insert(iscratch.data())) {
          base.push(iscratch.data()); nf.push(iscratch.data());
          if (base_set.size() > max_base_sets) { base_exceeded = true; base_complete = false; }
        }
      }
    }
    frontier_inter = std::move(nf);
  }
  if (frontier_inter.size() > 0 && !base_exceeded) base_complete = false;

  // Step 3: Connectivity
  bool connected = true; std::vector<std::vector<int>> components_vec;
  if (check_connected && n_elements >= 2) {
    int B = base.size();
    if (B == 0) { connected = true; components_vec.push_back(std::vector<int>()); for (int x=0;x<n_elements;x++) components_vec[0].push_back(x+1); }
    else {
      int Wb = (B+63)/64;
      std::vector<uint64_t> S(static_cast<size_t>(n_elements)*Wb, 0);
      for (int b=0;b<B;b++) for (int x=0;x<n_elements;x++) if (bs_test_bit(base.get(b),x)) S[static_cast<size_t>(x)*Wb+b/64] |= (1ULL<<(b%64));
      std::vector<int> comp_id(n_elements,-1); int nc=0;
      for (int s=0;s<n_elements;s++) { if (comp_id[s]!=-1) continue;
        std::queue<int> q; q.push(s); comp_id[s]=nc;
        while (!q.empty()) { int x=q.front(); q.pop();
          const uint64_t* Sx=&S[static_cast<size_t>(x)*Wb];
          for (int y=0;y<n_elements;y++) { if (comp_id[y]!=-1) continue;
            const uint64_t* Sy=&S[static_cast<size_t>(y)*Wb];
            bool xy=true,yx=true;
            for (int w=0;w<Wb;w++) { if ((Sx[w]&Sy[w])!=Sx[w]) xy=false; if ((Sy[w]&Sx[w])!=Sy[w]) yx=false; if (!xy&&!yx) break; }
            if (xy||yx) { comp_id[y]=nc; q.push(y); }
          }
        }
        nc++;
      }
      connected=(nc==1); components_vec.resize(nc);
      for (int x=0;x<n_elements;x++) components_vec[comp_id[x]].push_back(x+1);
    }
  } else if (n_elements==1) { connected=true; components_vec.push_back({1}); }

  // Step 4: Optional enumeration
  RtVec topology_vec(W); bool topo_complete=false; int iter_topo=0;
  bool axioms_ok=true; std::string axiom_fail; bool axioms_checked=false;
  if (enumerate_topology) {
    if (check_connected && !connected)
      Rcpp::warning("Enumerating topology of a space already determined to be disconnected (%d components).", static_cast<int>(components_vec.size()));
    RtHash tau_set(W, base.size()*4+16);
    tau_set.insert(empty_set.data()); topology_vec.push(empty_set.data());
    tau_set.insert(full_set.data()); topology_vec.push(full_set.data());
    bool exceeded=false;
    for (int i=0;i<base.size()&&!exceeded;i++) { if (tau_set.insert(base.get(i))) { topology_vec.push(base.get(i)); if (tau_set.size()>max_open_sets) exceeded=true; } }
    RtVec fu(W); std::vector<uint64_t> uscratch(W);
    if (!exceeded) for (int i=0;i<topology_vec.size();i++) fu.push(topology_vec.get(i));
    while (fu.size()>0 && !exceeded) { iter_topo++; RtVec nf(W); int ts=topology_vec.size();
      for (int fi=0;fi<fu.size()&&!exceeded;fi++) { const uint64_t* f=fu.get(fi);
        for (int ti=0;ti<ts&&!exceeded;ti++) { rtops::bor(uscratch.data(),f,topology_vec.get(ti),W);
          if (rtops::is_zero(uscratch.data(),W)) continue;
          if (rtops::equal(uscratch.data(),full_set.data(),W)) continue;
          if (tau_set.insert(uscratch.data())) { topology_vec.push(uscratch.data()); nf.push(uscratch.data());
            if (tau_set.size()>max_open_sets) exceeded=true; }
        }
      }
      fu=std::move(nf);
    }
    topo_complete=!exceeded&&(fu.size()==0);
    if (verify_axioms&&topo_complete) { axioms_checked=true; int tn=topology_vec.size();
      std::vector<uint64_t> ax(W);
      for (int i=0;i<tn&&axioms_ok;i++) for (int j=i+1;j<tn&&axioms_ok;j++) {
        rtops::band(ax.data(),topology_vec.get(i),topology_vec.get(j),W);
        if (!tau_set.contains(ax.data())) { axioms_ok=false; axiom_fail="Intersection of two open sets not in topology."; }
      }
    }
  }

  // Step 5: Pack
  Rcpp::List sl(subbase.size()); for (int i=0;i<subbase.size();i++) sl[i]=decode_bitset(subbase.get(i),W,n_elements);
  Rcpp::List bl(base.size()); for (int i=0;i<base.size();i++) bl[i]=decode_bitset(base.get(i),W,n_elements);
  Rcpp::List cl(components_vec.size()); for (size_t i=0;i<components_vec.size();i++) cl[i]=Rcpp::wrap(components_vec[i]);
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("subbase")=sl, Rcpp::Named("base")=bl,
    Rcpp::Named("base_complete")=base_complete,
    Rcpp::Named("connected")=check_connected?Rcpp::wrap(connected):Rcpp::wrap(NA_LOGICAL),
    Rcpp::Named("components")=cl, Rcpp::Named("iterations_base")=iter_base);
  if (enumerate_topology) {
    Rcpp::List tl(topology_vec.size()); for (int i=0;i<topology_vec.size();i++) tl[i]=decode_bitset(topology_vec.get(i),W,n_elements);
    result["topology"]=tl; result["n_open_sets"]=topology_vec.size();
    result["complete"]=topo_complete; result["iterations_topo"]=iter_topo;
    if (axioms_checked) { result["axioms_ok"]=axioms_ok; if (!axioms_ok) result["axiom_failure"]=axiom_fail; }
  } else { result["topology"]=R_NilValue; result["n_open_sets"]=NA_INTEGER; result["complete"]=false; result["iterations_topo"]=0; }
  return result;
}


// ======================================================================
// Section 8: Standalone exact connectivity check (complement-based)
// ======================================================================

// [[Rcpp::export]]
Rcpp::List is_connected_exact_cpp(Rcpp::List topology_sets, int n_elements) {
  if (n_elements < 1) {
    return Rcpp::List::create(
      Rcpp::Named("connected") = Rcpp::LogicalVector::create(NA_LOGICAL)
    );
  }

  int W = (n_elements + 63) / 64;
  std::vector<uint64_t> full_set(W);
  compute_full_set(full_set.data(), W, n_elements);

  RtHash tau_set(W, topology_sets.size() * 2);
  RtVec tau_vec(W);
  std::vector<uint64_t> scratch(W);

  for (int i = 0; i < topology_sets.size(); i++) {
    std::memset(scratch.data(), 0, W * 8);
    Rcpp::IntegerVector iv = topology_sets[i];
    for (int j = 0; j < iv.size(); j++) {
      int idx = iv[j] - 1;
      if (idx >= 0 && idx < n_elements) bs_set_bit(scratch.data(), idx);
    }
    if (tau_set.insert(scratch.data())) tau_vec.push(scratch.data());
  }

  bool connected = true;
  Rcpp::IntegerVector comp1, comp2;
  std::vector<uint64_t> complement(W);

  for (int i = 0; i < tau_vec.size(); i++) {
    const uint64_t* u = tau_vec.get(i);
    if (rtops::is_zero(u, W)) continue;
    if (rtops::equal(u, full_set.data(), W)) continue;

    rtops::bnot(complement.data(), u, full_set.data(), W);
    if (tau_set.contains(complement.data())) {
      connected = false;
      comp1 = decode_bitset(u, W, n_elements);
      comp2 = decode_bitset(complement.data(), W, n_elements);
      break;
    }
  }

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("connected") = connected
  );
  if (!connected) {
    result["component1"] = comp1;
    result["component2"] = comp2;
  }
  return result;
}
