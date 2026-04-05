#include <Rcpp.h>
#include <vector>
#include <stack>
#include <limits>

// ==========================================================================
// Horizontal Visibility Graph (HVG) — O(n) stack-based algorithm
//
// Two observations (t_a, y_a) and (t_b, y_b) with t_a < t_b are connected
// iff every intermediate observation (t_c, y_c) satisfies:
//   y_c < min(y_a, y_b)
//
// Algorithm: modified "next greater or equal element" using a monotone stack.
// The stack maintains a strictly decreasing sequence of y-values.
//
// Key correctness detail: when popping with <=, if ANY popped element has
// y == y[i], then an equal-valued point exists between the remaining stack
// top and i, blocking that edge (since y_c = y_i = min(y_top, y_i) violates
// the strict inequality). The `had_equal` flag handles this.
//
// Reference: Luque, B., Lacasa, L., Ballesteros, F., & Luque, J. (2009).
// Horizontal visibility graphs: exact results for random time series.
// Physical Review E, 80(4), 046103.
// ==========================================================================

// [[Rcpp::export]]
Rcpp::List hvg_cpp(Rcpp::NumericVector y) {
  int n = y.size();
  if (n < 2) {
    return Rcpp::List::create(
      Rcpp::Named("from") = Rcpp::IntegerVector(),
      Rcpp::Named("to")   = Rcpp::IntegerVector()
    );
  }

  std::vector<int> from_vec;
  std::vector<int> to_vec;
  from_vec.reserve(n * 2);
  to_vec.reserve(n * 2);

  std::stack<int> stk;

  for (int i = 0; i < n; i++) {
    bool had_equal = false;

    // Pop all elements with y <= y[i]. Each popped element j has its
    // "next greater or equal to the right" at i, meaning all values
    // between j and i are < y[j] <= y[i], so (j, i) is an HVG edge.
    while (!stk.empty() && y[stk.top()] <= y[i]) {
      int j = stk.top();
      stk.pop();
      if (y[j] == y[i]) had_equal = true;
      from_vec.push_back(j + 1);  // 1-based for R
      to_vec.push_back(i + 1);
    }

    // The remaining stack top (if any) has y > y[i]. It is connected
    // to i iff no intermediate point has y == y[i] (which would block
    // the horizontal line at height y[i] = min(y_top, y_i)).
    if (!stk.empty() && !had_equal) {
      from_vec.push_back(stk.top() + 1);
      to_vec.push_back(i + 1);
    }

    stk.push(i);
  }

  return Rcpp::List::create(
    Rcpp::Named("from") = Rcpp::wrap(from_vec),
    Rcpp::Named("to")   = Rcpp::wrap(to_vec)
  );
}


// ==========================================================================
// Natural Visibility Graph (NVG) — O(n log n) expected, O(n^2) worst case
//
// Two observations (t_a, y_a) and (t_b, y_b) with t_a < t_b are connected
// iff every intermediate observation (t_c, y_c) satisfies:
//   y_c < y_a + (y_b - y_a) * (t_c - t_a) / (t_b - t_a)
//
// Geometrically: a straight line between the two points clears all
// intermediate points (strict inequality).
//
// Algorithm: divide-and-conquer on the global maximum.
// Key property: the maximum point in [lo, hi] BLOCKS all cross-edges
// between the left and right halves (proof: the line between any left
// point l and right point r at position p has value at most
// max(y_l, y_r) <= y_p, violating the strict inequality).
// Therefore all edges are either:
//   (a) from the pivot to visible points (found by slope scan), or
//   (b) internal to the left/right halves (found recursively).
//
// The slope scan from pivot p looking left uses:
//   theta(j) = (y_p - y_j) / (p - j)
// Point j is visible from p iff theta(j) < min(theta(k)) for all k
// between j and p. This reduces to tracking the running minimum.
//
// Reference: Lacasa, L., Luque, B., Ballesteros, F., Luque, J., &
// Nuno, J.C. (2008). From time series to complex networks: The
// visibility graph. PNAS, 105(13), 4972-4975.
// ==========================================================================

namespace {

// Find index of maximum value in y[lo..hi] (inclusive).
// Returns leftmost maximum for determinism.
int find_max_idx(const double* y, int lo, int hi) {
  int idx = lo;
  for (int i = lo + 1; i <= hi; i++) {
    if (y[i] > y[idx]) {
      idx = i;
    }
  }
  return idx;
}

void nvg_recursive(const double* y, int lo, int hi,
                   std::vector<int>& from_vec, std::vector<int>& to_vec) {
  if (lo >= hi) return;

  // Base case: two adjacent points are always mutually visible
  if (hi - lo == 1) {
    from_vec.push_back(lo + 1);
    to_vec.push_back(hi + 1);
    return;
  }

  // Find pivot (global maximum in [lo, hi])
  int p = find_max_idx(y, lo, hi);

  // Left scan: from pivot looking left, find all visible points.
  // theta(j) = (y_p - y_j) / (p - j) represents the "drop rate".
  // A point j is visible iff its drop rate is less than all points
  // between j and p (i.e., less than the running minimum).
  {
    double min_theta = std::numeric_limits<double>::infinity();
    for (int j = p - 1; j >= lo; j--) {
      double theta = (y[p] - y[j]) / static_cast<double>(p - j);
      if (theta < min_theta) {
        from_vec.push_back(j + 1);
        to_vec.push_back(p + 1);
        min_theta = theta;
      }
    }
  }

  // Right scan: symmetric to left scan
  {
    double min_theta = std::numeric_limits<double>::infinity();
    for (int j = p + 1; j <= hi; j++) {
      double theta = (y[p] - y[j]) / static_cast<double>(j - p);
      if (theta < min_theta) {
        from_vec.push_back(p + 1);
        to_vec.push_back(j + 1);
        min_theta = theta;
      }
    }
  }

  // Recurse on left and right halves (no cross-edges exist)
  nvg_recursive(y, lo, p - 1, from_vec, to_vec);
  nvg_recursive(y, p + 1, hi, from_vec, to_vec);
}

} // anonymous namespace


// [[Rcpp::export]]
Rcpp::List nvg_cpp(Rcpp::NumericVector y) {
  int n = y.size();
  if (n < 2) {
    return Rcpp::List::create(
      Rcpp::Named("from") = Rcpp::IntegerVector(),
      Rcpp::Named("to")   = Rcpp::IntegerVector()
    );
  }

  std::vector<int> from_vec;
  std::vector<int> to_vec;
  // Expected number of edges: O(n log n) for NVG on random series
  from_vec.reserve(n * 3);
  to_vec.reserve(n * 3);

  nvg_recursive(y.begin(), 0, n - 1, from_vec, to_vec);

  return Rcpp::List::create(
    Rcpp::Named("from") = Rcpp::wrap(from_vec),
    Rcpp::Named("to")   = Rcpp::wrap(to_vec)
  );
}
