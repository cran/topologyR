#' Generate the Alexandrov topology from a directed visibility graph
#'
#' @description
#' Computes the Alexandrov topology induced by a directed acyclic graph (DAG)
#' via bitset-based reachability propagation. The open sets are the upsets of
#' the reachability relation: vertex \eqn{v}'s minimal open set is the set of
#' all vertices reachable from \eqn{v} via directed paths (including \eqn{v}
#' itself).
#'
#' The Alexandrov topology is a \strong{lower bound} on the resolution of
#' the Nada et al. (2018) topology: \eqn{\tau_A \subseteq \tau^+_{\mathrm{Nada}}}.
#' The difference in base sizes quantifies how much additional structure
#' the intersection/union closure of Nada et al. captures beyond the pure
#' order structure.
#'
#' @param out_adjacency A list of integer vectors, where element \code{i}
#'   contains the 1-based indices of all \strong{forward} (later-in-time)
#'   neighbors of node \code{i}. Must represent a DAG where all edges go
#'   from lower to higher index. Typically obtained from
#'   \code{\link{horizontal_visibility_graph}} or
#'   \code{\link{natural_visibility_graph}} with \code{directed = TRUE}.
#' @param n_elements Integer. Total number of elements in the ground set V.
#' @param max_open_sets Integer. Maximum number of open sets to enumerate
#'   (default: 0, meaning skip enumeration). The Alexandrov topology can have
#'   up to \eqn{2^n} open sets, so enumeration is only feasible for small
#'   \eqn{n}. Connectivity is always computed from the base without
#'   enumeration.
#' @param check_connected Logical. Whether to compute exact topological
#'   connectivity via the specialization preorder (default: \code{TRUE}).
#' @param verify_axioms Logical. Whether to verify topology axioms after
#'   enumeration (default: \code{FALSE}). Only meaningful when the topology
#'   is fully enumerated.
#' @return A \code{list} with the same structure as
#'   \code{\link{generate_topology}}:
#'   \describe{
#'     \item{subbase}{List of integer vectors. The upsets (= base for Alexandrov).}
#'     \item{base}{List of integer vectors. Same as subbase (upsets are already
#'       closed under intersection).}
#'     \item{base_complete}{Always \code{TRUE} (no iteration needed).}
#'     \item{connected}{Logical. Exact topological connectivity.}
#'     \item{components}{List of integer vectors. Exact connected components.}
#'     \item{topology}{List of integer vectors, or \code{NULL}. Open sets if
#'       enumerated.}
#'     \item{n_open_sets}{Integer, or \code{NA}. Number of open sets if enumerated.}
#'     \item{complete}{Logical. Whether enumeration finished.}
#'   }
#'
#' @details
#' The algorithm computes reachability in \eqn{O(n \cdot m / 64)} time by
#' processing vertices in reverse topological order (= reverse time order
#' for visibility graphs) and propagating reachability sets as multi-word
#' bitsets with OR operations. This reuses the same bitset infrastructure
#' as the Nada engine with compile-time template dispatch for
#' \eqn{n \leq 192}.
#'
#' The Alexandrov topology on a finite set equipped with a preorder is
#' canonical: there is a bijection between preorders and Alexandrov
#' topologies (Alexandrov, 1937). For directed visibility graphs, the
#' reachability preorder captures the temporal ordering structure of the
#' time series.
#'
#' @references
#' Alexandrov, P. (1937). Diskrete Raume. \emph{Matematicheskii Sbornik},
#' 2(44), 501-519.
#'
#' @seealso \code{\link{generate_topology}} for the Nada et al. topology,
#'   \code{\link{generate_bitopology}} for the complete bitopological analysis.
#'
#' @examples
#' series <- c(3, 1, 4, 1, 5)
#' g <- horizontal_visibility_graph(series, directed = TRUE)
#' alex <- generate_alexandrov_topology(g$out_adjacency, g$n)
#' alex$connected
#' length(alex$base)
#'
#' @export
generate_alexandrov_topology <- function(out_adjacency, n_elements,
                                         max_open_sets = 0L,
                                         check_connected = TRUE,
                                         verify_axioms = FALSE) {
  if (!is.list(out_adjacency)) {
    stop("'out_adjacency' must be a list of integer vectors.", call. = FALSE)
  }
  n_elements <- as.integer(n_elements)
  if (n_elements < 1L) {
    stop("'n_elements' must be >= 1.", call. = FALSE)
  }
  if (length(out_adjacency) != n_elements) {
    stop("Length of 'out_adjacency' must equal 'n_elements'.", call. = FALSE)
  }
  max_open_sets <- as.integer(max_open_sets)
  enumerate <- max_open_sets > 0L

  alexandrov_topology_engine(
    out_adjacency = out_adjacency,
    n_elements = n_elements,
    max_open_sets = max_open_sets,
    enumerate_topology = enumerate,
    check_connected = check_connected,
    verify_axioms = verify_axioms
  )
}


#' Generate a bitopological analysis from a time series
#'
#' @description
#' Performs a complete bitopological analysis of a time series by constructing
#' directed visibility graphs and computing four topological spaces:
#' \enumerate{
#'   \item \strong{Undirected Nada topology} \eqn{\tau}: from the symmetric
#'     (undirected) visibility graph (standard Nada et al. 2018).
#'   \item \strong{Forward Nada topology} \eqn{\tau^+}: from closed forward
#'     neighborhoods \eqn{N^+[v] = \{v\} \cup \{w : v \to w\}}.
#'   \item \strong{Backward Nada topology} \eqn{\tau^-}: from closed backward
#'     neighborhoods \eqn{N^-[v] = \{v\} \cup \{w : w \to v\}}.
#'   \item \strong{Alexandrov topology} \eqn{\tau_A}: from reachability upsets
#'     of the directed graph.
#' }
#'
#' The pair \eqn{(X, \tau^+, \tau^-)} forms a \strong{bitopological space}
#' (Kelly, 1963). The divergence between \eqn{\tau^+} and \eqn{\tau^-}
#' quantifies temporal irreversibility: in a reversible process,
#' \eqn{\tau^+ \cong \tau^-}; in an irreversible process (e.g., asymmetric
#' business cycles), they diverge.
#'
#' @param series Numeric vector representing the time series.
#' @param graph_type Character. Either \code{"hvg"} (Horizontal Visibility
#'   Graph) or \code{"nvg"} (Natural Visibility Graph). Default: \code{"hvg"}.
#' @param max_open_sets Integer. Maximum open sets to enumerate per topology
#'   (default: 0, skip enumeration). Required for pairwise connectedness check.
#' @param max_base_sets Integer. Maximum base elements during intersection
#'   closure for the Nada topologies (default: 100,000).
#' @param alexandrov Logical. Whether to also compute the Alexandrov topology
#'   for resolution comparison (default: \code{TRUE}).
#' @return A \code{list} with:
#'   \describe{
#'     \item{graph}{The directed visibility graph (output of
#'       \code{\link{horizontal_visibility_graph}} or
#'       \code{\link{natural_visibility_graph}} with \code{directed = TRUE}).}
#'     \item{undirected}{The undirected Nada topology
#'       (output of \code{\link{generate_topology}}).}
#'     \item{forward}{The forward Nada topology \eqn{\tau^+}.}
#'     \item{backward}{The backward Nada topology \eqn{\tau^-}.}
#'     \item{alexandrov}{The Alexandrov topology \eqn{\tau_A} (if requested).}
#'     \item{invariants}{Bitopological invariants (output of
#'       \code{\link{bitopology_invariants}}).}
#'   }
#'
#' @details
#' The existing undirected engine is called three times (undirected, forward,
#' backward) and the Alexandrov engine once. No new C++ code is needed for the
#' forward and backward topologies: passing directed adjacency lists to the
#' existing \code{\link{generate_topology}} produces the correct closed
#' neighborhoods \eqn{N^+[v]} or \eqn{N^-[v]} automatically (the engine
#' always adds the self-loop \eqn{v \in N[v]}).
#'
#' @references
#' Kelly, J. C. (1963). Bitopological spaces. \emph{Proceedings of the
#' London Mathematical Society}, 3(1), 71-89.
#'
#' Lacasa, L., & Toral, R. (2010). Description of stochastic and chaotic
#' series using visibility graphs. \emph{Physical Review E}, 82(3), 036120.
#'
#' @seealso \code{\link{bitopology_invariants}} for computing invariants
#'   from pre-computed topologies.
#'
#' @examples
#' series <- c(3, 1, 4, 1, 5, 9, 2, 6)
#' bt <- generate_bitopology(series, graph_type = "hvg")
#' bt$invariants$forward_components
#' bt$invariants$backward_components
#' bt$invariants$irreversibility_components
#'
#' @export
generate_bitopology <- function(series,
                                graph_type = c("hvg", "nvg"),
                                max_open_sets = 0L,
                                max_base_sets = 100000L,
                                alexandrov = TRUE) {
  graph_type <- match.arg(graph_type)

  if (!is.numeric(series)) {
    stop("'series' must be a numeric vector.", call. = FALSE)
  }
  if (anyNA(series)) {
    stop("'series' must not contain NA values.", call. = FALSE)
  }
  if (!all(is.finite(series))) {
    stop("'series' must not contain Inf or NaN values.", call. = FALSE)
  }
  if (length(series) < 2L) {
    stop("'series' must have at least 2 observations.", call. = FALSE)
  }

  # Build directed visibility graph
  if (graph_type == "hvg") {
    g <- horizontal_visibility_graph(series, directed = TRUE)
  } else {
    g <- natural_visibility_graph(series, directed = TRUE)
  }

  n <- g$n

  # Undirected topology (standard Nada et al.)
  tau_undirected <- generate_topology(g$adjacency, n,
                                      max_open_sets = max_open_sets,
                                      max_base_sets = max_base_sets)

  # Forward topology τ+ (from N+[v])
  tau_forward <- generate_topology(g$out_adjacency, n,
                                   max_open_sets = max_open_sets,
                                   max_base_sets = max_base_sets)

  # Backward topology τ- (from N-[v])
  tau_backward <- generate_topology(g$in_adjacency, n,
                                    max_open_sets = max_open_sets,
                                    max_base_sets = max_base_sets)

  # Alexandrov topology (optional)
  tau_alexandrov <- NULL
  if (alexandrov) {
    tau_alexandrov <- generate_alexandrov_topology(
      g$out_adjacency, n,
      max_open_sets = max_open_sets
    )
  }

  # Compute invariants
  inv <- bitopology_invariants(tau_forward, tau_backward, n,
                               tau_undirected, tau_alexandrov)

  list(
    graph = g,
    undirected = tau_undirected,
    forward = tau_forward,
    backward = tau_backward,
    alexandrov = tau_alexandrov,
    invariants = inv
  )
}


#' Compute bitopological invariants from forward and backward topologies
#'
#' @description
#' Given the forward topology \eqn{\tau^+} and backward topology \eqn{\tau^-}
#' (as returned by \code{\link{generate_topology}}), computes bitopological
#' invariants that quantify temporal irreversibility and structural asymmetry.
#'
#' @param forward The forward Nada topology \eqn{\tau^+} (output of
#'   \code{\link{generate_topology}} called with forward adjacency).
#' @param backward The backward Nada topology \eqn{\tau^-} (output of
#'   \code{\link{generate_topology}} called with backward adjacency).
#' @param n_elements Integer. Total number of elements in the ground set V.
#' @param undirected The undirected Nada topology (optional, for comparison).
#' @param alexandrov The Alexandrov topology (optional, for resolution
#'   comparison).
#' @return A \code{list} with:
#'   \describe{
#'     \item{forward_components}{Integer. Number of connected components in
#'       \eqn{\tau^+}.}
#'     \item{backward_components}{Integer. Number of connected components in
#'       \eqn{\tau^-}.}
#'     \item{forward_connected}{Logical. Whether \eqn{\tau^+} is connected.}
#'     \item{backward_connected}{Logical. Whether \eqn{\tau^-} is connected.}
#'     \item{forward_base_size}{Integer. Number of base elements in
#'       \eqn{\tau^+}.}
#'     \item{backward_base_size}{Integer. Number of base elements in
#'       \eqn{\tau^-}.}
#'     \item{irreversibility_components}{Numeric in the range 0 to 1. Normalized
#'       asymmetry of component counts:
#'       \eqn{|C^+ - C^-| / \max(C^+, C^-)}.
#'       Zero means equal component counts (necessary for reversibility).}
#'     \item{irreversibility_base}{Numeric in the range 0 to 1. Normalized asymmetry
#'       of base sizes.}
#'     \item{asymmetry_direction}{Integer. \eqn{C^- - C^+}. Positive means
#'       the forward topology is more connected (fewer components) than the
#'       backward topology, consistent with gradual expansions (forward
#'       visibility) and abrupt contractions (backward obstruction).}
#'     \item{pairwise}{A list with pairwise connectedness results. Contains
#'       \code{pairwise_connected} (logical or \code{NA}),
#'       \code{clopen_forward} (integer vector, the forward-open / backward-
#'       closed set if disconnected), and \code{clopen_backward} (its
#'       complement). Requires both topologies to be fully enumerated
#'       (\code{max_open_sets > 0} and \code{complete = TRUE}).}
#'     \item{resolution}{A list comparing the Nada and Alexandrov topologies
#'       (if \code{alexandrov} is provided). Contains base size comparisons
#'       and relative resolution gains.}
#'     \item{undirected_info}{A list with undirected topology summary (if
#'       \code{undirected} is provided).}
#'   }
#'
#' @details
#' \strong{Pairwise connectedness} (Kelly, 1963): the bitopological space
#' \eqn{(X, \tau^+, \tau^-)} is pairwise disconnected if there exists a
#' proper non-empty subset \eqn{A \subset X} that is simultaneously
#' \eqn{\tau^+}-open and \eqn{\tau^-}-closed (its complement is
#' \eqn{\tau^-}-open). This is checked by iterating over all
#' \eqn{\tau^+} open sets and testing whether their complement is
#' \eqn{\tau^-}-open.
#'
#' \strong{Irreversibility prediction} (Lacasa & Toral, 2010): for a time
#' series with asymmetric dynamics (e.g., gradual expansions and abrupt
#' recessions), the forward topology \eqn{\tau^+} should be more connected
#' than \eqn{\tau^-}, because forward visibility is less obstructed during
#' gradual rises than backward visibility after abrupt drops.
#'
#' \strong{Resolution gain}: the difference
#' \eqn{|\mathcal{B}_{\mathrm{Nada}}| - |\mathcal{B}_A|} quantifies how
#' much additional structure the Nada intersection/union closure captures
#' beyond the pure Alexandrov order structure. A large gain means the
#' Nada pipeline is extracting genuinely new topological information.
#'
#' @references
#' Kelly, J. C. (1963). Bitopological spaces. \emph{Proceedings of the
#' London Mathematical Society}, 3(1), 71-89.
#'
#' @seealso \code{\link{generate_bitopology}} for the complete pipeline.
#'
#' @examples
#' series <- c(3, 1, 4, 1, 5, 9, 2, 6)
#' g <- horizontal_visibility_graph(series, directed = TRUE)
#' tf <- generate_topology(g$out_adjacency, g$n, max_open_sets = 0L)
#' tb <- generate_topology(g$in_adjacency, g$n, max_open_sets = 0L)
#' inv <- bitopology_invariants(tf, tb, g$n)
#' inv$asymmetry_direction
#'
#' @export
bitopology_invariants <- function(forward, backward, n_elements,
                                  undirected = NULL, alexandrov = NULL) {
  if (!is.list(forward) || is.null(forward$components)) {
    stop("'forward' must be a topology result from generate_topology().",
         call. = FALSE)
  }
  if (!is.list(backward) || is.null(backward$components)) {
    stop("'backward' must be a topology result from generate_topology().",
         call. = FALSE)
  }
  n_elements <- as.integer(n_elements)

  # --- Component analysis ---
  fwd_n_comp <- length(forward$components)
  bwd_n_comp <- length(backward$components)
  fwd_connected <- isTRUE(forward$connected)
  bwd_connected <- isTRUE(backward$connected)

  # --- Base analysis ---
  fwd_base_size <- length(forward$base)
  bwd_base_size <- length(backward$base)

  # --- Irreversibility indices ---
  max_comp <- max(fwd_n_comp, bwd_n_comp)
  irreversibility_components <- if (max_comp > 0L) {
    abs(fwd_n_comp - bwd_n_comp) / max_comp
  } else {
    0
  }

  max_base <- max(fwd_base_size, bwd_base_size)
  irreversibility_base <- if (max_base > 0L) {
    abs(fwd_base_size - bwd_base_size) / max_base
  } else {
    0
  }

  # Positive = forward more connected (fewer components)
  asymmetry_direction <- bwd_n_comp - fwd_n_comp

  # --- Pairwise connectedness (Kelly, 1963) ---
  pairwise <- list(
    pairwise_connected = NA,
    clopen_forward = NULL,
    clopen_backward = NULL
  )

  fwd_enumerated <- !is.null(forward$topology) && isTRUE(forward$complete)
  bwd_enumerated <- !is.null(backward$topology) && isTRUE(backward$complete)

  if (fwd_enumerated && bwd_enumerated) {
    # Encode τ- open sets as character keys for O(1) lookup
    bwd_env <- new.env(hash = TRUE, parent = emptyenv())
    for (s in backward$topology) {
      if (length(s) == 0L || length(s) == n_elements) next
      key <- paste(sort(s), collapse = ",")
      bwd_env[[key]] <- TRUE
    }

    all_elements <- seq_len(n_elements)
    found_clopen <- FALSE

    for (U in forward$topology) {
      if (length(U) == 0L || length(U) == n_elements) next
      complement <- setdiff(all_elements, U)
      comp_key <- paste(sort(complement), collapse = ",")
      if (!is.null(bwd_env[[comp_key]])) {
        pairwise$pairwise_connected <- FALSE
        pairwise$clopen_forward <- sort(U)
        pairwise$clopen_backward <- sort(complement)
        found_clopen <- TRUE
        break
      }
    }

    if (!found_clopen) {
      pairwise$pairwise_connected <- TRUE
    }
  }

  # --- Resolution comparison with Alexandrov ---
  resolution <- NULL
  if (!is.null(alexandrov) && is.list(alexandrov)) {
    alex_base_size <- length(alexandrov$base)
    alex_n_comp <- length(alexandrov$components)

    resolution <- list(
      alexandrov_base_size = alex_base_size,
      alexandrov_components = alex_n_comp,
      alexandrov_connected = isTRUE(alexandrov$connected),
      nada_forward_base_gain = fwd_base_size - alex_base_size,
      nada_backward_base_gain = bwd_base_size - alex_base_size,
      nada_forward_relative_gain = (fwd_base_size - alex_base_size) /
        max(alex_base_size, 1L),
      nada_backward_relative_gain = (bwd_base_size - alex_base_size) /
        max(alex_base_size, 1L)
    )
  }

  # --- Undirected comparison ---
  undirected_info <- NULL
  if (!is.null(undirected) && is.list(undirected)) {
    undirected_info <- list(
      undirected_components = length(undirected$components),
      undirected_connected = isTRUE(undirected$connected),
      undirected_base_size = length(undirected$base)
    )
  }

  list(
    forward_components = fwd_n_comp,
    backward_components = bwd_n_comp,
    forward_connected = fwd_connected,
    backward_connected = bwd_connected,
    forward_base_size = fwd_base_size,
    backward_base_size = bwd_base_size,
    irreversibility_components = irreversibility_components,
    irreversibility_base = irreversibility_base,
    asymmetry_direction = asymmetry_direction,
    pairwise = pairwise,
    resolution = resolution,
    undirected_info = undirected_info
  )
}
