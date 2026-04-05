#' Generate a topology from a graph's neighborhood structure
#'
#' @description
#' Given a graph (represented by its adjacency list of neighborhoods), generates
#' the induced topology following Nada et al. (2018). The procedure is:
#' \enumerate{
#'   \item The neighborhoods \eqn{N(v)} form the \strong{subbase}.
#'   \item The base is the closure of the subbase under finite intersections
#'     (wavefront propagation in C++).
#'   \item Topological connectivity is determined \strong{exactly} via the
#'     specialization preorder on the base (Alexandrov theorem for finite
#'     spaces), without requiring topology enumeration.
#'   \item Optionally, the full topology (closure under arbitrary unions) is
#'     enumerated, subject to a \code{max_open_sets} limit.
#' }
#'
#' The connectivity result is \strong{exact}: the connected components returned
#' are the true topological connected components, not approximations.
#'
#' @param adjacency A list of integer vectors, where element \code{i}
#'   contains the 1-based indices of all neighbors of node \code{i}.
#'   Typically obtained from \code{\link{horizontal_visibility_graph}} or
#'   \code{\link{natural_visibility_graph}}.
#' @param n_elements Integer. Total number of elements in the ground set V.
#' @param max_open_sets Integer. Maximum number of open sets to enumerate
#'   (default: 1,000,000). If the topology exceeds this limit, enumeration
#'   stops and \code{complete = FALSE}. Connectivity is unaffected (computed
#'   from the base, not the topology). Set to 0 to skip enumeration entirely.
#' @param max_base_sets Integer. Maximum number of base elements during
#'   intersection closure (default: 100,000). If exceeded, closure stops
#'   early and \code{base_complete = FALSE}. Connectivity results remain
#'   valid but may be approximate if the base is incomplete.
#' @param check_connected Logical. Whether to compute exact topological
#'   connectivity via the specialization preorder (default: \code{TRUE}).
#' @param verify_axioms Logical. Whether to verify topology axioms after
#'   enumeration (default: \code{FALSE}). Only meaningful when the topology
#'   is fully enumerated (\code{complete = TRUE}).
#' @return A \code{list} with:
#'   \describe{
#'     \item{subbase}{List of integer vectors. The neighborhoods used as subbase.}
#'     \item{base}{List of integer vectors. The base (closure under intersections).}
#'     \item{base_complete}{Logical. Whether the base closure finished within
#'       the \code{max_base_sets} limit.}
#'     \item{connected}{Logical. Whether the topological space is connected.
#'       This is an \strong{exact} result based on the specialization preorder.
#'       \code{NA} if \code{check_connected = FALSE}.}
#'     \item{components}{List of integer vectors. The exact topological connected
#'       components. Each component is a vector of 1-based element indices.}
#'     \item{topology}{List of integer vectors, or \code{NULL}. The open sets,
#'       if enumeration was requested and feasible.}
#'     \item{n_open_sets}{Integer, or \code{NA}. Number of open sets if enumerated.}
#'     \item{complete}{Logical. Whether the topology enumeration finished
#'       within the \code{max_open_sets} limit.}
#'     \item{iterations_base}{Integer. Wavefront iterations for base closure.}
#'     \item{iterations_topo}{Integer. Wavefront iterations for union closure.}
#'     \item{axioms_ok}{Logical (only if \code{verify_axioms = TRUE} and
#'       \code{complete = TRUE}).}
#'     \item{axiom_failure}{Character (only if axiom verification fails).}
#'   }
#'
#' @references
#' Nada, S., El Atik, A. E. F., & Atef, M. (2018). New types of topological
#' structures via graphs. \emph{Mathematical Methods in the Applied Sciences},
#' 41(15), 5801-5810. \doi{10.1002/mma.4726}
#'
#' @examples
#' # From a small visibility graph
#' series <- c(3, 1, 4, 1, 5)
#' g <- horizontal_visibility_graph(series)
#' topo <- generate_topology(g$adjacency, g$n)
#' topo$connected
#' length(topo$components)
#'
#' @export
generate_topology <- function(adjacency, n_elements,
                              max_open_sets = 1000000L,
                              max_base_sets = 100000L,
                              check_connected = TRUE,
                              verify_axioms = FALSE) {
  if (!is.list(adjacency)) {
    stop("'adjacency' must be a list of integer vectors.", call. = FALSE)
  }
  n_elements <- as.integer(n_elements)
  if (n_elements < 1L) {
    stop("'n_elements' must be >= 1.", call. = FALSE)
  }
  if (length(adjacency) != n_elements) {
    stop("Length of 'adjacency' must equal 'n_elements'.", call. = FALSE)
  }
  max_open_sets <- as.integer(max_open_sets)
  max_base_sets <- as.integer(max_base_sets)
  enumerate <- max_open_sets > 0L

  result <- generate_topology_engine(
    adjacency = adjacency,
    n_elements = n_elements,
    max_open_sets = max_open_sets,
    max_base_sets = max_base_sets,
    enumerate_topology = enumerate,
    check_connected = check_connected,
    verify_axioms = verify_axioms
  )

  result
}


#' Check topological connectivity of an existing topology
#'
#' @description
#' Checks whether a topological space \eqn{(V, \tau)} is connected by verifying
#' the exact definition: no proper non-empty open set has an open complement.
#' Uses multi-word bitset representation and hash lookup.
#'
#' @param topology A list of integer vectors representing the open sets.
#' @param n_elements Integer. Total number of elements in V.
#' @return A \code{list} with:
#'   \describe{
#'     \item{connected}{Logical. \code{TRUE} if the space is topologically
#'       connected.}
#'     \item{component1}{Integer vector (only if disconnected). One component.}
#'     \item{component2}{Integer vector (only if disconnected). The complement.}
#'   }
#'
#' @examples
#' # A connected topology
#' tau <- list(integer(0), 1:3, c(1L, 2L), c(2L, 3L))
#' is_topology_connected_exact(tau, 3L)
#'
#' # A disconnected topology
#' tau2 <- list(integer(0), 1:4, c(1L, 2L), c(3L, 4L))
#' is_topology_connected_exact(tau2, 4L)
#'
#' @export
is_topology_connected_exact <- function(topology, n_elements) {
  if (!is.list(topology)) {
    stop("'topology' must be a list of integer vectors.", call. = FALSE)
  }
  n_elements <- as.integer(n_elements)
  if (n_elements < 1L) {
    stop("'n_elements' must be >= 1.", call. = FALSE)
  }

  topology <- lapply(topology, as.integer)
  is_connected_exact_cpp(topology, n_elements)
}
