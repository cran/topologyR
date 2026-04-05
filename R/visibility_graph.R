#' Construct the Horizontal Visibility Graph of a time series
#'
#' @description
#' Constructs the Horizontal Visibility Graph (HVG) from a time series.
#' Two observations \eqn{(t_a, y_a)} and \eqn{(t_b, y_b)} with \eqn{t_a < t_b}
#' are connected if and only if every intermediate observation \eqn{(t_c, y_c)}
#' satisfies \eqn{y_c < \min(y_a, y_b)}.
#'
#' The HVG is parameter-free and captures structural properties of the series
#' dynamics. It is computed in \eqn{O(n)} time using a stack-based algorithm.
#'
#' @param series Numeric vector representing the time series. The position
#'   index is used as the time coordinate.
#' @param directed Logical. If \code{TRUE}, also returns directed adjacency
#'   lists (\code{out_adjacency} and \code{in_adjacency}) for constructing
#'   forward and backward Nada topologies or Alexandrov topologies. Edges
#'   are directed from earlier to later time (default: \code{FALSE}).
#' @return A \code{list} with:
#'   \describe{
#'     \item{edges}{A \code{data.frame} with columns \code{from} and \code{to}
#'       (integer, 1-based indices) representing undirected edges.}
#'     \item{n}{Integer. Number of observations in the series.}
#'     \item{n_edges}{Integer. Number of edges in the graph.}
#'     \item{adjacency}{A list of integer vectors, where element \code{i}
#'       contains the indices of all neighbors of node \code{i}. This is
#'       the neighborhood structure used as subbase for topology induction.}
#'     \item{out_adjacency}{(Only if \code{directed = TRUE}.) A list of integer
#'       vectors, where element \code{i} contains the indices of all forward
#'       (later-in-time) neighbors. Used as subbase for the forward Nada
#'       topology \eqn{\tau^+} or as input for the Alexandrov topology.}
#'     \item{in_adjacency}{(Only if \code{directed = TRUE}.) A list of integer
#'       vectors, where element \code{i} contains the indices of all backward
#'       (earlier-in-time) neighbors. Used as subbase for the backward Nada
#'       topology \eqn{\tau^-}.}
#'   }
#'
#' @details
#' The algorithm uses a monotone stack with \eqn{O(n)} time and space
#' complexity. Each element is pushed and popped at most once.
#'
#' The adjacency list returned as \code{adjacency} provides the neighborhood
#' structure \eqn{N(v)} for each vertex, which serves as the subbase for
#' inducing a topology following Nada et al. (2018).
#'
#' When \code{directed = TRUE}, the function additionally returns directed
#' adjacency lists. By construction, all edges satisfy \code{from < to}
#' (earlier time to later time), so the directed graph is a DAG with the
#' time index as topological order.
#'
#' @references
#' Luque, B., Lacasa, L., Ballesteros, F., & Luque, J. (2009). Horizontal
#' visibility graphs: exact results for random time series.
#' \emph{Physical Review E}, 80(4), 046103.
#' \doi{10.1103/PhysRevE.80.046103}
#'
#' @examples
#' series <- c(3, 1, 4, 1, 5, 9, 2, 6)
#' g <- horizontal_visibility_graph(series)
#' g$n_edges
#' head(g$edges)
#'
#' # Directed graph for bitopological analysis
#' gd <- horizontal_visibility_graph(series, directed = TRUE)
#' gd$out_adjacency[[1]]  # forward neighbors of first observation
#'
#' @export
horizontal_visibility_graph <- function(series, directed = FALSE) {
  if (!is.numeric(series)) {
    stop("'series' must be a numeric vector.", call. = FALSE)
  }
  if (anyNA(series)) {
    stop("'series' must not contain NA values.", call. = FALSE)
  }
  if (!all(is.finite(series))) {
    stop("'series' must not contain Inf or NaN values.", call. = FALSE)
  }
  n <- length(series)
  if (n < 2L) {
    result <- list(
      edges = data.frame(from = integer(0), to = integer(0)),
      n = n,
      n_edges = 0L,
      adjacency = lapply(seq_len(n), function(i) integer(0))
    )
    if (directed) {
      result$out_adjacency <- lapply(seq_len(n), function(i) integer(0))
      result$in_adjacency <- lapply(seq_len(n), function(i) integer(0))
    }
    return(result)
  }

  raw <- hvg_cpp(as.double(series))

  edges <- data.frame(from = raw$from, to = raw$to)

  # Build undirected adjacency list
  adj <- vector("list", n)
  for (i in seq_len(n)) adj[[i]] <- integer(0)
  for (k in seq_along(raw$from)) {
    i <- raw$from[k]
    j <- raw$to[k]
    adj[[i]] <- c(adj[[i]], j)
    adj[[j]] <- c(adj[[j]], i)
  }

  result <- list(
    edges = edges,
    n = n,
    n_edges = nrow(edges),
    adjacency = adj
  )

  if (directed) {
    # By construction, raw$from < raw$to (earlier -> later time)
    out_adj <- vector("list", n)
    in_adj <- vector("list", n)
    for (i in seq_len(n)) {
      out_adj[[i]] <- integer(0)
      in_adj[[i]] <- integer(0)
    }
    for (k in seq_along(raw$from)) {
      i <- raw$from[k]  # earlier time index
      j <- raw$to[k]    # later time index
      out_adj[[i]] <- c(out_adj[[i]], j)
      in_adj[[j]] <- c(in_adj[[j]], i)
    }
    result$out_adjacency <- out_adj
    result$in_adjacency <- in_adj
  }

  result
}


#' Construct the Natural Visibility Graph of a time series
#'
#' @description
#' Constructs the Natural Visibility Graph (NVG) from a time series.
#' Two observations \eqn{(t_a, y_a)} and \eqn{(t_b, y_b)} with \eqn{t_a < t_b}
#' are connected if and only if every intermediate observation \eqn{(t_c, y_c)}
#' satisfies:
#' \deqn{y_c < y_a + (y_b - y_a) \cdot \frac{t_c - t_a}{t_b - t_a}}
#'
#' The NVG is parameter-free and richer in structure than the HVG (more edges),
#' preserving dynamical properties such as periodicity, chaoticity, and
#' long-range correlations. Computed in \eqn{O(n \log n)} expected time.
#'
#' @param series Numeric vector representing the time series. The position
#'   index is used as the time coordinate.
#' @param directed Logical. If \code{TRUE}, also returns directed adjacency
#'   lists (\code{out_adjacency} and \code{in_adjacency}) for constructing
#'   forward and backward Nada topologies or Alexandrov topologies. Edges
#'   are directed from earlier to later time (default: \code{FALSE}).
#' @return A \code{list} with the same structure as
#'   \code{\link{horizontal_visibility_graph}}.
#'
#' @details
#' The algorithm uses divide-and-conquer on the global maximum. The maximum
#' point blocks all cross-edges between the left and right halves, so the
#' problem decomposes cleanly. Edges from the pivot are found by a linear
#' slope scan. Expected time is \eqn{O(n \log n)}; worst case (monotone data)
#' is \eqn{O(n^2)}.
#'
#' When \code{directed = TRUE}, the function additionally returns directed
#' adjacency lists. By construction, all edges satisfy \code{from < to}
#' (earlier time to later time), so the directed graph is a DAG.
#'
#' @references
#' Lacasa, L., Luque, B., Ballesteros, F., Luque, J., & Nuno, J. C. (2008).
#' From time series to complex networks: The visibility graph.
#' \emph{Proceedings of the National Academy of Sciences}, 105(13), 4972-4975.
#' \doi{10.1073/pnas.0709247105}
#'
#' @examples
#' series <- c(3, 1, 4, 1, 5, 9, 2, 6)
#' g <- natural_visibility_graph(series)
#' g$n_edges
#' head(g$edges)
#'
#' # Directed graph for bitopological analysis
#' gd <- natural_visibility_graph(series, directed = TRUE)
#' gd$out_adjacency[[1]]  # forward neighbors of first observation
#'
#' @export
natural_visibility_graph <- function(series, directed = FALSE) {
  if (!is.numeric(series)) {
    stop("'series' must be a numeric vector.", call. = FALSE)
  }
  if (anyNA(series)) {
    stop("'series' must not contain NA values.", call. = FALSE)
  }
  if (!all(is.finite(series))) {
    stop("'series' must not contain Inf or NaN values.", call. = FALSE)
  }
  n <- length(series)
  if (n < 2L) {
    result <- list(
      edges = data.frame(from = integer(0), to = integer(0)),
      n = n,
      n_edges = 0L,
      adjacency = lapply(seq_len(n), function(i) integer(0))
    )
    if (directed) {
      result$out_adjacency <- lapply(seq_len(n), function(i) integer(0))
      result$in_adjacency <- lapply(seq_len(n), function(i) integer(0))
    }
    return(result)
  }

  raw <- nvg_cpp(as.double(series))

  edges <- data.frame(from = raw$from, to = raw$to)

  adj <- vector("list", n)
  for (i in seq_len(n)) adj[[i]] <- integer(0)
  for (k in seq_along(raw$from)) {
    i <- raw$from[k]
    j <- raw$to[k]
    adj[[i]] <- c(adj[[i]], j)
    adj[[j]] <- c(adj[[j]], i)
  }

  result <- list(
    edges = edges,
    n = n,
    n_edges = nrow(edges),
    adjacency = adj
  )

  if (directed) {
    out_adj <- vector("list", n)
    in_adj <- vector("list", n)
    for (i in seq_len(n)) {
      out_adj[[i]] <- integer(0)
      in_adj[[i]] <- integer(0)
    }
    for (k in seq_along(raw$from)) {
      i <- raw$from[k]
      j <- raw$to[k]
      out_adj[[i]] <- c(out_adj[[i]], j)
      in_adj[[j]] <- c(in_adj[[j]], i)
    }
    result$out_adjacency <- out_adj
    result$in_adjacency <- in_adj
  }

  result
}
