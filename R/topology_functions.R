#' Check if a topology is connected using undirected graph approach
#'
#' @description
#' Converts the topology into an undirected graph (two elements share an edge
#' if they appear together in any open set) and checks graph connectivity via
#' depth-first search.
#'
#' \strong{Note:} This is a necessary but not sufficient condition for
#' topological connectivity. For an exact check, use
#' \code{\link{is_topology_connected_exact}}.
#'
#' @param topology A list of integer vectors representing the open sets.
#' @return \code{logical} scalar. \code{TRUE} if the derived graph is connected,
#'   \code{FALSE} otherwise or if the topology is empty.
#' @examples
#' topology <- list(c(1, 2, 3), c(3, 4, 5))
#' is_topology_connected(topology)
#' @export
is_topology_connected <- function(topology) {
  if (!length(topology)) return(FALSE)

  elements <- unique(unlist(topology))
  if (length(elements) <= 1L) return(TRUE)

  n <- max(elements)
  edges <- matrix(0L, nrow = n, ncol = n)
  for (set in topology) {
    for (i in set) {
      for (j in set) {
        if (i != j && i <= n && j <= n) {
          edges[i, j] <- 1L
          edges[j, i] <- 1L
        }
      }
    }
  }

  visited <- rep(FALSE, n)
  stack <- elements[1L]
  while (length(stack) > 0L) {
    v <- stack[1L]
    stack <- stack[-1L]
    if (!visited[v]) {
      visited[v] <- TRUE
      neighbors <- which(edges[v, ] == 1L)
      stack <- c(stack, neighbors[!visited[neighbors]])
    }
  }

  all(visited[elements])
}


#' Check if a topology is connected using directed graph approach
#'
#' @description
#' Converts the topology into a directed graph (sequential edges between
#' sorted elements within each set) and checks reachability via DFS.
#'
#' \strong{Note:} This is a necessary but not sufficient condition for
#' topological connectivity. For an exact check, use
#' \code{\link{is_topology_connected_exact}}.
#'
#' @param topology A list of integer vectors representing the open sets.
#' @return \code{logical} scalar. \code{TRUE} if all elements are reachable
#'   from the minimum element via directed edges, \code{FALSE} otherwise.
#' @examples
#' topology <- list(c(1, 2, 3), c(3, 4, 5))
#' is_topology_connected2(topology)
#' @export
is_topology_connected2 <- function(topology) {
  if (!length(topology)) return(FALSE)

  elements <- unique(unlist(topology))
  if (length(elements) <= 1L) return(TRUE)

  n <- max(elements)
  edges <- matrix(0L, nrow = n, ncol = n)

  for (set in topology) {
    if (length(set) > 1L) {
      values <- sort(set)
      for (i in seq_len(length(values) - 1L)) {
        edges[values[i], values[i + 1L]] <- 1L
      }
    }
  }

  visited <- rep(FALSE, n)
  start <- min(elements)
  stack <- start
  while (length(stack) > 0L) {
    current <- stack[1L]
    stack <- stack[-1L]
    if (!visited[current]) {
      visited[current] <- TRUE
      neighbors <- which(edges[current, ] == 1L)
      stack <- c(stack, neighbors[!visited[neighbors]])
    }
  }

  all(visited[elements])
}


#' Check if all elements are covered by the topology
#'
#' @description
#' Checks whether every integer from 1 to the maximum value in the topology
#' appears in at least one open set. This verifies \strong{coverage}
#' (\eqn{\bigcup \tau = X}), not topological connectivity.
#'
#' \strong{Note:} A space can have full coverage and be maximally disconnected
#' (e.g., the discrete topology). For connectivity, use
#' \code{\link{is_topology_connected_exact}}.
#'
#' @param topology A list of integer vectors representing the open sets.
#' @return \code{logical} scalar. \code{TRUE} if all elements from 1 to the
#'   maximum are present in at least one set.
#' @examples
#' topology <- list(c(1, 2, 3), c(3, 4, 5))
#' is_topology_connected_manual(topology)
#' @export
is_topology_connected_manual <- function(topology) {
  all_elements <- unique(unlist(topology))
  if (length(all_elements) == 0L) return(FALSE)
  n <- max(all_elements)
  all(seq_len(n) %in% all_elements)
}


#' Analyze topology characteristics for different IQR factors
#'
#' @description
#' Analyzes how different IQR factors affect topology characteristics.
#' Helps determine the optimal factor by showing how the factor choice
#' impacts base size and set sizes.
#'
#' \strong{Note:} This function uses threshold-based neighborhoods (metric
#' approximation), not the graph-theoretic approach of Nada et al. (2018).
#' For the theoretically faithful approach, use visibility graphs with
#' \code{\link{generate_topology}}.
#'
#' @param data Numeric vector containing the data to analyze.
#' @param factors Numeric vector of factors to test
#'   (default: \code{c(1, 2, 4, 8, 16)}).
#' @param plot Logical, whether to return a plot object (default: \code{TRUE}).
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{factor}{Numeric. The IQR factor used.}
#'     \item{threshold}{Numeric. The calculated threshold (IQR/factor).}
#'     \item{base_size}{Integer. Number of sets in the base.}
#'     \item{max_set_size}{Integer. Size of the largest set.}
#'     \item{min_set_size}{Integer. Size of the smallest set.}
#'   }
#'   If \code{plot = TRUE}, the data.frame carries an attribute \code{"plot"}
#'   containing a \code{ggplot} object.
#'
#' @examples
#' data <- rnorm(50)
#' results <- analyze_topology_factors(data)
#' print(results)
#'
#' @export
analyze_topology_factors <- function(data, factors = NULL, plot = TRUE) {
  if (!is.numeric(data) || length(data) < 2L) {
    stop("'data' must be a numeric vector with at least 2 elements.",
         call. = FALSE)
  }
  if (anyNA(data)) {
    stop("'data' must not contain NA values.", call. = FALSE)
  }
  if (is.null(factors)) {
    factors <- c(1, 2, 4, 8, 16)
  }

  calculate_topology_with_factor <- function(data, factor) {
    threshold <- stats::IQR(data) / factor
    n <- length(data)
    subbase <- lapply(seq_len(n), function(i) {
      which(abs(data - data[i]) <= threshold)
    })
    base <- list(integer(0), seq_len(n))
    for (i in seq_len(n)) {
      for (j in seq(i, n)) {
        inter <- intersect(subbase[[i]], subbase[[j]])
        if (length(inter) > 0L) {
          base <- c(base, list(inter))
        }
      }
    }
    unique(base)
  }

  results <- lapply(factors, function(f) {
    topo <- calculate_topology_with_factor(data, f)
    sizes <- vapply(topo, length, integer(1))
    list(
      factor = f,
      threshold = stats::IQR(data) / f,
      base_size = length(topo),
      max_set_size = max(sizes),
      min_set_size = min(sizes)
    )
  })

  results_df <- do.call(rbind, lapply(results, data.frame))

  if (plot) {
    p <- ggplot2::ggplot(results_df, ggplot2::aes(x = factor)) +
      ggplot2::geom_line(ggplot2::aes(y = base_size, color = "Base Size")) +
      ggplot2::geom_line(ggplot2::aes(y = max_set_size,
                                      color = "Maximum Set Size")) +
      ggplot2::geom_line(ggplot2::aes(y = min_set_size,
                                      color = "Minimum Set Size")) +
      ggplot2::scale_x_log10() +
      ggplot2::labs(title = "Effect of IQR Factor on Topology",
                    x = "IQR Factor", y = "Size") +
      ggplot2::theme_minimal()
    attr(results_df, "plot") <- p
  }

  results_df
}


#' Calculate multiple threshold methods for topology analysis
#'
#' @description
#' Computes five different threshold methods for defining neighborhoods
#' in threshold-based topology construction.
#'
#' \strong{Note:} Threshold-based methods produce metric topologies, not
#' the graph-induced topologies of Nada et al. (2018). See
#' \code{\link{generate_topology}} for the theoretically faithful approach.
#'
#' @param data Numeric vector to calculate thresholds for.
#' @return A named \code{list} with five threshold values:
#'   \describe{
#'     \item{mean_diff}{Mean of absolute differences between adjacent sorted values.}
#'     \item{median_diff}{Median of absolute differences between adjacent sorted values.}
#'     \item{sd}{Standard deviation of the data.}
#'     \item{iqr}{IQR divided by 4.}
#'     \item{dbscan}{k-th nearest neighbor distance (k = ceiling(log(n))).}
#'   }
#' @examples
#' calculate_thresholds(rnorm(100))
#' @export
calculate_thresholds <- function(data) {
  if (!is.numeric(data) || length(data) < 2L) {
    stop("'data' must be a numeric vector with at least 2 elements.",
         call. = FALSE)
  }
  if (anyNA(data)) {
    stop("'data' must not contain NA values.", call. = FALSE)
  }
  sorted_data <- sort(data)
  diffs <- abs(diff(sorted_data))
  k <- ceiling(log(length(data)))
  dist_sorted <- sort(as.numeric(stats::dist(matrix(data, ncol = 1))))

  list(
    mean_diff = mean(diffs),
    median_diff = stats::median(diffs),
    sd = stats::sd(data),
    iqr = stats::IQR(data) / 4,
    dbscan = dist_sorted[min(k * length(data), length(dist_sorted))]
  )
}


#' Calculate topology base size for a given threshold
#'
#' @param data Numeric vector to analyze.
#' @param threshold Numeric value for the threshold parameter. Must be
#'   non-negative.
#' @return \code{integer} scalar. The number of sets in the topological base.
#' @examples
#' calculate_topology(rnorm(30), threshold = 0.5)
#' @export
calculate_topology <- function(data, threshold) {
  if (!is.numeric(data) || length(data) < 2L) {
    stop("'data' must be a numeric vector with at least 2 elements.",
         call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || threshold < 0) {
    stop("'threshold' must be a single non-negative number.", call. = FALSE)
  }
  n <- length(data)
  subbase <- lapply(seq_len(n), function(i) {
    which(abs(data - data[i]) <= threshold)
  })

  base <- list(integer(0), seq_len(n))
  for (i in seq_len(n)) {
    for (j in seq(i, n)) {
      inter <- intersect(subbase[[i]], subbase[[j]])
      if (length(inter) > 0L) {
        base <- c(base, list(inter))
      }
    }
  }
  base <- unique(base)
  length(base)
}


#' Visualize and compare different threshold methods
#'
#' @param data Numeric vector to analyze.
#' @param plot Logical indicating whether to return plot objects
#'   (default: \code{TRUE}).
#' @return A \code{data.frame} with columns \code{method}, \code{threshold},
#'   and \code{base_size}. If \code{plot = TRUE}, carries an attribute
#'   \code{"plots"} containing a list of three \code{ggplot} objects.
#' @importFrom ggplot2 ggplot aes geom_bar geom_point geom_text theme_minimal labs
#' @examples
#' \donttest{
#' data <- rnorm(50)
#' results <- visualize_topology_thresholds(data)
#' }
#' @export
visualize_topology_thresholds <- function(data, plot = TRUE) {
  if (!is.numeric(data) || length(data) < 2L) {
    stop("'data' must be a numeric vector with at least 2 elements.",
         call. = FALSE)
  }
  if (anyNA(data)) {
    stop("'data' must not contain NA values.", call. = FALSE)
  }

  thresholds <- calculate_thresholds(data)
  base_sizes <- vapply(thresholds, function(t) calculate_topology(data, t),
                        integer(1))

  df <- data.frame(
    method = names(thresholds),
    threshold = unlist(thresholds),
    base_size = base_sizes,
    stringsAsFactors = FALSE
  )

  if (plot) {
    p1 <- ggplot2::ggplot(df, ggplot2::aes(x = method, y = threshold)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Threshold comparison by method",
                    x = "Method", y = "Threshold")

    p2 <- ggplot2::ggplot(df, ggplot2::aes(x = method, y = base_size)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Base size comparison by method",
                    x = "Method", y = "Base size")

    p3 <- ggplot2::ggplot(df, ggplot2::aes(x = threshold, y = base_size,
                                            label = method)) +
      ggplot2::geom_point() +
      ggplot2::geom_text(hjust = -0.1, vjust = 0) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Relationship between threshold and base size",
                    x = "Threshold", y = "Base size")

    attr(df, "plots") <- list(threshold = p1, base_size = p2, scatter = p3)
  }

  df
}


#' Create the discrete topology (completely disconnected)
#'
#' @description
#' Generates the discrete topology on \code{n} elements, where the base
#' consists of all singleton sets \eqn{\{\{1\}, \{2\}, \ldots, \{n\}\}}.
#' The full topology is the power set \eqn{2^V}, but only the base is
#' returned explicitly (enumerating \eqn{2^n} sets is impractical for
#' large \code{n}).
#'
#' @param data Numeric vector containing the data points.
#' @return A \code{list} with:
#'   \describe{
#'     \item{subbase}{List of singleton sets.}
#'     \item{base}{Same as subbase (singletons are already a base).}
#'     \item{topology_type}{Character: \code{"discrete"}.}
#'     \item{n}{Integer. Number of elements.}
#'     \item{connected}{Logical. Always \code{FALSE} for \eqn{n \geq 2}
#'       (the discrete topology is maximally disconnected).}
#'   }
#' @examples
#' result <- simplest_topology(c(1, 2, 3, 4, 5))
#' result$connected
#' @export
simplest_topology <- function(data) {
  if (!is.numeric(data) || length(data) < 1L) {
    stop("'data' must be a numeric vector with at least 1 element.",
         call. = FALSE)
  }
  n <- length(data)
  singletons <- lapply(seq_len(n), function(i) i)

  list(
    subbase = singletons,
    base = singletons,
    topology_type = "discrete",
    n = n,
    connected = n <= 1L
  )
}


#' Create a complete-graph topology from data
#'
#' @description
#' Constructs a topology from the complete graph on the data points.
#' In the complete graph, every vertex is a neighbor of every other vertex,
#' so each neighborhood is \eqn{N(v) = V \setminus \{v\}}. The resulting
#' topology is generated using the fixed-point closure algorithm.
#'
#' \strong{Note:} The complete graph produces trivial neighborhoods
#' (all vertices except self). For meaningful topological structure,
#' use visibility graphs with \code{\link{generate_topology}}.
#'
#' @param data Numeric vector containing the data points. Length must be
#'   between 2 and 64.
#' @param verify_axioms Logical. Whether to verify topology axioms
#'   (default: \code{FALSE}).
#' @return A \code{list} with the same structure as \code{\link{generate_topology}}.
#' @examples
#' result <- complete_topology(c(1, 2, 3, 4, 5))
#' result$n_open_sets
#' result$connected
#' @export
complete_topology <- function(data, verify_axioms = FALSE) {
  if (!is.numeric(data) || length(data) < 2L) {
    stop("'data' must be a numeric vector with at least 2 elements.",
         call. = FALSE)
  }
  if (anyNA(data)) {
    stop("'data' must not contain NA values.", call. = FALSE)
  }
  n <- length(data)

  vertices <- seq_len(n)
  # Complete graph: each vertex neighbors all others
  adjacency <- lapply(vertices, function(v) setdiff(vertices, v))

  generate_topology(
    adjacency = adjacency,
    n_elements = n,
    check_connected = TRUE,
    verify_axioms = verify_axioms
  )
}
