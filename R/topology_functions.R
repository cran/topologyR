#' @title Topology Analysis Functions for Time Series Data
#' @description Functions for analyzing topological properties of time series data
#' @author José Mauricio Gómez Julián
#' @details
#' This module provides three main functions for analyzing the connectivity
#' of topological structures, particularly focused on economic time series:
#' - is_topology_connected: Uses an undirected graph approach
#' - is_topology_connected2: Uses a directed graph approach
#' - is_topology_connected_manual: Uses a manual checking approach
#'
#' Created as part of research on economic cycle analysis.

#' Check if a topology is connected using undirected graph approach
#'
#' @param topology A list of sets representing the topology
#' @return \code{logical} scalar. Returns \code{TRUE} if the topology is 
#'   connected (all elements can be reached from any starting point through 
#'   the undirected graph representation), \code{FALSE} otherwise. Returns 
#'   \code{FALSE} if the topology is empty.
#' @examples
#' topology <- list(c(1,2,3), c(3,4,5))
#' is_topology_connected(topology)
#' @export
is_topology_connected <- function(topology) {
  if (!length(topology)) return(FALSE)

  elements <- unique(unlist(topology))
  n <- max(elements)  # Changed from length(elements)

  # Create undirected graph from topology
  edges <- matrix(0, nrow = n, ncol = n)
  for (set in topology) {
    for (i in set) {
      for (j in set) {
        if (i != j && i <= n && j <= n) {  # Added boundary check
          edges[i, j] <- 1
          edges[j, i] <- 1
        }
      }
    }
  }

  # Perform DFS to check connectivity
  visited <- rep(FALSE, n)
  stack <- c(elements[1])

  while (length(stack) > 0) {
    v <- stack[1]
    stack <- stack[-1]
    if (!visited[v]) {
      visited[v] <- TRUE
      neighbors <- which(edges[v, ] == 1)
      stack <- c(stack, elements[neighbors])
    }
  }

  all(visited[elements])  # Only check elements that are present
}

#' Check if a topology is connected using directed graph approach
#'
#' @param topology A list of sets representing the topology
#' @return \code{logical} scalar. Returns \code{TRUE} if the topology is 
#'   connected through directed graph traversal (considering sequential 
#'   connections between elements), \code{FALSE} otherwise. Returns 
#'   \code{FALSE} if the topology is empty.
#' @examples
#' topology <- list(c(1,2,3), c(3,4,5))
#' is_topology_connected2(topology)
#' @export
is_topology_connected2 <- function(topology) {
  # Handle empty topology case
  if (!length(topology)) return(FALSE)

  # Get unique elements and determine matrix size
  elements <- unique(unlist(topology))
  n <- max(elements)

  # Initialize adjacency matrix
  edges <- matrix(0, nrow = n, ncol = n)

  # Build directed edges
  for (set in topology) {
    if (length(set) > 1) {
      # Get values in current set
      values <- sort(set)
      # Create edges between consecutive elements
      for (i in 1:(length(values)-1)) {
        current_val <- values[i]
        next_val <- values[i+1]
        edges[current_val, next_val] <- 1
      }
    }
  }

  # Initialize DFS variables
  visited <- rep(FALSE, n)
  start <- min(elements)
  stack <- c(start)

  # Perform DFS
  while (length(stack) > 0) {
    current <- stack[1]
    stack <- stack[-1]

    if (!visited[current]) {
      visited[current] <- TRUE
      # Find all unvisited neighbors
      neighbors <- which(edges[current,] == 1)
      stack <- c(stack, neighbors[!visited[neighbors]])
    }
  }

  # Check if all elements in topology are visited
  return(all(visited[elements]))
}

#' Check if a topology is connected using manual checking approach
#'
#' @param topology A list of sets representing the topology
#' @return \code{logical} scalar. Returns \code{TRUE} if all elements from 
#'   1 to the maximum value in the topology are present in at least one set, 
#'   \code{FALSE} otherwise.
#' @examples
#' topology <- list(c(1,2,3), c(3,4,5))
#' is_topology_connected_manual(topology)
#' @export
is_topology_connected_manual <- function(topology) {
  # Get all unique elements in the topology
  all_elements <- unique(unlist(topology))

  # Get elements from original set (assuming they go from 1 to n)
  n <- max(all_elements)
  original_elements <- 1:n

  # Check if each element appears in at least one set
  for (element in original_elements) {
    if (!any(sapply(topology, function(set) element %in% set))) {
      return(FALSE)
    }
  }

  return(TRUE)
}

#' Calculate topology characteristics for different IQR factors
#'
#' @description
#' This function analyzes how different IQR (Interquartile Range) factors affect
#' the topology's characteristics. It helps users determine the optimal factor
#' for their specific data by showing how the factor choice impacts the base size
#' and set sizes in the resulting topology.
#'
#' @param data Numeric vector containing the data to analyze
#' @param factors Numeric vector of factors to test (default: c(1, 2, 4, 8, 16))
#' @param plot Logical, whether to display a plot (default: TRUE)
#' @return A \code{data.frame} with the following columns:
#'   \describe{
#'     \item{factor}{Numeric. The IQR factor used for threshold calculation.}
#'     \item{threshold}{Numeric. The calculated threshold value (IQR/factor).}
#'     \item{base_size}{Integer. Number of sets in the topological base.}
#'     \item{max_set_size}{Integer. Size of the largest set in the base.}
#'     \item{min_set_size}{Integer. Size of the smallest set in the base.}
#'   }
#'   If \code{plot = TRUE}, also generates a line plot showing the relationship
#'   between IQR factors and base sizes as a side effect.
#'
#' @details
#' The function works by:
#' 1. Calculating different thresholds using IQR/factor
#' 2. Creating a subbase using these thresholds
#' 3. Generating the base from intersections of subbase elements
#' 4. Analyzing the resulting topology's characteristics
#'
#' A larger factor results in a smaller threshold, which typically leads to
#' a finer topology with more distinct sets but smaller set sizes.
#'
#' @examples
#' # Generate sample data
#' data <- rnorm(100)
#'
#' # Analyze topology with default factors
#' results <- analyze_topology_factors(data)
#' print(results)
#'
#' # Use custom factors
#' custom_results <- analyze_topology_factors(data, factors = c(2, 4, 6))
#' print(custom_results)
#'
#' @export
analyze_topology_factors <- function(data, factors = NULL, plot = TRUE) {
  # Set default factors if not provided
  if (is.null(factors)) {
    factors <- c(1, 2, 4, 8, 16)
  }

  # Inner function to calculate topology for a single factor
  calculate_topology_with_factor <- function(data, factor) {
    threshold <- stats::IQR(data) / factor
    n <- length(data)

    # Calculate subbase using threshold
    subbase <- lapply(1:n, function(i) {
      which(abs(data - data[i]) <= threshold)
    })

    # Generate base (including empty set and full set)
    base <- list(integer(0), 1:n)
    for (i in 1:n) {
      for (j in i:n) {
        intersection <- intersect(subbase[[i]], subbase[[j]])
        if (length(intersection) > 0) {
          base <- c(base, list(intersection))
        }
      }
    }
    unique(base)
  }

  # Calculate results for each factor
  results <- lapply(factors, function(f) {
    topology <- calculate_topology_with_factor(data, f)
    list(
      factor = f,
      threshold = stats::IQR(data) / f,
      base_size = length(topology),
      max_set_size = max(sapply(topology, length)),
      min_set_size = min(sapply(topology, length))
    )
  })

  # Convert results to data frame
  results_df <- do.call(rbind, lapply(results, data.frame))

  # Optional plotting
  if (plot) {
    p <- ggplot2::ggplot(results_df, ggplot2::aes(x = factor)) +
      ggplot2::geom_line(ggplot2::aes(y = base_size, color = "Base Size")) +
      ggplot2::geom_line(ggplot2::aes(y = max_set_size, color = "Maximum Set Size")) +
      ggplot2::geom_line(ggplot2::aes(y = min_set_size, color = "Minimum Set Size")) +
      ggplot2::scale_x_log10() +
      ggplot2::labs(title = "Effect of IQR Factor on Topology",
           x = "IQR Factor", y = "Size") +
      ggplot2::theme_minimal()

    # Print the plot
    print(p)
  }

  # Return the results data frame
  return(results_df)
}

#' Calculate multiple threshold methods for topology analysis
#'
#' @param data Numeric vector to calculate thresholds for
#' @return A named \code{list} containing five threshold values:
#'   \describe{
#'     \item{mean_diff}{Numeric. Threshold based on mean of absolute differences 
#'       between adjacent sorted values.}
#'     \item{median_diff}{Numeric. Threshold based on median of absolute differences 
#'       between adjacent sorted values.}
#'     \item{sd}{Numeric. Threshold based on standard deviation of the data.}
#'     \item{iqr}{Numeric. Threshold based on IQR divided by 4.}
#'     \item{dbscan}{Numeric. Threshold based on DBSCAN-like density estimation 
#'       using k-th nearest neighbor distance.}
#'   }
#' @export
calculate_thresholds <- function(data) {
  list(
    mean_diff = mean(abs(diff(sort(data)))),
    median_diff = stats::median(abs(diff(sort(data)))),
    sd = stats::sd(data),
    iqr = stats::IQR(data) / 4,
    dbscan = {
      k <- ceiling(log(length(data)))
      sort(stats::dist(matrix(data, ncol=1)))[k * length(data)]
    }
  )
}

#' Calculate topology base size for a given threshold
#'
#' @param data Numeric vector to analyze
#' @param threshold Numeric value for the threshold parameter
#' @return \code{integer} scalar. The number of sets in the topological base. 
#'   This value represents the complexity of the topology - larger values 
#'   indicate more complex topological structure.
#' @export
calculate_topology <- function(data, threshold) {
  n <- length(data)
  subbase <- lapply(1:n, function(i) {
    which(abs(data - data[i]) <= threshold)
  })

  base <- list(integer(0), 1:n)
  for (i in 1:n) {
    for (j in i:n) {
      intersection <- intersect(subbase[[i]], subbase[[j]])
      if (length(intersection) > 0) {
        base <- c(base, list(intersection))
      }
    }
  }
  base <- unique(base)

  length(base)  # Returns the base size as a complexity measure
}

#' Visualize Topology Thresholds
#'
#' @title Visualize and Compare Different Threshold Methods
#' @description This function creates a comprehensive visualization of different
#'   threshold methods used in topology analysis. It generates three distinct plots
#'   to help understand the relationships between different threshold methods and
#'   their effects on topological structure.
#'
#' @param data Numeric vector to analyze
#' @param plot Logical indicating whether to display plots (default: TRUE)
#' @return A \code{data.frame} with one row per threshold method containing:
#'   \describe{
#'     \item{method}{Character. Name of the threshold calculation method.}
#'     \item{threshold}{Numeric. The calculated threshold value.}
#'     \item{base_size}{Integer. Number of sets in the resulting topological base.}
#'   }
#'   If \code{plot = TRUE}, generates three plots as side effects:
#'   (1) Bar chart comparing threshold values by method,
#'   (2) Bar chart comparing base sizes by method,
#'   (3) Scatter plot showing threshold vs base size relationship.
#' @importFrom ggplot2 ggplot aes geom_bar geom_point geom_text theme_minimal labs
#'
#' @examples
#' \donttest{
#' # Generate sample data
#' data <- rnorm(100)
#'
#' # Visualize threshold comparisons
#' results <- visualize_topology_thresholds(data)
#' }
#' @export
visualize_topology_thresholds <- function(data, plot = TRUE) {

  # Calculate thresholds
  thresholds <- calculate_thresholds(data)

  # Calculate base sizes (using the existing calculate_topology function)
  base_sizes <- sapply(thresholds, function(t) calculate_topology(data, t))

  # Create data frame
  df <- data.frame(
    method = names(thresholds),
    threshold = unlist(thresholds),
    base_size = base_sizes
  )

  if (plot) {
    # Threshold comparison plot
    p1 <- ggplot2::ggplot(df, ggplot2::aes(x = method, y = threshold)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Threshold comparison by method",
           x = "Method", y = "Threshold")

    # Base size comparison plot
    p2 <- ggplot2::ggplot(df, ggplot2::aes(x = method, y = base_size)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Base size comparison by method",
           x = "Method", y = "Base size")

    # Threshold vs base size plot
    p3 <- ggplot2::ggplot(df, ggplot2::aes(x = threshold, y = base_size, label = method)) +
      ggplot2::geom_point() +
      ggplot2::geom_text(hjust = -0.1, vjust = 0) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Relationship between threshold and base size",
           x = "Threshold", y = "Base size")

    # Print plots
    print(p1)
    print(p2)
    print(p3)
  }

  # Return data frame for further analysis
  return(df)
}

#' Create a topology with completely disconnected sets
#'
#' @description
#' This function generates a topology where each set in the topology
#' is a singleton (contains only one element), resulting in a completely
#' disconnected topological structure. Each vertex exists in isolation,
#' with no meaningful connections between sets.
#'
#' @details
#' The subbase contains individual elements. The base consists of singleton sets.
#' The topology is formed by these singleton sets. No meaningful topological 
#' relationships are established between elements.
#'
#' @param datos Numeric vector containing the data points to analyze
#' @return A \code{list} with four components:
#'   \describe{
#'     \item{R}{List. The equivalence relation defined on the data vertices.}
#'     \item{subbase}{List. Sets of individual vertices forming the subbase.}
#'     \item{base}{List. Singleton sets forming the base of the topology.}
#'     \item{topology}{List. The complete topology consisting of disconnected singleton sets.}
#'   }
#' @examples
#' data <- c(1, 2, 3, 4, 5)
#' result <- simplest_topology(data)
#' @export
simplest_topology <- function(datos) {
  # Step 1: Get number of vertices
  n <- length(datos)
  # Step 2: Create complete graph from data
  vertices <- seq_len(n)
  adj_matrix <- matrix(1, nrow = n, ncol = n)
  diag(adj_matrix) <- 0
  # Step 3: Assign data values as edge weights
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        adj_matrix[i, j] <- abs(datos[i] - datos[j])
      }
    }
  }
  # Step 4: Calculate vertex degrees
  degree <- rowSums(adj_matrix > 0, na.rm = TRUE)
  # Step 5: Define relation R and classes
  R <- list()
  class <- list()
  for (i in 1:n) {
    for (j in 1:n) {
      if (degree[i] != 0) {
        class[[paste(i, j)]] <- vertices[j]
        R[[paste(vertices[i], vertices[j])]] <- list(x = paste0(degree[i], "_", datos[i]),
                                                     y = paste0(degree[j], "_", datos[j]))
      }
    }
  }
  # Step 6: Get subbase
  subbase <- list()
  for (i in 1:n) {
    for (j in 1:n) {
      if (!is.null(class[[paste(i, j)]])) {
        subbase <- append(subbase, class[[paste(i, j)]])
      }
    }
  }
  subbase <- unique(subbase)
  # Step 7: Get base
  base <- subbase
  for (vi in vertices) {
    for (vj in vertices) {
      base <- append(base, intersect(vi, vj))
    }
  }
  base <- unique(base)
  # Step 8: Get topology
  topology <- base
  for (vi in vertices) {
    for (vj in vertices) {
      union_set <- union(vi, vj)
      if (length(union_set) > 0) {
        topology <- append(topology, union_set)
      }
    }
  }
  topology <- unique(topology)

  list(R = R, subbase = subbase, base = base, topology = topology)
}

#' Create a complete topology from data points with neighborhood structure
#'
#' @param datos Numeric vector containing the data points to analyze
#' @return A \code{list} with four components:
#'   \describe{
#'     \item{R}{List. The equivalence relation defined on the data, where each
#'       element represents a relationship between vertices based on their degrees
#'       and values.}
#'     \item{subbase}{List. The neighborhoods of each vertex, where each element
#'       contains the indices of neighboring vertices in the graph.}
#'     \item{base}{List. The base generated from intersections of neighborhoods,
#'       including the empty set, full set, and all non-empty intersections.}
#'     \item{topology}{List. The complete topology generated by taking unions of 
#'       base elements, satisfying the axioms of topological spaces.}
#'   }
#' @examples
#' data <- c(1, 2, 3, 4, 5)
#' result <- complete_topology(data)
#' @export
complete_topology <- function(datos) {
  # Step 1: Get number of vertices
  n <- length(datos)
  # Step 2: Create complete graph from data
  vertices <- seq_len(n)
  adj_matrix <- matrix(1, nrow = n, ncol = n)
  diag(adj_matrix) <- 0
  # Step 3: Assign data values as edge weights
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        adj_matrix[i, j] <- abs(datos[i] - datos[j])
      }
    }
  }
  # Step 4: Calculate vertex degrees
  degree <- rowSums(adj_matrix > 0, na.rm = TRUE)
  # Step 5: Define relation R and classes
  R <- list()
  class <- list()
  for (i in 1:n) {
    for (j in 1:n) {
      if (degree[i] != 0) {
        class[[paste(i, j)]] <- vertices[j]
        R[[paste(vertices[i], vertices[j])]] <- list(x = paste0(degree[i], "_", datos[i]),
                                                     y = paste0(degree[j], "_", datos[j]))
      }
    }
  }
  # Step 6: Get subbase (neighborhoods of each vertex)
  subbase <- list()
  for (i in 1:n) {
    neighbors <- which(adj_matrix[i, ] > 0)
    subbase[[i]] <- vertices[neighbors]
  }
  # Step 7: Get base
  base <- list(integer(0), vertices)
  for (i in 1:n) {
    for (j in i:n) {
      intersection <- intersect(subbase[[i]], subbase[[j]])
      if (length(intersection) > 0) {
        base <- append(base, list(intersection))
      }
    }
  }
  base <- unique(base)
  # Step 8: Get topology
  topology <- base
  for (i in 1:length(base)) {
    for (j in i:length(base)) {
      union_set <- union(base[[i]], base[[j]])
      if (length(union_set) > 0) {
        topology <- append(topology, list(union_set))
      }
    }
  }
  topology <- unique(topology)

  list(R = R, subbase = subbase, base = base, topology = topology)
}