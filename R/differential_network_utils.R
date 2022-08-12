# differential_network_utils.R
#
# This script contains helper functions for simulating differential networks 
# for the epiregulon project. 
# 
# Author: Timothy Keyes 
# Version 2022.11.08



#' Simulate a random weight matrix.
#' 
#' Simulate a random U or V weight matrix with sparsity encoded by a random indicator
#' matrix with entries set to 0 or 1 based on a bernoulli distribution. Entries
#' in the matrix itself are normally distributed.
#'
#' @param num_rows The number of rows in the simulated weight matrix 
#' @param num_cols The number of columns in the simulated weight matrix
#' @param bernoulli_param A numeric vector specifying the the probability of 
#' any entry in the indicator matrix being equal to 1. If bernouilli_param is 
#' length 1, the same probability will be applied across all rows of the indicator
#' matrix. If bernoulli_param is of length num_rows, each row will be given a 
#' unique amount of sparsity. 
#'
#' @return A [num_rows x num_cols] matrix. 
#' 
#'
simulate_weight <- function(num_rows, num_cols, bernoulli_param = 0.25) { 
  # create a random sequence of normally distributed edges
  random_sequence_normal <- rnorm(n = num_rows * num_cols)
  
  # if bernoulli param is length 1, 
  # randomly pick bernoulli_param * num_rows * num_cols edges, and set all 
  # others to 0 for sparsity
  if (length(bernoulli_param) == 1) { 
    random_sequence_bernoulli <- 
      rbinom(n = num_rows * num_cols, size = 1, prob = bernoulli_param)
    
    # create indicator matrix
    indicator_matrix <- 
      matrix(data = random_sequence_bernoulli, nrow = num_rows, ncol = num_cols)
    
  } else if (length(bernoulli_param == num_rows)) { 
    # create indicator matrix with each row having a unique sparsity
    indicator_matrix <- 
      bernoulli_param |> 
      matrix(nrow = num_rows, ncol = 1) |> 
      apply(MARGIN = 1, FUN = function(x) rbinom(n = num_cols, size = 1, prob = x)) |> 
      t()
    
  } else { 
      stop("bernoulli_param must be length 1 or length num_rows.")
  }
  
  # create base matrix 
  matrix_base <- 
    matrix(data = random_sequence_normal, nrow = num_rows, ncol = num_cols)
  
  # return result 
  matrix_final <- matrix_base * indicator_matrix
  
  return(matrix_final)
}

#' Simulate a random U matrix. 
#' 
#' Simulate a random U weight matrix with sparsity encoded by a random indicator
#' matrix with entries set to 0 or 1 based on a bernoulli distribution. Entries
#' in the matrix itself are normally distributed.
#' @param num_tfs The number of TFs in the simulated weight matrix 
#' @param num_res The number of REs in the simulated weight matrix 
#' @param bernoulli_param The probability of any entry in the indicator 
#' matrix being equal to 1. 
#'
#' @return A [num_tfs x num_res] matrix.
#' 
#'
simulate_u <- function(num_tfs, num_res, bernoulli_param = 0.25) { 
  
  if (length(bernoulli_param) != 1 & length(bernoulli_param) != num_tfs) { 
    stop("bernoulli_param must be either length 1 or length num_tfs")
  }
  
  u_matrix <- 
    simulate_weight(
      num_rows = num_tfs, 
      num_cols = num_res, 
      bernoulli_param = bernoulli_param
    )
  
  rownames(u_matrix) <- paste0("TF", 1:num_tfs)
  colnames(u_matrix) <- paste0("RE", 1:num_res)
  
  return(u_matrix)
}

#' Simulate a random V matrix. 
#' 
#' Simulate a random V weight matrix with sparsity encoded by a random indicator
#' matrix with entries set to 0 or 1 based on a bernoulli distribution. Entries
#' in the matrix itself are normally distributed.

#' @param num_res The number of REs in the simulated weight matrix 
#' @param num_tgs The number of TGs in the simulated weight matrix 
#' @param bernoulli_param The probability of any entry in the indicator 
#' matrix being equal to 1. 
#'
#' @return A [num_res x num_tgs] matrix.
#'
simulate_v <- function(num_res, num_tgs, bernoulli_param = 0.25) { 
  
  if (length(bernoulli_param) != 1 & length(bernoulli_param) != num_res) { 
    stop("bernoulli_param must be either length 1 or length num_res")
  }
  
  v_matrix <- 
    simulate_weight(
      num_rows = num_res, 
      num_cols = num_tgs, 
      bernoulli_param = bernoulli_param
    )
  
  rownames(v_matrix) <- paste0("RE", 1:num_res)
  colnames(v_matrix) <- paste0("TG", 1:num_tgs)
  
  return(v_matrix)
}



#' Permute edges of a network/graph adjacency matrix 
#'
#' @param weight_matrix adjacency matrix of the network in which to permute 
#' edges randomly
#' @param permutation_prop the proportion of total edges in the network that you
#' want to permute
#'
#' @return an adjacency matrix of the same size as `weight_matrix` with a
#' `permutation_prop` proportion of non-zero edges swapped with weight 0 edges
#' in the original `weight_matrix`. 
#'
permute_edges <- 
  function(weight_matrix, permutation_prop = 0.1, permutation_num = NULL) { 
    
    # find positive and null edges
    positive_indices <- which(abs(weight_matrix) > 0)
    null_indices <- 
      setdiff(1:(nrow(weight_matrix) * ncol(weight_matrix)), positive_indices)
    
    # select positive and null edges to sample
    if (!is.null(permutation_num)) { 
      num_edges_to_sample <- permutation_num
    } else {
      num_edges_to_sample <- floor((length(positive_indices) * permutation_prop)/2)
    }
    positive_edges_to_permute <- 
      sample(positive_indices, size = num_edges_to_sample)
    null_edges_to_permute <- 
      sample(null_indices, size = num_edges_to_sample)
    
    # permute edges 
    positive_weights <- weight_matrix[positive_edges_to_permute]
    weight_matrix[positive_edges_to_permute] <- 0
    weight_matrix[null_edges_to_permute] <- positive_weights
    
    return(weight_matrix)
  }

add_edges <- function(weight_matrix, addition_prop = 0.1, addition_num = NULL) { 
  # find positive and null edges
  positive_indices <- which(abs(weight_matrix) > 0)
  null_indices <- 
    setdiff(1:(nrow(weight_matrix) * ncol(weight_matrix)), positive_indices)
  
  # find number of edges to add
  if (!is.null(addition_num)) { 
    num_edges_to_add <- addition_num
  } else {
    num_edges_to_add <- addition_prop * length(positive_indices)
  }
  num_edges_to_add <- min(num_edges_to_add, length(null_indices))
  
  # sample null edges to pick which ones to turn into positive values 
  null_edges_to_add <- 
    sample(null_indices, size = num_edges_to_add)
  
  # change weights 
  weight_matrix[null_edges_to_add] <- rnorm(n = num_edges_to_add)
  
  return(weight_matrix)
}


remove_edges <- function(weight_matrix, remove_prop = 0.1, remove_num = NULL) { 
  # find positive and null edges
  positive_indices <- which(abs(weight_matrix) > 0)
  
  # find number of edges to remove 
  if (!is.null(remove_num)) { 
    num_edges_to_remove <- remove_num 
  } else {
    num_edges_to_remove <- remove_prop * length(positive_indices)
  }
  num_edges_to_remove <- min(num_edges_to_remove, length(positive_indices))
  
  # sample null edges to pick which ones to turn into positive values 
  edges_to_remove <- sample(positive_indices, size = num_edges_to_remove)
  
  # change weights 
  weight_matrix[edges_to_remove] <- 0
  
  return(weight_matrix)
}

# add_dropout <- function(weight_matrix, dropout_prop = 0.05) { 
#   # find positive and null edges
#   positive_indices <- which(abs(weight_matrix) > 0)
#   
#   # find number of edges to remove 
#   num_edges_to_remove <- dropout_prop * length(positive_indices)
#   num_edges_to_remove <- min(num_edges_to_remove, length(positive_indices))
#   
#   # sample null edges to pick which ones to turn into positive values 
#   edges_to_remove <- sample(positive_indices, size = num_edges_to_remove)
#   
#   # change weights 
#   weight_matrix[edges_to_remove] <- 0
#   
#   return(weight_matrix)
# }


add_noise <- 
  function(
    weight_matrix, 
    addition_prop = 0.1, 
    remove_prop = 0.1, 
    permutation_prop = 0.1, 
    addition_num = NULL, 
    remove_num = NULL, 
    permutation_num = NULL
  ) { 
    result <- 
      weight_matrix |> 
      add_edges(addition_prop = addition_prop, addition_num  = addition_num) |> 
      remove_edges(remove_prop = remove_prop, remove_num = remove_num) |> 
      permute_edges(
        permutation_prop = permutation_prop, 
        permutation_num = permutation_num
      )
    
    return(result)
  }


#' Permute edges of a network/graph edge list in tidy format 
#'
#' @param weight_matrix an edge list in tidy format with 4 columns (tf, re, 
#' weight, and edge_id).  
#' @param permutation_prop the proportion of total edges in the network that you
#' want to permute
#'
#' @return a tidy edge list for the network with a 
#' [permutation_percentage x nrow(weight_matrix)] number of non-zero edges set to 
#' 0 and permutation_percentage x nrow(weight_matrix) number of zero-edges given
#' random weights (sampled from the non-zero edges that were removed from the graph). 
#'
permute_edges_legacy <- function(weight_matrix, permutation_prop = 0.1) { 
  positive_edges <- 
    weight_matrix |> 
    filter(abs(weight) > 0)
  
  sparse_edges <- 
    weight_matrix |> 
    filter(weight == 0)
  
  num_edges_to_permute <- 
    (nrow(weight_matrix) * permutation_prop) |> 
    floor()
  
  # should replace be TRUE or FALSE?
  positive_indices <- 
    sample(1:nrow(positive_edges), size = num_edges_to_permute, replace = TRUE)
  positive_weights <- positive_edges$weight[positive_indices]
  sparse_indices <- sample(1:nrow(sparse_edges), size = num_edges_to_permute)
  
  # replace edges 
  positive_edges$weight[positive_indices] <- 0
  sparse_edges$weight[sparse_indices] <- sample(positive_weights)
  
  # bind together final result 
  result <- 
    bind_rows(positive_edges, sparse_edges) |> 
    mutate(edge_id = factor(edge_id, levels = weight_matrix$edge_id)) |> 
    arrange(edge_id) |> 
    mutate(edge_id = as.character(edge_id))
  
  return(result)
}

#' Add a true differential signal to a network edge list in tidy format. 
#'
#' @param weight_matrix An edge list in tidy format with 4 columns (tf, re, weight, and edge_id).  
#' @param num_tfs_to_modulate The number of TFs that should be altered such that they have TRUE differential signal. 
#' @param prop_edges_to_enrich The number of additional edges to add to each 
#' modulated TF expressed as a proportion of its initial number of edges.
#' @param prop_edges_to_deplete The number of edges to remove from each 
#' modulated TF expressed as a proportion of its initial number of edges.
#' @param prop_edges_to_permute The number of edges to permute from each 
#' modulated TF expressed as a proportion of its initial number of edges.
#'
#' @return A tidy edge list for the network with random edges added/removed 
#' according to the specified values of prop_edges_to_enrich 
#' and prop_edges_to_deplete.
#' 
#'
add_signal_legacy <- 
  function(
    weight_matrix, 
    num_tfs_to_modulate = 3L, 
    prop_edges_to_enrich = 0.1, 
    prop_edges_to_deplete = 0, 
    prop_edges_to_permute = 0
  ) { 
    weight_matrix <- 
      weight_matrix |> 
      mutate(tf = factor(tf, levels = paste0("TF", 1:nrow(weight_matrix))))
    
    num_edges_to_change <-
      weight_matrix |> 
      filter(abs(weight) > 0) |> 
      count(tf) |> 
      transmute(
        tf, 
        num_to_add = n * prop_edges_to_enrich, 
        num_to_remove = n * prop_edges_to_deplete
      ) |> 
      slice_head(n = num_tfs_to_modulate)
    
    # initialize result     
    #result <- weight_matrix
    modulated_tfs <- 
      weight_matrix |> 
      filter(tf %in% dplyr::pull(num_edges_to_change, tf))
    unmodulated_tfs <- 
      weight_matrix |> 
      filter(!(tf %in% dplyr::pull(num_edges_to_change, tf)))
    
    result <- modulated_tfs
    
    for (i in 1:num_tfs_to_modulate) { 
      my_tf <- num_edges_to_change$tf[[i]]
      num_edges_to_add <- num_edges_to_change$num_to_add[[i]]
      num_edges_to_remove <- num_edges_to_change$num_to_remove[[i]]
      
      #tf_indices <- which(weight_matrix$tf == my_tf)
      # sparse_edge_indices <- 
      #   tf_indices[which(weight_matrix$weight[tf_indices] == 0)]
      # nonzero_edge_indices <- 
      #   tf_indices[which(abs(weight_matrix$weight[tf_indices]) > 0)]
      
      tf_indices <- which(modulated_tfs$tf == my_tf)
      sparse_edge_indices <- 
        tf_indices[which(modulated_tfs$weight[tf_indices] == 0)]
      nonzero_edge_indices <- 
        tf_indices[which(abs(modulated_tfs$weight[tf_indices]) > 0)]
      
      edges_to_add <- sample(sparse_edge_indices, size = num_edges_to_add)
      edges_to_remove <- 
        sample(nonzero_edge_indices, size = num_edges_to_remove)
      
      result$weight[edges_to_add] <- rnorm(n = length(edges_to_add))
      result$weight[edges_to_remove] <- 0
    }
    
    result <- dplyr::bind_rows(result, unmodulated_tfs)
    
    return(result)
    
  }

#' Enrich the number of nonzero edges for a TF  
#' 
#' This is a helper function for `add_signal()`
#'
#' @param vec a numeric vector
#' @param prop_edges_to_enrich the proportion of zero-weight edges in vec to 
#' change to a random value sampled from a standard normal distribution
#'
#' @return a new vector of the same length as `vec`, with additional nonzero 
#' edges added. 
#'
enrich_tf <- function(vec, prop_edges_to_enrich) { 
  num_positive_edges <- 
    (abs(vec) > 0) |> 
    sum()
  
  num_edges_to_add <- 
    num_positive_edges |> 
    {function(x) x * prop_edges_to_enrich}() |> 
    floor()
  
  num_edges_to_add <- min(num_edges_to_add, length(vec) - num_positive_edges)
  
  # choose entries to enrich
  indices_to_enrich <- 
    sample(which(abs(vec) == 0), size = num_edges_to_add) |> 
    as.numeric()
  
  # add edges
  result <- vec
  result[indices_to_enrich] <- rnorm(n = num_edges_to_add)
  
  return(result)
  
}


#' Deplete the number of nonzero edges for a TF  
#' 
#' This is a helper function for `add_signal()`
#'
#' @param vec a numeric vector
#' @param prop_edges_to_deplete the proportion of zero-weight edges in vec to 
#' change to a random value sampled from a standard normal distribution
#'
#' @return a new vector of the same length as `vec`, with some nonzero 
#' edges removed.
#'
deplete_tf <- function(vec, prop_edges_to_deplete) { 
  num_positive_edges <- 
    (abs(vec) > 0) |> 
    sum()
  
  num_edges_to_remove <- 
    num_positive_edges |> 
    {function(x) x * prop_edges_to_deplete}() |> 
    floor()
  
  num_edges_to_remove <- min(num_edges_to_remove, num_positive_edges)
  
  # choose entries to enrich
  indices_to_remove <- 
    sample(which(abs(vec) > 0), size = num_edges_to_remove) |> 
    as.numeric()
  
  # add edges
  result <- vec
  result[indices_to_remove] <- 0
  
  return(result)
  
}


#' Add a true differential signal to a network edge list in tidy format. 
#'
#' @param weight_matrix an adjacency matrix of the network in which to modulate
#' specific TFs 
#' @param num_tfs_to_modulate The number of TFs that should be altered such 
#' that they have REAL differential signal. Defaults to 1/3 TFs in the dataset.
#' @param prop_edges_to_enrich The number of additional edges to add to each 
#' modulated TF expressed as a proportion of its initial number of edges.
#' @param prop_edges_to_deplete The number of edges to remove from each 
#' modulated TF expressed as a proportion of its initial number of edges.
#'
#' @return An adjacency matrix for the network with random edges added/removed 
#' according to the specified values of prop_edges_to_enrich 
#' and prop_edges_to_deplete.
#' 
#'
add_signal <- 
  function(
    weight_matrix, 
    num_tfs_to_modulate = 0.3 * nrow(weight_matrix), 
    prop_edges_to_enrich = 0.1, 
    prop_edges_to_deplete = 0, 
    prop_edges_to_permute = 0
  ) { 
    result <- weight_matrix
    # add edges
    if (prop_edges_to_enrich > 0) { 
      modulated_tfs <- 
        weight_matrix[1:num_tfs_to_modulate, ] |> 
        apply(
          MARGIN = 1, 
          FUN = 
            function(x) enrich_tf(vec = x, prop_edges_to_enrich = prop_edges_to_enrich)
        ) |> 
        t()
      
      result[1:num_tfs_to_modulate, ] <- modulated_tfs
    }
    
    # remove edges
    if (prop_edges_to_deplete > 0) { 
      modulated_tfs <- 
        result[1:num_tfs_to_modulate, ] |> 
        apply(
          MARGIN = 1, 
          FUN = 
            function(x) deplete_tf(vec = x, prop_edges_to_deplete = prop_edges_to_deplete)
        ) |> 
        t()
      
      result[1:num_tfs_to_modulate, ] <- modulated_tfs
    }
    
    # permute edges
    if (prop_edges_to_permute > 0) {
      modulated_tfs <- 
        result[1:num_tfs_to_modulate, ] |> 
        apply(
          MARGIN = 1, 
          FUN = 
            function(x) permute_edges(weight_matrix = as.matrix(x), permutation_prop = prop_edges_to_permute)
        ) |> 
        t()
      
      result[1:num_tfs_to_modulate, ] <- modulated_tfs
    }
    
    return(result)
    
  }


#' Title
#'
#' @param num_tfs 
#' @param num_res 
#' @param num_tgs 
#' @param bernoulli_param_u 
#' @param bernoulli_param_v 
#'
#' @return
#' @export
#'
#' @examples
simulate_u_and_v <- 
  function(
    num_tfs, 
    num_res, 
    num_tgs, 
    bernoulli_param_u = 0.25, 
    bernoulli_param_v = 0.25
  ) { 
  
  u <- 
    simulate_u(
      num_tfs = num_tfs, 
      num_res = num_res, 
      bernoulli_param = bernoulli_param_u
    )
  
  v <- 
    simulate_v(
      num_res = num_res, 
      num_tgs = num_tgs, 
      bernoulli_param = bernoulli_param_v
    )
  
  result <- list(u = u, v = v)
  return(result)
  }

# returns a tidy edge list for the weight matrix u 
#' Title
#'
#' @param u 
#'
#' @return
#' @export
#'
#' @examples
tidy_u <- function(u) {
  u_tidy <- 
    u |> 
    as_tibble(rownames = "tf") |> 
    pivot_longer(
      cols = -tf, 
      names_to = "re", 
      values_to = "weight"
    ) |> 
    mutate(
      edge_type = "u", 
      edge_id = paste0("u",as.character(1:n()))
    )
  
  return(u_tidy)
}

# returns a tidy edge list for the weight matrix v
#' Title
#'
#' @param v 
#'
#' @return
#' @export
#'
tidy_v <- function(v) { 
  v_tidy <- 
    v |> 
    as_tibble(rownames = "re") |> 
    pivot_longer(
      cols = -re, 
      names_to = "tg", 
      values_to = "weight"
    ) |> 
    mutate(
      edge_type = "v", 
      edge_id = paste0("v", as.character(1:n()))
    )
  
  return(v_tidy)
}

#' Title
#'
#' @param u 
#' @param v 
#'
#' @return
#'
build_graph <- function(u, v) { 
  
  # create tibble listing all nodes 
  ## avoid collision of TF protein and gene names 
  if (!missing(v)) { 
    v$tg <- paste0(v$tg, "_gene")
  }
  
  if (missing(v)) { 
    u_nodes <- 
      u |> 
      dplyr::select(tf, re) |> 
      pivot_longer(
        cols = everything(), 
        names_to = "node_type", 
        values_to = "node_name"
      ) |> 
      distinct(node_name, node_type)
    
    v_nodes <- NULL
    
    edges <- 
      u |> 
      select(tf, re, contains("weight")) |> 
      rename(from = tf, to = re) |> 
      filter(abs(weight) > 0)
    
  } else { 
    u_nodes <- 
      u |> 
      distinct(tf) |> 
      rename(node_name = tf) |> 
      mutate(node_type = "tf")
    
    v_nodes <- 
      v |> 
      distinct(tg) |> 
      rename(node_name = tg) |> 
      mutate(node_type = "tg")
    
    u_edges <- 
      u |> 
      select(tf, re, u_weight = weight) |> 
      filter(abs(u_weight) > 0)
    
    v_edges <- 
      v |> 
      select(re, tg, v_weight = weight) |> 
      filter(abs(v_weight) > 0)
    
    edges <- 
      u_edges |> 
      left_join(v_edges, by = "re") |> 
      drop_na(tg) |> 
      transmute(
        from = tf, 
        to = tg, 
        weight = u_weight * v_weight, 
        edge_id = re
      )
  }
  
  nodes <- bind_rows(u_nodes, v_nodes)
  
  result <-
    tidygraph::tbl_graph(
      nodes = nodes,
      edges = edges,
      node_key = "node_name"
    )
  
  return(result)
}

build_graph_matrix <- function(u, v) {
  
  # create tibble listing all nodes 
  ## avoid collision of TF protein and gene names 
  if (!missing(v)) { 
    colnames(v) <- paste0(colnames(v), "_gene")
  }
  
  if (missing(v)) { 
    u_nodes <- 
      u |> 
      tidy_u() |> 
      dplyr::select(tf, re) |> 
      pivot_longer(
        cols = everything(), 
        names_to = "node_type", 
        values_to = "node_name"
      ) |> 
      distinct(node_name, node_type)
    
    v_nodes <- NULL
    
    edges <- 
      u |> 
      tidy_u() |> 
      select(tf, re, contains("weight")) |> 
      rename(from = tf, to = re) |> 
      filter(abs(weight) > 0)
    
  } else { 
    uv <- u %*% v
    
    
    
    u_nodes <- 
      tibble::tibble(node_name = rownames(uv)) |>  
      mutate(node_type = "tf")
    
    v_nodes <- 
      tibble::tibble(node_name = colnames(uv)) |>  
      mutate(node_type = "tg")
    
    # gives the sum of all RE->TG edges for each TF-> TG pair instead of 
    # maintaining each RE (as an edge) separately. 
    edges <- 
      uv |> 
      tidy_u() |> 
      rename(tg = re) |> 
      select(-edge_type, -edge_id) |> 
      filter(abs(weight) > 0)
    
  }
  
  nodes <- bind_rows(u_nodes, v_nodes)
  
  result <-
    tidygraph::tbl_graph(
      nodes = nodes,
      edges = edges,
      node_key = "node_name"
    )
  
  return(result)
  }

#' Title
#'
#' @param u 
#' @param v 
#'
#' @return
#' @export
#'
#' @examples
build_graph_legacy <- function(u, v) { 
  
  # create tibble listing all nodes 
  ## avoid collision of TF protein and gene names 
  if (!missing(v)) { 
    v$tg <- paste0(v$tg, "_gene")
  }
  
  if (!missing(u)) { 
  u_nodes <- 
    u |> 
    dplyr::select(tf, re) |> 
    pivot_longer(
      cols = everything(), 
      names_to = "node_type", 
      values_to = "node_name"
    ) |> 
    distinct(node_name, node_type)
  } else { 
    u_nodes <- NULL
  }
  
  if (!missing(v)) { 
    v_nodes <- 
      v |> 
      dplyr::select(re, tg) |> 
      pivot_longer(
        cols = everything(), 
        names_to = "node_type", 
        values_to = "node_name"
      ) |> 
      distinct(node_name, node_type)
  } else { 
    v_nodes <- NULL
  }
  
  nodes <- bind_rows(u_nodes, v_nodes)
  
  # create tibbles representing the edges of the graph
  if (!missing(u)) { 
    u_edges <- 
      u |> 
      select(tf, re, contains("weight")) |> 
      rename(from = tf, to = re) 
  } else { 
    u_edges <- NULL
  }
  
  if (!missing(v)) { 
    v_edges <- 
      v |> 
      select(re, tg, contains("weight")) |> 
      rename(from = re, to = tg)
    
  } else { 
    v_edges <- NULL
  }
  
  result <-
    tidygraph::tbl_graph(
      nodes = nodes,
      edges = dplyr::bind_rows(u_edges, v_edges),
      node_key = "node_name"
    )

  return(result)
}
  
#' Title
#'
#' @param u_1 
#' @param u_2 
#' @param from_col 
#' @param to_col 
#' @param metadata_cols 
#'
#' @return
#'
build_difference_graph <- function(u_1, u_2, v_1, v_2, metadata_cols) { 
  
  if (missing(u_1) | missing(u_2)) { 
    stop("both u_1 and u_2 must be supplied.")
  }
  
  # if you have both v_1 and v_2
  if (!missing(v_1) & !missing(v_2)) { 
    graph_1 <- build_graph(u = u_1, v = v_1)
    graph_2 <- build_graph(u = u_2, v = v_2)
    
    edges_1 <- 
      graph_1 %E>% 
      as_tibble() |> 
      transmute(
        from = as_tibble(graph_1)$node_name[from], 
        to = as_tibble(graph_1)$node_name[to], 
        weight_1 = weight, 
        edge_id
      )
    
    edges_2 <- 
      graph_2 %E>% 
      as_tibble() |> 
      transmute(
        from = as_tibble(graph_2)$node_name[from], 
        to = as_tibble(graph_2)$node_name[to], 
        weight_2 = weight, 
        edge_id
      )
    
    edges_combined <- 
      full_join(edges_1, edges_2, by = c("from", "to", "edge_id")) |> 
      relocate(starts_with("weight"), .after = last_col()) |> 
      mutate(
        weight_1 = replace_na(weight_1, replace = 0), 
        weight_2 = replace_na(weight_2, replace = 0), 
        weight_diff = weight_2 - weight_1
      ) |> 
      filter(abs(weight_diff) > 0) |> 
      rename(weight = weight_diff, tf = from, re = to)
    
    result <- 
      edges_combined |> 
      build_graph() |> 
      mutate(node_type = if_else(node_type == "re", "tg", node_type))
    
    
    # if you only have one of v_1 or v_2
  } else if (xor(missing(v_1), missing(v_2))) {
    stop("If v_1 or v_2 is supplied, the other must be supplied as well.")
    
    # you have neither v_1 nor v_2
  } else { 
    
    u_1 <- rename(u_1, weight_1 = weight)
    u_2 <- rename(u_2, weight_2 = weight)
    
    metadata_colnames <- 
      u_1 |> 
      select(tf, re, {{ metadata_cols }}) |> 
      colnames()
    
    u_combined <- 
      full_join(u_1, u_2, by = metadata_colnames) |> 
      relocate(starts_with("weight"), .after = last_col())
    
    result <- 
      u_combined |> 
      mutate(
        weight_1 = replace_na(weight_1, replace = 0), 
        weight_2 = replace_na(weight_2, replace = 0), 
        weight_diff = weight_2 - weight_1
      ) |> 
      filter(abs(weight_diff) > 0) |> 
      rename(weight = weight_diff) |> 
      build_graph() |> 
      activate("edges") |> 
      rename(weight_diff = weight) |> 
      activate("nodes")
  }
    
  return(result)
  
}

#' Title
#'
#' @param graph_tbl 
#' @param weights 
#'
#' @return
#'
calculate_tf_centralities <- function(graph_tbl, weights = NULL) {
  result <- 
    graph_tbl |> 
    mutate(
      degree = centrality_degree(weights = {{ weights }}), 
      #betweenness = centrality_betweenness(weights = weights), 
      #closeness = centrality_closeness(weights = weights), 
      #eigen = centrality_eigen(weights = weights), 
      #hub = centrality_hub(weights = weights, scale = FALSE), 
      #katz = centrality_katz(), 
      #pagerank = centrality_pagerank(weights = weights), 
      #authority = centrality_authority(weights = NULL)
    ) |> 
    filter(node_type == "tf") |> 
    as_tibble() |> 
    arrange(-degree)
  
  return(result)
}


make_predictions <- 
  function(centralities_1, centralities_2, centralities_diff) { 
    # calculate predictions 
    if (missing(centralities_diff)) { 
      
      if (missing(centralities_1) | missing(centralities_2)) { 
        stop(
          "Either centralities_diff or both centralities_1 
          and centralities_2 must be specified"
        )
      } else {
        prediction <- 
          centralities_1 |> 
          rename(degree_1 = degree) |> 
          left_join(
            centralities_2 |> 
              rename(degree_2 = degree), 
            by = c("node_type", "node_name")
          ) |> 
          transmute(
            node_name, 
            degree_diff = degree_2 - degree_1
          ) |> 
          arrange(-abs(degree_diff)) |> 
          #mutate(prediction = c(rep("yes", 3), rep("no", nrow(centralities_1) - 3))) |> 
          mutate(
            .pred = scales::rescale(abs(degree_diff)), 
            p_value = 1 - .pred
          )
      }
    } else { 
      prediction <- 
        centralities_diff |> 
        arrange(-degree) |> 
        #mutate(prediction = c(rep("yes", 3), rep("no", nrow(centralities_diff) - 3))) |> 
        mutate(.pred = scales::rescale(degree)) 
    }
    
    return(prediction)
  }
  

#' Title
#'
#' @param prediction Needs to have node_name column. 
#' @param truth Needs to have node_name column.
#'
#' @return
#'
assess_differential_prediction <- function(prediction, truth) { 
  dat <- 
    prediction |> 
    left_join(truth, by = "node_name") |>
    mutate(is_differential = factor(is_differential, levels = c("yes", "no"))) 
  
  roc_auc <- 
    dat |> 
    yardstick::roc_auc(.pred, truth = is_differential) |> 
    pull(.estimate)
  
  roc <- 
    dat |> 
    yardstick::roc_curve(.pred, truth = is_differential)
  
  roc_plot <- 
    roc |> 
    autoplot() + 
    labs(caption = paste0("ROC AUC: ", round(roc_auc, 3)))
  
  return(list(roc_auc = roc_auc, roc = roc, roc_plot = roc_plot))
}



#' Title
#'
#' @param num_tfs 
#' @param num_res 
#' @param num_tgs 
#' @param permutation_prop 
#' @param prop_edges_to_enrich 
#' @param prop_edges_to_deplete 
#' @param difference_method 
#'
#' @return
#' @export
#'
run_full_simulation <- 
  function(
    num_tfs = 10, 
    num_res = 300, 
    num_tgs = 1000, 
    num_tfs_to_modulate = floor(0.3 * num_tfs), 
    permutation_prop = 0.05, 
    prop_edges_to_enrich = 0.3, 
    prop_edges_to_deplete = 0, 
    prop_edges_to_permute = 0,
    difference_method = c("centrality", "edge")
  ) { 
    u <- simulate_u(num_tfs = num_tfs, num_res = num_res)
    
    v_tidy <-
      simulate_v(num_res = num_res, num_tgs = num_tgs) |>
      tidy_v()
    
    # create null distribution for each TF using u_tidy (and v_tidy)?
    
    u_1 <- 
      u |> 
      permute_edges(permutation_prop = permutation_prop) |> 
      tidy_u() |> 
      filter(abs(weight) > 0) 
    
    u_2 <- 
      u |> 
      add_signal(
        num_tfs_to_modulate = num_tfs_to_modulate, 
        prop_edges_to_enrich = prop_edges_to_enrich, 
        prop_edges_to_deplete = prop_edges_to_deplete, 
        prop_edges_to_permute = prop_edges_to_permute
      ) |> 
      permute_edges(permutation_prop = permutation_prop) |> 
      tidy_u() |> 
      filter(abs(weight) > 0)
    
    # build graphs 
    graph_1 <- build_graph(u = u_1, v = v_tidy)
    graph_2 <- build_graph(u = u_2, v = v_tidy)
    
    graph_diff <- 
      build_difference_graph(
        u_1 = u_1, 
        u_2 = u_2, 
        v_1 = v_tidy,
        v_2 = v_tidy
      )
    
    # calculate centralities 
    centralities_1 <- 
      graph_1 |> 
      calculate_tf_centralities()
    
    centralities_2 <- 
      graph_2 |> 
      calculate_tf_centralities()
    
    centralities_diff <- 
      graph_diff |> 
      calculate_tf_centralities()
    
    # calculate predictions 
    if (difference_method == "centrality") { 
      prediction <- 
        centralities_1 |> 
        rename(degree_1 = degree) |> 
        left_join(
          centralities_2 |> 
            rename(degree_2 = degree), 
          by = c("node_type", "node_name")
        ) |> 
        transmute(
          node_name, 
          degree_diff = degree_2 - degree_1
        ) |> 
        arrange(-abs(degree_diff)) |> 
        #mutate(prediction = c(rep("yes", num_tfs_to_modulate), rep("no", nrow(centralities_1) - num_tfs_to_modulate))) |> 
        mutate(.pred = scales::rescale(abs(degree_diff)))
    } else { 
      prediction <- 
        centralities_diff |> 
        arrange(-degree) |> 
        #mutate(prediction = c(rep("yes", 3), rep("no", nrow(centralities_diff) - 3))) |> 
        mutate(.pred = scales::rescale(degree)) 
    }
    
    # evaluate predictions 
    truth <- 
      tibble(
        node_name = paste0("TF", 1:nrow(centralities_1)), 
        is_differential = 
          c(
            rep("yes", num_tfs_to_modulate), 
            rep("no", nrow(centralities_1) - num_tfs_to_modulate)
          )
      )
    
    result <- 
      prediction |> 
      assess_differential_prediction(truth = truth)
    
    result$prediction <- prediction
    
    return(result)
  }


compute_null_distributions <- 
  function(
    weight_matrix, 
    num_replicates = 100, 
    addition_prop = 0.1, 
    remove_prop = 0.1, 
    permutation_prop = 0.1, 
    addition_num = NULL, 
    remove_num = NULL, 
    permutation_num = NULL, 
    difference_method = c("centrality", "edge"), 
    weights = NULL
  ) { 
    u_list <- 
      map(
        .x = 1:num_replicates, 
        .f = ~
          add_noise(
            weight_matrix = weight_matrix, 
            addition_prop = addition_prop, 
            remove_prop = remove_prop, 
            permutation_prop = permutation_prop, 
            addition_num = addition_num, 
            remove_num = remove_num, 
            permutation_num = permutation_num
          )
      ) |> 
      map(.f = tidy_u) |> 
      map(.f = ~ filter(.x, abs(weight) > 0))
    
    u_1 <- 
      weight_matrix |> 
      add_noise(
        addition_prop = addition_prop, 
        remove_prop = remove_prop, 
        permutation_prop = permutation_prop, 
        addition_num = addition_num, 
        remove_num = remove_num, 
        permutation_num = permutation_num
      ) |> 
      tidy_u() |> 
      filter(abs(weight) > 0)
    
    if (difference_method == "centrality") { 
      centralities_1 <- 
        u_1 |> 
        build_graph() |> 
        calculate_tf_centralities(weights = weights)
      
      null_tibble <- 
        u_list |> 
        map(.f = build_graph) |> 
        map(.f = calculate_tf_centralities, weights = weights) |> 
        map2_dfr(
          .y = 1:num_replicates, 
          .f = ~ 
            make_predictions(centralities_1 = centralities_1, centralities_2 = .x) |> 
            mutate(replicate = .y)
        )  |>
        mutate(degree_diff = abs(degree_diff)) |> 
        select(-replicate, -.pred, -p_value) |> 
        group_by(node_name) |> 
        nest(data = degree_diff) |> 
        ungroup() |> 
        transmute(
          node_name,
          null_distribution = map(.x = data, .f = tibble::deframe)
        )
    } else if (difference_method == "edge") { 
      null_tibble <- 
        u_list |> 
        map(
          .f = ~ 
            build_difference_graph(
              u_1 = u_1, 
              u_2 = .x, 
              from_col = tf, 
              to_col = re, 
              metadata_cols = c(edge_type, edge_id)
            )
        ) |> 
        map(.f = calculate_tf_centralities) |> 
        map2_dfr(
          .y = 1:num_replicates, 
          .f = ~ 
            make_predictions(centralities_diff = .x) |> 
            mutate(replicate = .y)
        ) |> 
        select(node_name, degree) |> 
        group_by(node_name) |> 
        nest(data = degree) |> 
        ungroup() |> 
        transmute(
          node_name,
          null_distribution = map(.x = data, .f = tibble::deframe)
        )
    } else { 
      stop("difference_method must be centrality or edge.")
    }
    
    
    
    return(null_tibble)
  }


compare_networks <- 
  function(
    edge_list_1, 
    edge_list_2, 
    difference_method = c("centrality", "edge"), 
    null_distributions
  ) { 
    
    u_1 <- 
      edge_list_1 |> 
      filter(abs(weight) > 0)
    
    u_2 <- 
      edge_list_2 |> 
      filter(abs(weight) > 0)
    
    if (difference_method == "centrality") { 
    graph_1 <- 
      u_1 |> 
      build_graph()
    
    graph_2 <- 
      u_2 |> 
      build_graph()
    
    centralities_1 <- calculate_tf_centralities(graph_tbl = graph_1)
    centralities_2 <- calculate_tf_centralities(graph_tbl = graph_2)
    
    predictions <- 
      make_predictions(
        centralities_1 = centralities_1, 
        centralities_2 = centralities_2
      )
    
    result <- 
      predictions |> 
      left_join(null_distributions, by = "node_name") |> 
      mutate(
        p_value_null = 
          map2_dbl(
            .x = abs(degree_diff), 
            .y = null_distribution, 
            .f = ~ mean(.x < abs(.y))
          ), 
        p_value_null_adj = p.adjust(p_value_null, method = "fdr")
      ) |> 
      select(-.pred, -null_distribution)
    
    } else if (difference_method == "edge") { 
      diff_graph <- 
        build_difference_graph(
          u_1 = u_1, 
          u_2 = u_2, 
          from_col = tf, 
          to_col = re, 
          metadata_cols = c(edge_type, edge_id)
        )
      
      centralities_diff <- calculate_tf_centralities(graph_tbl = diff_graph)
      
      predictions <- 
        make_predictions(centralities_diff = centralities_diff) |> 
        mutate(p_value = 1 - .pred)
      
      result <- 
        predictions |> 
        left_join(null_distributions, by = "node_name") |> 
        mutate(
          p_value_null = 
            map2_dbl(
              .x = degree, 
              .y = null_distribution, 
              .f = ~ mean(.x < .y)
            ), 
          p_value_null_adj = p.adjust(p_value_null, method = "fdr")
        ) |> 
        select(-.pred, -null_distribution)
    }
    
    return(result)
    
  }

plot_roc_curve <- function(predictions, pred_col, truth) { 
  truth <- 
    truth |> 
    mutate(is_differential = factor(is_differential, levels = c("yes", "no")))
  
  roc_auc <- 
    predictions |> 
    left_join(truth, by = "node_name") |> 
    yardstick::roc_auc({{ pred_col }}, truth = is_differential) |> 
    pull(.estimate) |> 
    round(3)
  
  roc_curve <- 
    predictions |> 
    left_join(truth, by = "node_name") |> 
    yardstick::roc_curve({{ pred_col }}, truth = is_differential) |> 
    autoplot() + 
    labs(caption = paste0("ROC AUC: ", roc_auc))
  
  return(roc_curve)
}

plot_null_distributions <- function(null_distributions, predictions, tf_numbers = 1:10) { 

  plot_tib <- 
    predictions |> 
    left_join(null_distributions, by = "node_name") |> 
    mutate(
      node_name = factor(node_name, levels = paste0("TF", 1:nrow(null_distributions)))
      ) |> 
    arrange(node_name) |> 
    slice(tf_numbers)
  
  plot_list <- 
    map2(
      .x = plot_tib$null_distribution,
      .y = plot_tib$degree, 
      .f = ~ ggplot(aes(x = .x), data = NULL) + 
        geom_histogram(bins = 30) + 
        geom_vline(xintercept = .y) + 
        theme_bw()
    ) 
  
  return(plot_list)
}


rename_u <- function(u, tf_names, re_names) { 
  rownames(u) <- tf_names
  colnames(u) <- re_names
  return(u)
}

rename_v <- function(v, re_names, tg_names) { 
  rownames(v) <- re_names
  colnames(v) <- tg_names
  return(v)
}

convert_u_to_edge_list <- function(u) { 
  result <- 
    u |> 
    as.matrix() |> 
    tidy_u() |> 
    filter(abs(weight) > 1e-6)
  
  return(result)
}


convert_v_to_edge_list <- function(v) { 
  result <- 
    v |> 
    as.matrix() |> 
    tidy_v() |> 
    filter(abs(weight) > 1e-6)
  
  return(result)
}

clean_u <- function(u, tf_names, re_names) { 
  result <- 
    u |> 
    rename_u(tf_names = tf_names, re_names = re_names) |> 
    convert_u_to_edge_list()
  return(result)
}

clean_v <- function(v, re_names, tg_names) { 
  result <- 
    v |> 
    rename_v(re_names = re_names, tg_names = tg_names) |> 
    convert_v_to_edge_list()
  return(result)
}


# TO DO: document
normalize_centralities <- 
  function(diff_centralities, cluster_1_centralities, cluster_2_centralities) { 
    cluster_1_centralities <- 
      cluster_1_centralities |> 
      rename(degree_1 = degree) |> 
      select(node_name, degree_1)
    
    cluster_2_centralities <- 
      cluster_2_centralities |> 
      rename(degree_2 = degree) |> 
      select(node_name, degree_2)
    
    result <- 
      diff_centralities |> 
      select(node_name, degree) |> 
      left_join(cluster_1_centralities, by = "node_name") |> 
      left_join(cluster_2_centralities, by = "node_name") |> 
      mutate(
        mean_centrality = (degree_1 + degree_2) / 2.0,  
        normalized_degree = degree / mean_centrality
      ) |> 
      select(node_name, normalized_degree) |> 
      arrange(-normalized_degree)
    
    return(result)
  }

normalize_predictions <- normalize_centralities

# TO DO: write a function to threshold the differential edges based on size

compare_differential_cluster_networks <- 
  function(diff_graph, weights = NULL, cluster_1_centralities, cluster_2_centralities) { 
    result <- 
      diff_graph |> 
      calculate_tf_centralities(weights = {{ weights }}) |> 
      normalize_centralities(
        cluster_1_centralities = cluster_1_centralities, 
        cluster_2_centralities = cluster_2_centralities
      ) |> 
      mutate(rank = rank(-normalized_degree))
    
    return(result)
  
  }

