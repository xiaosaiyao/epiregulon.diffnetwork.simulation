library(Matrix)

# TF expression X_ct simulated from Poisson
# nTFs is total number of TFs to simulate
# nCells is total number of cells to simulate
# groups is total number of groups or "cell types" to simulate, each will have
# if groups = 2, the number of cells per group is nCells/2
# if rates (Poisson parameter) can be specified to a groups x nTFs data frame that specifies the
# rate values for each TF in each group. 
# otherwise, defaults to making a single TF differ drastically in the rate
# parameter for Group 1, and the rest the same
# plz make n_cells divisible by groups
simulateTF <- function(n_tfs = 5, n_cells = 100, groups = 2, rates = NULL,
                       baseline_rate = 0.5, perturbed_rate = 2) {
  if (is.null(rates)) {
    rates <- matrix(
      baseline_rate,
      nrow = groups, ncol = n_tfs,
      dimnames = list(paste0('Group', 1:groups),
                      paste0('TF', 1:n_tfs))
    )
    if (groups > 1) {
      # create artificial differences 
      rates[1,1] <- perturbed_rate
    }
  }
  # expand rates matrix to generate for each cell
  n_cells_per_group <- n_cells/groups
  rates_cells <- matrix(nrow = n_cells, ncol = n_tfs)
  for (g in 1:groups) {
    start_idx <- (g-1)*n_cells_per_group+1
    end_idx <- start_idx + n_cells_per_group-1
    rates_cells[start_idx:end_idx,] <- rep(rates[g,], each = n_cells_per_group)
  }
  # create TF expression matrix
  expr <- matrix(
    rpois(n_cells*n_tfs, rate =  rates_cells),
    nrow = n_cells, ncol = n_tfs,
    dimnames = list(
      paste0('cell',1:n_cells),
      paste0('tf',1:n_tfs)
    )
  )
  return(expr)
}


# simulate lambda_cr if not inputting from peakVI output
# either beta- or gamma-distributed
simulateLambda <- function(n_cells = 100, n_peaks = 300, option='beta') {
  return(
    matrix(
      data = ifelse(
        option=='beta', rbeta(n_cells*n_peaks, 1, 1), 
        rgamma(n_cells*n_peaks, 1, 1)
      ),
      nrow = n_cells, ncol = n_peaks,
      dimnames = list(
        paste0('cell',1:n_cells),
        paste0('peak',1:n_peaks)
      )
    )
  )
}

# simulate peaks
simulatePeak <- function(n_cells = 100, n_peaks = 300, option = 'bernoulli', lambda = NULL) {
  if (is.null(lambda)) {
    lambda <- simulateLambda(n_cells, n_peaks)
  }
  return(
    matrix(
      rbinom(n_cells*n_peaks, 1, lambda),
      nrow = n_cells, ncol = n_peaks,
      dimnames = list(
        paste0('cell', 1:n_cells),
        paste0('peak', 1:n_peaks)
      )
    )
  )
}

# simulate indicator matrix if not inputting from ChIP-seq/TF motifs for
# TF>RE or inputting positional information (REs within 100kb of gene) RE>TG
# density is what fraction of the matrix is non-zero
simulateIndicator <- function(nrow, ncol, density=0.1) {
  return(
    rsparsematrix(nrow=nrow, ncol=ncol, density=density, 
                  rand.x = function(n)  rbinom(n, size=1, prob=1))
  )
}

# simulate binding affinity B
simulateB <- function(n_cells = 100, n_peaks = 300, binary = T, indicator_mat = NULL) {
  if (is.null(indicator_mat)) {
    indicator_mat <- simulateIndicator(n_cells, n_peaks, density=0.1)
  }
  B <- indicator_mat
  if (!binary) {
    idx <- which(B>0)
    B[idx] <- runif(length(idx), 0, 1)
  } 
  return(B)
}

# simulate the RE->TG weights V
simulateV <- function(n_peaks = 300, n_genes = 50, binary = T,  indicator_mat = NULL) {
  if (is)
}

simulateGene <- function()

# Target gene expression
# n is the number of genes to simulate the expression of
simulateTG <- function(RE, n, v = NULL) {
  lambdas <- RE$lambdas
  nCells <- nrow(lambdas)
  nREs <- ncol(lambdas)
  if (!is.null(v)) {
    if (length(dim(v))!=2) {
      stop('if specifying, v must be a matrix of weights of RE x TG')
    }
    if (ncol(v) != n) {
      stop('Length of v columns must be equal to number of generated TGs n')
    }
    if (nrow(v) != nREs) {
      stop('Length of v rows must be equal to number of REs in RE object')
    }
  } else {
    v <- matrix(runif(nREs*n, 0, 0.1), nrow = nREs, ncol = n,
                dimnames = list(paste0('Peak', seq(1, nREs)),
                                paste0('Gene', seq(1, n))))
  }
  means <- lambdas%*%v
  lambdas <- matrix(rgamma(nCells*n, shape = exp(means)), nrow = nCells, ncol = n,
                    dimnames = list(paste0('Cell', seq(1, nCells)),
                                    paste0('Gene', seq(1, n))))
  y_cg <- matrix(rpois(n = nCells*n, lambda = lambdas), nrow = nCells, ncol = n,
                 dimnames = list(paste0('Cell', seq(1, nCells)),
                                 paste0('Gene', seq(1, n))))
  dropout <- matrix(rbinom(n = nCells*n, size = 1, prob = 0.5), nrow = nCells, ncol = n)
  y_cg[dropout==1] <- 0
  return(list(weights = v, lambdas = lambdas, counts = y_cg))
}
