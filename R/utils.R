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
  expr <- Matrix::Matrix(
    nrow = n_cells, ncol = n_tfs,
    data = rpois(n_cells*n_tfs, lambda =  rates_cells),
    sparse=T,
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
  if (option=='beta') {
    vals <- rbeta(n=n_cells*n_peaks, shape1=1,shape2=8)
  } else if (option=='gamma') {
    vals <- rgamma(n=n_cells*n_peaks, shape=1,rate=1)
  } else {
    stop('option must be either beta or gamma')
  }
  return(
    matrix(
      data = vals, nrow = n_cells, ncol = n_peaks,
      dimnames = list(paste0('cell',1:n_cells), paste0('peak',1:n_peaks)
      )
    )
  )
}

# simulate peaks
simulatePeak <- function(n_cells = 100, n_peaks = 300, option = 'bernoulli', lambda = NULL) {
  if (is.null(lambda)) {
    if (option == 'bernoulli') {
      lambda <- simulateLambda(n_cells, n_peaks, option='beta')
    } else if (option == 'poisson') {
      lambda <- simulateLambda(n_cells, n_peaks, option='gamma')
    } else {
      stop('option must be bernoulli or poisson')
    }
  }
  vals <- rbinom(n_cells*n_peaks, 1, lambda)
  M <- Matrix::Matrix(nrow = n_cells, ncol = n_peaks,
                      data = vals, sparse=T)
  return(M)
}

# simulate indicator matrix if not inputting from ChIP-seq/TF motifs for
# TF>RE or inputting positional information (REs within 100kb of gene) RE>TG
# density is what fraction of the matrix is non-zero
# if the indicator matrix is banded, ignores density
simulateIndicator <- function(nrow, ncol, density=NULL, banded = F, band_size=NULL) {
  if (banded) {
    if (is.null(band_size)) {
      band_size <-5
    }
    M <- bandSparse(nrow=nrow, ncol=ncol, k=c(0:(band_size-1)))
  } else {
    if (is.null(density)) {
      density <- 0.1 # default density value
    }
    M <- rsparsematrix(nrow=nrow, ncol=ncol, density=density, 
                       rand.x = function(n)  rbinom(n, size=1, prob=1))
  }
  return(M)
}

# generic simulate weights function
# can be used for U, V, and B
simulateWeights <- function(nrow, ncol, binary = T, density = NULL, 
                            indicator_mat = NULL, banded = FALSE, band_size =  NULL) {
  if (is.null(indicator_mat)) {
    indicator_mat <- simulateIndicator(nrow, ncol, density=density, banded = banded)
  }
  M <- indicator_mat
  if (!binary) {
    idx <- which(M>0)
    M[idx] <- runif(length(idx), 0, 1)
  } 
  return(M)
}

# simulate binding affinity B
simulateB <- function(n_tfs = 5, n_peaks = 300, binary = T, indicator_mat = NULL) {
  return(simulateWeights(n_tfs, n_peaks, binary, indicator_mat))
}

# simulate the RE->TG weights V
simulateV <- function(n_peaks = 300, n_genes = 50, binary = T,  indicator_mat = NULL) {
  return(simulateWeights(n_peaks, n_genes, binary, indicator_mat))
}

# simulate the regulatory effect of TF t on RE r, alpha
simulateAlpha <- function(n_tfs = 5, n_peaks  = 300) {
  # simulate TF centers
  means <- rnorm(n_tfs, 2, 2)
  # simulate TF regulator effects
  alphas <- matrix(nrow = n_tfs, ncol = n_peaks)
  for (i in 1:n_tfs) {
    alphas[i,] <- rnorm(n_peaks, mean=means[i], 1)
  }
  return(alphas)
}

# in case want to simulate U separately without gene expression info
simulateU <-  function(n_tfs = 5, n_peaks = 300, alphas = NULL, B = NULL) {
  if (is.null(alphas)) {
    alphas <- simulateAlpha(n_tfs, n_peaks)
  }
  if (is.null(B)) {
    B <- simulateB(n_tfs, n_peaks)
  }
  return(alphas*B)
}


# simulate target gene expression
# X_tf is tf expression
# returns X_cg, target gene expression
simulateGene <- function(alphas, V, X_tf, B, lambdas) {
  binding_prob <-  1-exp(-X_tf*(lambdas%*%t(B)))
  binding <- Matrix::Matrix(nrow = nrow(binding_prob), ncol = ncol(binding_prob),
                            data = rbinom(n=length(binding_prob), size=1, prob=as.matrix(binding_prob)))
  X_cg <- binding %*% alphas %*% V
  return(X_cg)
}


