# TF activity simulated from Uniform(0,1)
# 1 is active, 0 is inactive
# nTFs is total number of TFs to simulate
# nCells is total number of cells to simulate
# groups is total number of groups or "cell types" to simulate, each will have
# a separate baseline TF activity profile. if groups = 2, the number of cells
# per group is nCells/2
simulateTF <- function(nTFs, nCells, groups = 2, sd = 0.1) {
  baseline_activity <- matrix(runif(nTFs*2), nrow = groups, ncol = nTFs,
                              dimnames = list(paste0('Group', 1:groups),
                                              paste0('TF', seq(1, nTFs))))
  baseline_activity_cells <- rnorm(
    nTFs*nCells,
    mean = rep(baseline_activity, each = nCells/groups),
    sd = sd
  )
  baseline_activity_cells[baseline_activity_cells<0] <- 0
  baseline_activity_cells[baseline_activity_cells>1] <- 1
  activity <- matrix(baseline_activity_cells,
                     nrow = nCells, ncol = nTFs,
                     dimnames = list(paste0('Cell', seq(1, nCells)),
                                     paste0('TF', seq(1, nTFs))))
  return(list(baseline_activity = baseline_activity,
              activity = activity,
              groups = rep(seq(1,groups), each = nCells/groups)))
}

# Regulatory element (ATAC peak) accessibility
# n is total number of ATAC peaks
# each peak is associated with at least 1 TF
# u is the matrix of weights for each TF x RE interaction
# N is total count per cell (sequencing depth)
# 
simulateRE <- function(TF, n, u = NULL, N = NULL) {
  x_ct <- TF$activity
  nTFs <- ncol(x_ct)
  nCells <- nrow(x_ct)
  if (!is.null(u)) {
    if (length(dim(u))!=2) {
      stop('if specifying, u must be a matrix of weights of TF x RE')
    }
    if (ncol(u) != n) {
      stop('Length of u columns must be equal to number of generated REs n')
    }
    if (nrow(u) != nTFs) {
      stop('Length of u rows must be equal to number of TFs in TF object')
    }
  } else {
    u <- matrix(runif(nTFs*n, 0, 0.05), nrow = nTFs, ncol = n,
                dimnames = list(paste0('TF', seq(1, nTFs)),
                                paste0('Peak', seq(1, n))))
  }
  shapes <- x_ct%*%u
  lambdas <- matrix(rgamma(nCells*n, shape = exp(shapes)), nrow = nCells, ncol = n,
                    dimnames = list(paste0('Cell', seq(1, nCells)),
                                    paste0('Peak', seq(1, n))))
  y_cr <- matrix(rpois(n = nCells*n, lambda = lambdas), nrow = nCells, ncol = n,
                 dimnames = list(paste0('Cell', seq(1, nCells)),
                                 paste0('Peak', seq(1, n))))
  dropout <- matrix(rbinom(n = nCells*n, size = 1, prob = 0.5), nrow = nCells, ncol = n)
  y_cr[dropout==1] <- 0
  return(list(weights = u, lambdas = lambdas, counts = y_cr))
}

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
