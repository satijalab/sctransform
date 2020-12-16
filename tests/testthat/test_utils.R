context("Rcpp utility functions")

test_that('row_mean_grouped runs and returns expected output', {
  skip_on_cran()
  suppressWarnings(RNGversion(vstr = "3.5.0"))
  set.seed(42)

  grouping <- as.factor(sample(c('a','b','c'), size = ncol(pbmc), replace = TRUE))
  means <- sctransform:::row_mean_grouped_dgcmatrix(matrix = pbmc, group = grouping, shuffle = FALSE)
  means_agg <- t(apply(pbmc, 1, function(x) {
    aggregate(x = x, by = list(group = grouping), FUN = mean)$x
  }))
  colnames(means_agg) <- levels(grouping)

  expect_equal(means, means_agg)

  gmeans <- sctransform:::row_gmean_grouped_dgcmatrix(matrix = pbmc, group = grouping, eps = 1, shuffle = FALSE)
  gmeans_agg <- sapply(levels(grouping), function(g) {
    sctransform:::row_gmean(pbmc[, grouping == g])
  })

  expect_equal(gmeans, gmeans_agg)
  
  # very sparse input matrix
  mat <- Matrix::rsparsematrix(100, 1000, density = 0.01)
  grouping <- as.factor(sample(c('a','b','c'), size = ncol(mat), replace = TRUE))
  means <- sctransform:::row_mean_grouped_dgcmatrix(matrix = mat, group = grouping, shuffle = FALSE)
  means_agg <- t(apply(mat, 1, function(x) {
    aggregate(x = x, by = list(group = grouping), FUN = mean)$x
  }))
  colnames(means_agg) <- levels(grouping)
  
  expect_equal(means, means_agg)
  
  mat[mat < 0] <- 0
  means <- sctransform:::row_gmean_grouped_dgcmatrix(matrix = mat, group = grouping, eps = 1, shuffle = FALSE)
  means_agg <- t(apply(mat, 1, function(x) {
    aggregate(x = x, by = list(group = grouping), FUN = function(y) {
      expm1(mean(log1p(y)))
    })$x
  }))
  colnames(means_agg) <- levels(grouping)
  
  expect_equal(means, means_agg)
})

test_that('row_nonzero_count runs and returns expected output', {
  skip_on_cran()
  suppressWarnings(RNGversion(vstr = "3.5.0"))
  set.seed(42)
  
  nzc <- sctransform:::row_nonzero_count_dgcmatrix(pbmc)
  nzc2 <- Matrix::rowSums(pbmc > 0)
  expect_equal(nzc, nzc2)
  
  grouping <- as.factor(sample(c('a','b','c'), size = ncol(pbmc), replace = TRUE))
  nzc <- sctransform:::row_nonzero_count_grouped_dgcmatrix(pbmc, grouping)
  f2 <- function(mat, grp) {
    ret <- sapply(levels(grp), function(g) {
      rowSums(mat[, grp == g, drop = FALSE] > 0)
    })
    colnames(ret) <- levels(grp)
    ret
  }
  nzc2 <- f2(pbmc, grouping)
  expect_equal(nzc, nzc2)
})
