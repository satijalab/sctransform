context("Rcpp utility functions")

test_that('row_mean_grouped runs and returns expected output', {
  skip_on_cran()
  suppressWarnings(RNGversion(vstr = "3.5.0"))
  set.seed(42)

  grouping <- as.factor(sample(c('a','b','c'), size = ncol(pbmc), replace = TRUE))
  means <- sctransform:::row_mean_grouped(pbmc, grouping)
  means_agg <- t(apply(pbmc, 1, function(x) {
    aggregate(x = x, by = list(group = grouping), FUN = mean)$x
  }))
  colnames(means_agg) <- levels(grouping)

  expect_equal(means, means_agg)

  gmeans <- sctransform:::row_gmean_grouped(pbmc, grouping, eps = 1)
  gmeans_agg <- sapply(levels(grouping), function(g) {
    sctransform:::row_gmean(pbmc[, grouping == g])
  })

  expect_equal(gmeans, gmeans_agg)
})
