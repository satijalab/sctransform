context("vst function")

test_that('vst runs and returns expected output', {
  skip_on_cran()
  suppressWarnings(RNGversion(vstr = "3.5.0"))
  set.seed(42)
  vst_out <- vst(pbmc, return_gene_attr = TRUE, return_cell_attr = TRUE)
  expect_equal(c(910, 283), dim(vst_out$y))
  ga <- vst_out$gene_attr[order(-vst_out$gene_attr$residual_variance), ]
  expect_equal(c("GNLY", "NKG7", "GZMB", "LYZ", "S100A9"), rownames(ga)[1:5])
  expect_equal(c(27.8, 26.9, 18.7, 18.2, 16.8), ga$residual_variance[1:5], tolerance = 1e-01)
})

test_that('vst runs with multicore futures', {
  skip_on_cran()
  suppressWarnings(RNGversion(vstr = "3.5.0"))
  set.seed(42)

  options(future.fork.enable = TRUE)
  options(future.globals.maxSize = 10 * 1024 ^ 3)
  future::plan(strategy = 'multicore', workers = 2)

  vst_out <- vst(pbmc, return_gene_attr = TRUE, return_cell_attr = TRUE)
  expect_equal(c(910, 283), dim(vst_out$y))
  ga <- vst_out$gene_attr[order(-vst_out$gene_attr$residual_variance), ]
  expect_equal(c("GNLY", "NKG7", "GZMB", "LYZ", "S100A9"), rownames(ga)[1:5])
  expect_equal(c(27.8, 26.9, 18.7, 18.2, 16.8), ga$residual_variance[1:5], tolerance = 1e-01)
})

test_that('vst with batch variable works', {
  skip_on_cran()
  suppressWarnings(RNGversion(vstr = "3.5.0"))
  set.seed(42)

  ca <- data.frame(batch = sample(x = c('A', 'B'), size = ncol(pbmc), replace = TRUE))
  rownames(ca) <- colnames(pbmc)

  vst_out <- vst(pbmc, batch_var = 'batch', cell_attr = ca, return_gene_attr = TRUE, return_cell_attr = TRUE)
  expect_equal(c(910, 283), dim(vst_out$y))
  ga <- vst_out$gene_attr[order(-vst_out$gene_attr$residual_variance), ]
  expect_equal(c("GNLY", "NKG7", "GZMB", "LYZ", "S100A9"), rownames(ga)[1:5])
  expect_equal(c(27.11, 26.11, 17.89, 17.66, 16), ga$residual_variance[1:5], tolerance = 1e-01)
})

test_that('vst with pre-calculated cell attributes works', {
  skip_on_cran()
  suppressWarnings(RNGversion(vstr = "3.5.0"))

  ca <- data.frame(umi = colSums(pbmc))
  set.seed(42)
  vst_out <- vst(pbmc, cell_attr = ca, latent_var = 'log_umi', return_gene_attr = TRUE, return_cell_attr = TRUE)
  expect_equal(c(910, 283), dim(vst_out$y))
  ga <- vst_out$gene_attr[order(-vst_out$gene_attr$residual_variance), ]
  expect_equal(c("GNLY", "NKG7", "GZMB", "LYZ", "S100A9"), rownames(ga)[1:5])
  expect_equal(c(27.8, 26.9, 18.7, 18.2, 16.8), ga$residual_variance[1:5], tolerance = 1e-01)

  ca <- data.frame(log_umi = log10(colSums(pbmc)))
  set.seed(42)
  vst_out <- vst(pbmc, cell_attr = ca, latent_var = 'umi', return_gene_attr = TRUE, return_cell_attr = TRUE)
  expect_equal(c(910, 283), dim(vst_out$y))
  ga <- vst_out$gene_attr[order(-vst_out$gene_attr$residual_variance), ]
  expect_equal(c("GNLY", "NKG7", "LYZ", "S100A9", "GZMB"), rownames(ga)[1:5])
  expect_equal(c(26.5, 20.5, 19.6, 16.6, 16.5), ga$residual_variance[1:5], tolerance = 1e-01)

  ca <- data.frame(log_umi = log10(colSums(pbmc)))
  set.seed(42)
  vst_out <- vst(pbmc, cell_attr = ca, latent_var = 'log_umi', return_gene_attr = TRUE, return_cell_attr = TRUE)
  expect_equal(c(910, 283), dim(vst_out$y))
  ga <- vst_out$gene_attr[order(-vst_out$gene_attr$residual_variance), ]
  expect_equal(c("GNLY", "NKG7", "GZMB", "LYZ", "S100A9"), rownames(ga)[1:5])
  expect_equal(c(27.8, 26.9, 18.7, 18.2, 16.8), ga$residual_variance[1:5], tolerance = 1e-01)
})
