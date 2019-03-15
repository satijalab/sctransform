context("vst function")

test_that('vst runs and returns expected output', {
  skip_on_cran()
  options(mc.cores = 2)
  set.seed(42)
  vst_out <- vst(pbmc, return_gene_attr = TRUE, return_cell_attr = TRUE)
  expect_equal(c(910, 283), dim(vst_out$y))
  ga <- vst_out$gene_attr[order(-vst_out$gene_attr$residual_variance), ]
  expect_equal(c("GNLY", "NKG7", "GZMB", "LYZ", "S100A9"), rownames(ga)[1:5])
  expect_equal(c(27.8, 26.9, 18.7, 18.2, 16.8), ga$residual_variance[1:5], tolerance = 1e-01)
})
