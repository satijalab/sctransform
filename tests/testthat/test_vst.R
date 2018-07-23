context("vst function")

test_that('vst runs and returns expected output', {
  skip_on_cran()
  options(mc.cores = 2)
  set.seed(42)
  vst_out <- vst(pbmc, return_gene_attr = TRUE, return_cell_attr = TRUE)
  expect_equal(c(9727, 402), dim(vst_out$y))
  ga <- vst_out$gene_attr[order(-vst_out$gene_attr$residual_variance), ]
  expect_equal(c("PPBP", "GNLY", "IGLL5", "PF4", "GZMB"), rownames(ga)[1:5])
  expect_equal(c(40.30323, 32.55276, 20.85651, 20.02775, 18.82615), ga$residual_variance[1:5], tolerance = 1e-05)
})
