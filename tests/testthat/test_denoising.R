context("denoising")

test_that('denoising runs and returns expected output', {
  skip_on_cran()
  options(mc.cores = 2)
  set.seed(42)
  vst_out <- vst(pbmc, return_cell_attr = TRUE)
  y_smooth <- smooth_via_pca(vst_out$y, do_plot = FALSE)
  expect_equal(c(9727, 402), dim(y_smooth))
  expect_equal(c(-0.02481276, -0.26294031, -0.63786577,
                 0.13792283, -0.31721807, 0.14007577,
                 1.77285804, -0.40266367, 2.53638759), as.numeric(y_smooth[1:3, 1:3]))
  umi_denoised <- denoise(vst_out)
  expect_equal(c(0, 0, 0, 0, 0, 2, 1, 0, 11), as.numeric(umi_denoised[1:3, 1:3]))
})
