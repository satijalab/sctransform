context("correcting")

test_that('correcting runs and returns expected output', {
  skip_on_cran()
  options(mc.cores = 2)
  suppressWarnings(RNGversion(vstr = "3.5.0"))
  set.seed(42)
  vst_out <- vst(pbmc, return_cell_attr = TRUE, res_clip_range = c(-Inf, Inf))
  y_smooth <- smooth_via_pca(vst_out$y, do_plot = FALSE)
  expect_equal(c(910, 283), dim(y_smooth))
  expect_equal(c(0.05809, -0.00707,  0.55978,  0.27358,
                 -0.01979, 0.83436, 0.03495, -0.09587,
                 -0.88417), as.numeric(y_smooth[1:3, 1:3]), tolerance = 1e-5)
  umi_corrected <- correct(vst_out)
  expect_equal(c(0, 1, 28, 1, 1, 37, 0, 0, 7), as.numeric(umi_corrected[1:3, 1:3]))
  umi_corrected <- correct(vst_out, data = y_smooth)
  expect_equal(c(0, 0, 30, 0, 0, 34, 0, 0, 9), as.numeric(umi_corrected[1:3, 1:3]))
})
