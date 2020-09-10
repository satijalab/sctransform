context("correcting")

test_that('correcting runs and returns expected output', {
  skip_on_cran()
  suppressWarnings(RNGversion(vstr = "3.5.0"))
  set.seed(42)
  vst_out <- vst(pbmc, return_cell_attr = TRUE, res_clip_range = c(-Inf, Inf))
  y_smooth <- smooth_via_pca(vst_out$y, do_plot = FALSE)
  expect_equal(c(910, 283), dim(y_smooth))
  expect_equal(c(0.0868, 0.0380,  0.6062,  0.3123,
                 0.0101, 0.8751, -0.0557, -0.2222,
                 -0.9911), as.numeric(y_smooth[1:3, 1:3]), tolerance = 1e-3)
  umi_corrected <- correct(vst_out)
  expect_equal(c(0, 1, 28, 1, 1, 37, 0, 0, 7), as.numeric(umi_corrected[1:3, 1:3]))
  umi_corrected <- correct(vst_out, data = y_smooth)
  expect_equal(c(0, 0, 31, 0, 0, 34, 0, 0, 8), as.numeric(umi_corrected[1:3, 1:3]))
})
