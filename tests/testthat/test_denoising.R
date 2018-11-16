context("denoising")

test_that('denoising runs and returns expected output', {
  skip_on_cran()
  options(mc.cores = 2)
  set.seed(42)
  vst_out <- vst(pbmc, return_cell_attr = TRUE, res_clip_range = c(-Inf, Inf))
  y_smooth <- smooth_via_pca(vst_out$y, do_plot = FALSE)
  expect_equal(c(910, 283), dim(y_smooth))
  expect_equal(c(0.06362666, -0.01528845,  0.84735003,  0.37468399,
                 -0.03071016, 1.34827719, -0.01143410, -0.29820178,
                 -0.64414609), as.numeric(y_smooth[1:3, 1:3]))
  umi_denoised <- denoise(vst_out)
  expect_equal(c(0, 1, 27, 1, 1, 36, 0, 0, 7), as.numeric(umi_denoised[1:3, 1:3]))
  umi_denoised <- denoise(vst_out, data = y_smooth)
  expect_equal(c(0, 0, 28, 0, 0, 34, 0, 0, 9), as.numeric(umi_denoised[1:3, 1:3]))
})
