context("generate function")

test_that('generate runs and returns expected output', {
  skip_on_cran()
  suppressWarnings(RNGversion(vstr = "3.5.0"))
  set.seed(42)
  vst_out <- vst(pbmc, return_cell_attr = TRUE)

  generated_data <- generate(vst_out)
  expect_equal(c(0, 0, 0, 8, 1), generated_data['ERP29', 1:5])

  genes <- sample(x = rownames(vst_out$model_pars_fit), size = 100)
  generated_data <- generate(vst_out = vst_out, genes = genes)
  expect_equal(c(100, 283), dim(generated_data))
  expect_equal(genes, rownames(generated_data))
})
