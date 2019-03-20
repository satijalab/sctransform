context("differential expression")

# test_that('compare expression runs and returns expected output', {
#   skip_on_cran()
#   options(mc.cores = 2)
#   set.seed(42)
#   vst_out <- vst(pbmc, return_cell_attr = TRUE)
#   # create fake clusters
#   clustering <- 1:ncol(pbmc) %/% 100
#   res <- compare_expression(x = vst_out, umi = pbmc, group = clustering, val1 = 0, val2 = 3)
#   expect_equal(c("AKAP17A", "LRBA", "SEC23A", "RRP8", "TRNT1"), rownames(res)[1:5])
#   expect_equal(c(-27.35713, -27.05464, -26.62938, -26.41430, -26.25116), res$log_fc[1:5], tolerance = 1e-05)
#   res <- compare_expression(x = vst_out, umi = pbmc, group = clustering, val1 = 0, val2 = 3, method = 't_test')
#   expect_equal(c("TMSB4X", "AKAP17A", "CALM3", "TOMM40", "HSPB11"), rownames(res)[1:5])
#   expect_equal(c(-0.6481318, -0.5870122, -0.7482577, -0.5022045, -0.5954648), res$log_fc[1:5], tolerance = 1e-05)
# })
