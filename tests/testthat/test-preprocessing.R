library(upbm)
context("pre-processing")

## ---------------

data(hoxc9alexa, package = "upbmData")
data(hoxc9cy3, package = "upbmData")
data(pbm_8x60k_v1, package = "upbmAux")
data(refcy3_8x60k_v1, package = "upbmAux")

## subset data
adata <- hoxc9alexa[, (colData(hoxc9alexa)$condition %in% c("HOXC9-REF", "HOXC9-R193K")) &
                      (colData(hoxc9alexa)$id %in% c(226, 249)) &
                      (colData(hoxc9alexa)$pmt == 450)]
cdata <- hoxc9cy3[, colData(hoxc9cy3)$id_idx %in% colData(adata)$id_idx]

## ---------------

test_that("background subtraction works", {
    expect_silent(bsi <- backgroundSubtract(adata, "fore", "back"))
    expect_equal(bsi, backgroundSubtract(adata))
    expect_equal(as.matrix(assay(bsi)[1:5, ]),
                 as.matrix(assay(adata, "fore")[1:5, ]) -
                 as.matrix(assay(adata, "back")[1:5, ]))
})

## ---------------

test_that("Cy3 empirical fitting works", {
    expect_silent(cdata_e <- cy3FitEmpirical(cdata, refcy3_8x60k_v1))
    expect_s4_class(cdata_e, "PBMExperiment")
    expect_equal(colnames(cdata_e), colnames(cdata))
    expect_equal(dim(cdata_e), dim(cdata))
    expect_true(all(c("ratio", "lowq") %in% assayNames(cdata_e)))
})

## ---------------

test_that("Cy3 model-based fitting works", {
    expect_silent(cdata_m <- cy3FitModel(cdata))
    expect_s4_class(cdata_m, "PBMExperiment")
    expect_equal(colnames(cdata_m), colnames(cdata))
    expect_equal(dim(cdata_m), dim(cdata))
    expect_true(all(c("ratio", "lowq") %in% assayNames(cdata_m)))
})

## ---------------

test_that("Cy3 normalization works", {
    expect_silent(cdata_e <- cy3FitEmpirical(cdata, refcy3_8x60k_v1))
    expect_silent(acy3 <- cy3Normalize(adata, cdata_e))
    expect_s4_class(acy3, "PBMExperiment")
    expect_equal(dim(acy3), dim(adata))
    expect_equal(assayNames(acy3)[1], "normalized")
})

## ---------------

test_that("spatial adjustment works", {
    expect_silent(asp <- spatiallyAdjust(adata))
    expect_s4_class(asp, "PBMExperiment")
    expect_equal(dim(asp), dim(pbmFilterProbes(adata)))
    expect_equal(assayNames(asp)[1], "normalized")
})

## ---------------

test_that("cross-replicate normalization works", {
    expect_silent(anwr <- normalizeWithinReplicates(adata, method = "tmm", group = "id",
                                                    stratify = "condition", baseline = "HOXC9-REF"))
    expect_s4_class(anwr, "PBMExperiment")
    expect_equal(dim(anwr), dim(adata))

    expect_silent(anwr <- normalizeWithinReplicates(adata, method = "quantile", group = "id",
                                                    stratify = "condition", baseline = "HOXC9-R193K"))
    expect_s4_class(anwr, "PBMExperiment")
    expect_equal(dim(anwr), dim(adata))
})

## ---------------

test_that("Cross-replicate normalization works", {
    expect_silent(anar <- normalizeAcrossReplicates(adata))
    expect_s4_class(anar, "PBMExperiment")
    expect_equal(dim(anar), dim(adata))
})
