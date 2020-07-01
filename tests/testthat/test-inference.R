library(upbm)
context("upbm-tests")

## ---------------

data(hoxc9alexa, package = "upbmData")
data(pbm_8x60k_v1, package = "upbmAux")
data(refcy3_8x60k_v1, package = "upbmAux")

## subset data
adata <- hoxc9alexa[, (colData(hoxc9alexa)$condition %in% c("HOXC9-REF", "HOXC9-R193K")) &
                      (colData(hoxc9alexa)$id %in% c(226, 249)) &
                      (colData(hoxc9alexa)$pmt == 450)]

adata <- normalizeWithinReplicates(adata)
adata <- normalizeAcrossReplicates(adata)

## ---------------

test_that("probe and kmer fits and tests return expected object types", {
    ## test fit probe models
    expect_silent(pf <- probeFit(adata))
    expect_s4_class(pf, "PBMExperiment")
    expect_equal(nrow(pf), nrow(pbmFilterProbes(adata)))
    expect_equal(ncol(pf), length(unique(colData(adata)$condition)))
    expect_equal(assayNames(pf), c("beta", "sd", "df"))

    ## test fit kmer models
    nkmers <- 100
    expect_silent(kf <- kmerFit(pf, uniqueKmers(8)[1:nkmers], method = "dl2"))
    expect_s4_class(kf, "SummarizedExperiment")
    expect_equal(nrow(kf), nkmers)

    expect_silent(kf2 <- kmerFit(pf, uniqueKmers(8)[1:nkmers], method = "dl"))
    expect_s4_class(kf2, "SummarizedExperiment")
    expect_equal(nrow(kf), nkmers)

    ## test kmer contrast test
    expect_silent(kdelta <- kmerTestContrast(kf))
    expect_s4_class(kdelta, "SummarizedExperiment")
    expect_true(all(paste0("contrast", c("Z", "P", "Q")) %in% assayNames(kdelta)))

    ## test kmer affinity test
    expect_silent(kaff <- kmerTestAffinity(kf))
    expect_s4_class(kaff, "SummarizedExperiment")
    expect_true(all(paste0("affinity", c("Z", "P", "Q")) %in% assayNames(kaff)))

    ## test kmer specificity test
    expect_silent(kspec <- kmerTestSpecificity(kf))
    expect_s4_class(kspec, "SummarizedExperiment")
    expect_true(all(paste0("specificity", c("Z", "P", "Q")) %in% assayNames(kspec)))
})
