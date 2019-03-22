library(upbm)
context("upbm-classes")

## ---------------

data(hoxc9alexa, package = "upbmData")
data(pbm_8x60k_v1, package = "upbmAux")
data(refcy3_8x60k_v1, package = "upbmAux")

## design examples
simple_design <-
    data.frame(probeID = LETTERS[1:10],
               Sequence = rep("ACGT", 10))
dupd_design <- 
    data.frame(probeID = c(LETTERS[1:10], LETTERS[1:5]),
               Type = c(rep("keep", 10), rep("skip", 5)), 
               Sequence = rep("ACGT", 15))

## SummarizedExperiment examples
se_assays <- SimpleList(a1 = DataFrame(a = 1:3, b = 2:4))
se_rowdata <- DataFrame(id = LETTERS[1:3])
se_coldata <- DataFrame(lab = letters[1:2])

se_empty <- SummarizedExperiment()
se_partial <- SummarizedExperiment(assays = se_assays)
se_complete <- SummarizedExperiment(assays = se_assays,
                                    rowData = se_rowdata,
                                    colData = se_coldata)

## ---------------

test_that("PBMDesign constructor works", {
    ## test behavior for empty constructor
    expect_warning(pd_empty <- PBMDesign())
    expect_is(pd_empty, "PBMDesign")

    ## test behavior for wrong input
    expect_warning(pd_empty <- PBMDesign(1:3), "Ignoring specified")
    expect_is(pd_empty, "PBMDesign")

    ## test behavior for bad input, data.frame no cols
    expect_error(PBMDesign(data.frame()))

    ## test behavior for expected input, data.frame with cols
    expect_silent(pd_df <- PBMDesign(simple_design))
    expect_is(pd_df, "PBMDesign")

    ## test behavior for other parameters
    expect_silent({
        pd_df <- PBMDesign(simple_design,
                           probeFilter = list(probeID = function(x) {
                               rep(TRUE, length(x))
                           }),
                           probeTrim = c(1, 36))
    })
    expect_equal(pd_df@design, simple_design)
    expect_equal(pd_df@probeTrim, c(1, 36))

    ## test behavior with duplicated probeIDs
    expect_error(PBMDesign(dupd_design),
                 "probe IDs must be unique after filtering")
    expect_silent(pd_dupd <- PBMDesign(dupd_design,
                                       probeFilter = list(Type = function(x) {
                                           x == "keep"
                                       }))
                  )
    expect_is(pd_dupd, "PBMDesign")
    expect_equal(pd_dupd@design, dupd_design)
})

## ---------------

test_that("PBMExperiment constructor works", {
    ## test behavior for empty constructor
    expect_silent(pe <- PBMExperiment())
    expect_is(pe, "PBMExperiment")

    ## test behavior for constructing from SummarizedExperiment - empty
    expect_silent(pe_empty <- PBMExperiment(se_empty))
    expect_is(pe_empty, "PBMExperiment")

    ## test behavior for constructing from SummarizedExperiment - with assays only
    expect_silent(pe_partial <- PBMExperiment(se_partial))
    expect_is(pe_partial, "PBMExperiment")
    expect_equal(assays(pe_partial), assays(se_partial))

    ## test behavior for constructing with SE components
    expect_silent(pe_sespecify <- PBMExperiment(assays = se_assays))
    expect_is(pe_sespecify, "PBMExperiment")
    expect_equal(pe_sespecify, pe_partial)

    ## test behavior for constructing with only PBDesign
    ## - if PBMDesign specified, need to have non-empty assays
    expect_error(pe_pd <- PBMExperiment(pbmDesign = pbm_8x60k_v1))

    ## test PBMExperiment -> SummarizedBechmark method
    expect_silent(se_from_pe <- SummarizedExperiment(hoxc9alexa))
    expect_is(se_from_pe, "SummarizedExperiment")
    expect_equal(assayNames(se_from_pe), assayNames(hoxc9alexa))
    expect_equal(dim(se_from_pe), dim(hoxc9alexa))

    ## test behavior for constructing with SE and PBDesign
    expect_warning({pe_pd <- PBMExperiment(se_from_pe, pbmDesign = pbm_8x60k_v1)},
                   "will be overwritten")
    expect_is(pe_pd, "PBMExperiment")
    expect_equal(pe_pd, hoxc9alexa)

    ## test behavior for constructing with SE and PBMDesign components
    expect_silent(pe_pdspecify <- PBMExperiment(se_from_pe,
                                                probeFilter = pbm_8x60k_v1@probeFilter,
                                                probeTrim = pbm_8x60k_v1@probeTrim,
                                                probeCols = colnames(pbm_8x60k_v1@design)))
    expect_is(pe_pdspecify, "PBMExperiment")
    expect_equal(pe_pdspecify, hoxc9alexa)

    ## test behavior for constructuing with SE components and PBMdesign
    expect_silent(pe_sespecifypd <- PBMExperiment(assays = assays(hoxc9alexa),
                                                  pbmDesign = pbm_8x60k_v1))
    expect_is(pe_sespecifypd, "PBMExperiment")
})

## ---------------

test_that("Setters and getters for new classes work", {
    test_pd <- pbm_8x60k_v1

    ## test all getters
    expect_equal(probeFilter(test_pd), test_pd@probeFilter)
    expect_equal(probeTrim(test_pd), test_pd@probeTrim)
    expect_equal(design(test_pd), test_pd@design)

    expect_equal(probeFilter(hoxc9alexa), hoxc9alexa@probeFilter)
    expect_equal(probeTrim(hoxc9alexa), hoxc9alexa@probeTrim)
    expect_equal(probeCols(hoxc9alexa), hoxc9alexa@probeCols)

    expect_equal(test_pd, PBMDesign(hoxc9alexa))

    ## test probeFilter replace
    expect_silent({
        probeFilter(test_pd) <-
            list(probeID = function(x) {
                grepl("^dBr", x)
            },
            Row = function(x) {
                rep(TRUE, length(x))
            })
    })
    expect_equal(names(probeFilter(test_pd)), c("probeID", "Row"))

    ## test probeTrim replace
    expect_silent({ probeTrim(test_pd) <- vector("numeric") })
    expect_equal(probeTrim(test_pd), vector("numeric"))

    ## test design replace
    expect_silent({ design(test_pd) <- cbind(design(test_pd), newcol = 1) })
    expect_equal(names(design(test_pd)), c(names(design(pbm_8x60k_v1)), "newcol"))

    ## test PBMDesign replace
    test_pe <- hoxc9alexa
    expect_silent({ PBMDesign(test_pe) <- test_pd })
    expect_equal(probeCols(test_pe), colnames(design(test_pd)))
    expect_equal(probeFilter(test_pe), probeFilter(test_pd))
    expect_equal(probeTrim(test_pe), probeTrim(test_pd))
})

