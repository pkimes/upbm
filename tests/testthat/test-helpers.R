##devtools::load_all()
library(upbm)
context("helpers")

## ---------------

data(refcy3_8x60k_v1, package = "upbmAux")
data(pbm_8x60k_v1, package = "upbmAux")
data(hoxc9alexa, package = "upbmData")

## ---------------

test_that("tidy methods can be called on PBMExperiment objects", {
    expect_silent({ tpe <- broom::tidy(hoxc9alexa) })
    expect_is(tpe, "data.frame")
    expect_true(all(names(rowData(hoxc9alexa)) %in% names(tpe)))
    expect_false(any(names(colData(hoxc9alexa)) %in% names(tpe)))
    
    expect_silent({ tpel <- broom::tidy(hoxc9alexa, long = TRUE) })
    expect_is(tpel, "data.frame")
    expect_true(all(names(rowData(hoxc9alexa)) %in% names(tpel)))
    expect_true(all(names(colData(hoxc9alexa)) %in% names(tpel)))
    
    expect_error({ broom::tidy(hoxc9alexa, "bad-assay-name", long = TRUE) },
                 "invalid")
})

## ---------------

test_that("k-mer summaries can be computed", {
    expect_silent(kmers <- uniqueKmers(8L))
    expect_is(kmers, "character")
    expect_length(kmers, 32896)
    expect_error(uniqueKmers(11L))

    expect_silent(ksum <- summarizeKmers(hoxc9alexa[1:1000, 1:2], kmers = kmers[1:10]))
    expect_equal(ncol(ksum), 2)
    expect_true(all(c("medianIntensity", "meanIntensity", "madIntensity",
                      "sdIntensity", "log2meanIntensity", "log2madIntensity",
                      "log2sdIntensity", "naProbes", "q25Intensity") %in% assayNames(ksum)))
})
