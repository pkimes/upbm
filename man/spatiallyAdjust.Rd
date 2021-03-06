% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatiallyAdjust.R
\name{spatiallyAdjust}
\alias{spatiallyAdjust}
\title{Perform spatial adjustment}
\usage{
spatiallyAdjust(
  pe,
  assay = SummarizedExperiment::assayNames(pe)[1],
  k = 15L,
  returnBias = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{pe}{a PBMExperiment object containing PBM intensity data.}

\item{assay}{a string name of the assay to adjust. (default = \code{SummarizedExperiment::assayNames(pe)[1]})}

\item{k}{an integer specifying the size of the region to use to for computing
local bias. Must be odd. (default = 15L)}

\item{returnBias}{a logical whether to include the spatial bias as an
additional 'assay' (called 'spatialbias') in the returned
PBMExperiment object. (default = TRUE)}

\item{verbose}{a logical value whether to print verbose output
during analysis. (default = FALSE)}
}
\value{
Original PBMExperiment object with assay containing spatially adjusted intensities
(\code{"normalized"}) and a new column added to the colData,
\code{"spatialMedian"}, containing the global median intensity of the original
probe-level data used to compute spatial bias for each sample.
If specified, the estimated spatial bias will also be included in an additional
assay (\code{"spatialbias"}). If assays with the same name are already included
in the object, they will be overwritten.
}
\description{
Given a PBMExperiment containing probe intensity data with array coordinates
specified in rowData as \code{"Column"} and \code{"Row"}, this function
computes the local intensity bias for each probe, and returns the
same PBMExperiment object with an additional array containing spatially
corrected (scaled) probe intensities. Optionally, the per-probe bias is also
returned as an array in the bias corrected PBMExperiment object.

The spatial bias at each probe is computed as the ratio between the median
intensity in a \code{k} by \code{k}
region surrounding the probe and the median intensity of all probes
across the array. This approach is taken directly from the Universal PBM
Analysis Suite described in Berger and Bulyk (Nature Protocols, 2008).

Only probes passing all \code{probeFilter} rules defined in the original PBMExperiment
object are spatially adjusted and included in the returned object.
}
\references{
\itemize{
\item Berger, M. F., & Bulyk, M. L. (2009). Universal protein-binding microarrays for the comprehensive characterization of the DNA-binding specificities of transcription factors. Nature Protocols, 4(3), 393-411.
}
}
\author{
Patrick Kimes
}
