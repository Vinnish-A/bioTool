% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrich.R
\name{dotPlotPlus}
\alias{dotPlotPlus}
\title{Enhanced Dot Plot for ActSup-Enrichment Results}
\usage{
dotPlotPlus(enrichLst_, nTerms_ = 7)
}
\arguments{
\item{enrichLst_}{A list containing enrichment results for activated (`act`) and suppressed (`sup`) genes.}

\item{nTerms_}{An integer specifying the maximum number of terms to display. Defaults to `7`.}
}
\value{
A `ggplot` object showing the dot plot.
}
\description{
This function generates a dot plot visualization of enrichment analysis results,
highlighting key terms and their significance.
}
\examples{
enrich_actsup = list(
  act = enrichment_of(c("TP53", "BRCA1")),
  sup = enrichment_of(c("MYC"))
)
dotPlotPlus(enrich_actsup, nTerms_ = 5)

}
