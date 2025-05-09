% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R
\name{diff_deseq2}
\alias{diff_deseq2}
\title{Differential Expression Analysis using DESeq2}
\usage{
diff_deseq2(count_matrix, phen, design)
}
\arguments{
\item{count_matrix}{A numeric matrix of raw count data, where rows are genes and columns are samples.}

\item{group}{A factor vector indicating the group assignments for the samples.}
}
\value{
A tibble containing the differential expression results with columns:
  - `symbol`: Gene identifiers.
  - `logFC`: Log2 fold change.
  - `lfcSE`: Standard error of the log2 fold change.
  - `stat`: Wald statistic.
  - `pvalue`: Raw p-value.
  - `padj`: Adjusted p-value (FDR).
}
\description{
This function performs differential expression analysis using the DESeq2 package.
}
