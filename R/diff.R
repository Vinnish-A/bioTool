
# diff --------------------------------------------------------------------


## check_input ----

checkGroup = function(group_) {

  assert_that(is.factor(group_))
  assert_that(length(levels(group_)) == 2)

}

## bulk ----

#' Convert Raw Counts to TPM (Transcripts Per Million)
#'
#' This function converts raw gene count matrix into TPM values, using gene length information and effective length from a reference dataset.
#' The reference dataset is organism-specific (human or mouse), and the function uses the specified keytype (e.g., Ensembl, Entrez, or Symbol) to map gene IDs.
#'
#' @param mat_ A numeric matrix of raw gene counts, where rows are genes and columns are samples.
#' @param org_ A character string indicating the organism, either 'hsa' for human or 'mmu' for mouse (default is 'hsa').
#' @param keytype_ A character string specifying the keytype to map gene IDs. Options are 'ensembl', 'entrez', or 'symbol' (default is 'symbol').
#' @return A numeric matrix of TPM values (genes in rows, samples in columns).
#'
#' @examples
#' mat_count = read.csv(fileWhenTest('bulk/count.csv', 'bioTool'), row.names = 'symbol')
#' mat_tpm = count2tpm(mat_count)
#'
#' @export
count2tpm = function(mat_, org_ = 'hsa', keytype_ = 'symbol') {

  org_ = tolower(org_)
  keytype_ = tolower(keytype_)

  match.arg(org_, c('hsa', 'mmu'))
  match.arg(keytype_, c('ensembl', 'entrez', 'symbol'))

  ref_ = readRDS(fileWhenTest(file.path('reference/', ifelse(org_ == 'hsa', 'range_hg38.rds', 'range_vm32.rds')), 'bioTool'))
  ref_ = tibble(id = ref_[[keytype_]], len = ref_[['eff_length']]) |>
    filter(id %in% rownames(mat_)) |>
    filter(!is.na(len)) |>
    distinct(id, .keep_all = T) |>
    arrange(id)

  count2tpmIn(mat_[ref_$id, ], ref_$len)

}

#' Helper Function to Convert Raw Counts to TPM
#'
#' This helper function performs the actual conversion of raw counts to TPM (Transcripts Per Million) values by
#' normalizing the counts based on gene length and scaling across samples.
#'
#' @param mat_ A numeric matrix of raw gene counts, where rows are genes and columns are samples.
#' @param len_ A numeric vector containing the effective length of each gene (in kilobases) corresponding to the rows of the count matrix.
#' @return A numeric matrix of TPM values (genes in rows, samples in columns).
#'
#' @noRd
#'
#' @keywords internal
count2tpmIn = function(mat_, len_) {

  mat_ = mat_ / c(len_/1000)
  mat_ = 1e6 * t(t(mat_) / colSums(mat_))

  return(mat_)

}

#' Log-Normalize Gene Expression Data
#'
#' This function performs log-normalization on a gene expression matrix, adding a constant value (default is 1) to avoid taking the logarithm of zero.
#' The normalization is done element-wise, applying the specified logarithm function (default is `log2`).
#'
#' @param mat_ A numeric matrix of gene expression values (genes in rows, samples in columns).
#' @param log_ A function for the logarithmic transformation (default is `log2`).
#' @param sf_ A scaling factor (default is 1), though it is not used in the current implementation.
#' @return A numeric matrix of log-transformed gene expression values.
#'
#' @examples
#' mat_count = read.csv(fileWhenTest('bulk/count.csv', 'bioTool'), row.names = 'symbol')
#' mat_tpm = count2tpm(mat_count)
#' mat_tpm = logNormalize(mat_tpm)
#'
#' @export
logNormalize = function(mat_, log_ = log2, sf_ = 1) {

  if (is.matrix(mat_)) {

    log_(mat_ + sf_)

  } else {

    mat_ |>
      auto_apply(1, \(vec__) log_(vec__ + 1)) |>
      t()

  }

}

merge_expression = function(df, by, symbol = 'symbol') {

  match.arg(by, c('max', 'min', 'mean'))

  symbol_duplicated = df[[symbol]][duplicated(df[[symbol]])]

  if (length(symbol_duplicated) != 0) {

    df_remain = df[!df[[symbol]] %in% symbol_duplicated, ]
    df_patch = df[df[[symbol]] %in% symbol_duplicated, ]
    lst_patch = split(df_patch[, -1], df_patch[[symbol]])

    if (by == 'max') {

      df_patched = lst_patch |>
        map(~ .x |> auto_lapply(max)) |>
        bind_rows()

    } else if (by == 'min') {

      df_patched = lst_patch |>
        map(~ .x |> auto_lapply(min)) |>
        bind_rows()

    } else if (by == 'mean') {

      df_patched = lst_patch |>
        map(~ .x |> auto_lapply(mean)) |>
        bind_rows()

    }

    df_patched = df_patched |>
      mutate(symbol = names(lst_patch)) |>
      relocate(symbol)

    colnames(df_patched) = colnames(df_remain)

    return(bind_rows(df_remain, df_patched))

  } else {

    return(df)

  }

}

## DA ----

#'
#' This function performs differential expression analysis using the Limma package.
#'
#' @param normal_matrix A numeric matrix of normalized expression values, where rows are genes and columns are samples.
#' @param group A factor vector indicating the group assignments for the samples.
#'
#' @return A tibble containing the differential expression results with columns:
#'   - `symbol`: Gene identifiers.
#'   - `logFC`: Log2 fold change.
#'   - `AveExpr`: Average expression level.
#'   - `t`: Moderated t-statistic.
#'   - `B`: Beta-statistic.
#'   - `p.value`: Raw p-value.
#'   - `padj`: Adjusted p-value (FDR).
#'
#' @export
diff_limma = function(normal_matrix, phen, design) {

  design = model.matrix(design, data = phen)
  fit = limma::lmFit(normal_matrix, design)
  fit = limma::eBayes(fit)
  result = limma::topTable(fit, coef = nth(colnames(fit$coefficients), -1), adjust.method = 'fdr', number = Inf)
  result = result |>
    as_tibble(rownames = 'symbol') |>
    rename(pvalue = P.Value, padj = adj.P.Val) |>
    relocate(B, .after = t)

  return(result)

}

#' Differential Expression Analysis using DESeq2
#'
#' This function performs differential expression analysis using the DESeq2 package.
#'
#' @param count_matrix A numeric matrix of raw count data, where rows are genes and columns are samples.
#' @param group A factor vector indicating the group assignments for the samples.
#'
#' @return A tibble containing the differential expression results with columns:
#'   - `symbol`: Gene identifiers.
#'   - `logFC`: Log2 fold change.
#'   - `lfcSE`: Standard error of the log2 fold change.
#'   - `stat`: Wald statistic.
#'   - `pvalue`: Raw p-value.
#'   - `padj`: Adjusted p-value (FDR).
#'
#' @export
diff_deseq2 = function(count_matrix, phen, design) {

  dds = DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = phen,
    design = design
  )

  dds = DESeq2::DESeq(dds)
  deseqPairs = DESeq2::resultsNames(dds)
  res = DESeq2::results(dds, name = deseqPairs[[length(deseqPairs)]])

  res_df = as.data.frame(res)
  res_df = res_df[order(res_df$padj), ]

  result = res_df |>
    as_tibble(rownames = 'symbol') |>
    dplyr::rename(logFC = log2FoldChange)

  return(result)

}


#' Differential Expression Analysis using EdgeR
#'
#' This function performs differential expression analysis using the EdgeR package.
#'
#' @param count_matrix A numeric matrix of raw count data, where rows are genes and columns are samples.
#' @param group A factor vector indicating the group assignments for the samples.
#'
#' @return A tibble containing the differential expression results with columns:
#'   - `symbol`: Gene identifiers.
#'   - `logFC`: Log2 fold change.
#'   - `logCPM`: Log counts per million.
#'   - `LR`: Likelihood ratio statistic.
#'   - `pvalue`: Raw p-value.
#'   - `padj`: Adjusted p-value (FDR).
#'
#' @export
diff_edger = function(count_matrix, phen, design) {

  y = edgeR::DGEList(counts = count_matrix)

  keep = edgeR::filterByExpr(y)
  y = y[keep, , keep.lib.sizes = FALSE]

  y = edgeR::calcNormFactors(y)
  design = model.matrix(design, data = phen)
  y = edgeR::estimateDisp(y, design)
  fit = edgeR::glmFit(y, design)
  lrt = edgeR::glmLRT(fit, coef = colnames(design)[[ncol(design)]])
  res = edgeR::topTags(lrt, n = Inf)

  res_df = as.data.frame(res)
  result = res_df |>
    as_tibble(rownames = 'symbol') |>
    rename(pvalue = PValue, padj = FDR)

  return(result)

}

diff_series = function(series, design, algorithm) {

  algorithm = tolower(algorithm)
  match.arg(algorithm, c('limma', 'deseq2', 'edger'))

  fun = switch(
    algorithm,
    limma = diff_limma,
    deseq2 = diff_deseq2,
    edger = diff_edger
  )

  fun(series$expr |> column_to_rownames('sample') |> t(), series$phen, design)

}


## test ----

test_diff = lst(
  limma = substitute({
    mat_count = read.csv(fileWhenTest('bulk/count.csv', 'bioTool'), row.names = 'symbol')
    mat_tpm = count2tpm(mat_count) |> logNormalize()
    group = factor(substr(colnames(mat_tpm), 1, 1), levels = c('C', 'T'))
    deg_limma = diff_limma(mat_tpm, group)
  })
)

