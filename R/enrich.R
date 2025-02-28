
# enrich ------------------------------------------------------------------


## dataset ----

#' Retrieve GO Data
#'
#' This function retrieves Gene Ontology (GO) data for a specified organism.
#'
#' @param organism_ A string specifying the organism. `'hsa'` for human and `'mmu'` for mouse. Defaults to `'hsa'`.
#' @param keyType_ A string specifying the key type for gene identifiers. Defaults to `'SYMBOL'`.
#'
#' @return A list of GO data, including mappings of terms to genes and term names.
#'
#' @examples
#' GO = getGO()
#'
#' @export
getGO = function(organism_ = 'hsa', keyType_ = 'SYMBOL') {

  checkReliance('clusterProfiler')

  if (organism_ == 'hsa') {
    org_ = 'org.Hs.eg.db'
  } else {
    org_ = 'org.Mm.eg.db'
  }

  GO_ = clusterProfiler:::get_GO_data(org_, "ALL", keyType_)
  GO_$NAME2PATHID = setNames(names(GO_$PATHID2NAME), GO_$PATHID2NAME)

  res_ = mget(ls(envir = GO_), envir = GO_)

  return(res_)

}

#' Retrieve Genes Associated with a GO Term
#'
#' This function retrieves the genes associated with a specified GO term.
#'
#' @param name_ A string specifying the name of the GO term.
#' @param GO_ A list of GO data, typically obtained from `getGO()`.
#'
#' @return A vector of gene identifiers associated with the GO term.
#'
#' @examples
#' GO = getGO()
#' genes = symbolIn("cell cycle", GO_ = GO)
#'
#' @export
symbolIn = function(name_, GO_) {

  res_ = GO_$PATHID2EXTID[[GO_$NAME2PATHID[name_]]]

  if (is.null(res_)) stop('term ', name_, ' not found')

  return(res_)

}

#' Convert a List to GMT Format
#'
#' This function converts a list of terms and associated genes into a GMT-like format.
#'
#' @param lst_ A named list where names are terms and values are vectors of associated genes.
#'
#' @return A tibble with two columns: `term` and `gene`.
#'
#' @examples
#' lst = list(Term1 = c("GeneA", "GeneB"), Term2 = c("GeneC"))
#' gmt = lst2gmt(lst)
#'
#' @export
lst2gmt = function(lst_) {

  tibble(
    term = rep(names(lst_), sapply(lst_, length)),
    gene = unlist(lst_)
  )

}


## easy_enrich ----

#' Perform Easy Enrichment Analysis
#'
#' This function performs GO or KEGG enrichment analysis for a given set of genes.
#'
#' @param genes_ A vector of gene symbols to analyze.
#' @param toGo_ A string specifying the type of analysis: `'GO'` or `'KEGG'`. Defaults to `'GO'`.
#' @param cutoff_ A numeric value specifying the FDR cutoff. Defaults to `0.05`.
#' @param organism_ A string specifying the organism. `'hsa'` for human and `'mmu'` for mouse. Defaults to `'hsa'`.
#'
#' @return An enrichment analysis result object.
#'
#' @import tidyr
#'
#' @examples
#' enrichment = enrichment_of(c("TP53", "BRCA1"))
#'
#' @export
enrichment_of = function(genes_, toGo_ = "GO", cutoff_ = 0.05, organism_ = 'hsa') {

  checkReliance('clusterProfiler')

  if (organism_ == 'hsa') {
    org_ = 'org.Hs.eg.db'
  } else {
    org_ = 'org.Mm.eg.db'
  }

  if (toGo_ == "KEGG") {
    kegg_category_ = clusterProfiler:::kegg_category_data()

    table_kegg_ = kegg_category_ |>
      mutate(id = paste0(organism_, id)) |>
      pull(name, id)

    gid_ = clusterProfiler::bitr(genes_, 'SYMBOL', 'ENTREZID', OrgDb = org_) |> drop_na()
    genes_ = gid_$ENTREZID
    res_ = clusterProfiler::enrichKEGG(genes_, pAdjustMethod = "fdr", qvalueCutoff = cutoff_, organism = organism_, use_internal_data = TRUE)

    res_@result$Description = table_kegg_[res_@result$Description]

  } else if (toGo_ == "GO") {
    res_ = clusterProfiler::enrichGO(genes_, org_, pAdjustMethod = "fdr", qvalueCutoff = cutoff_, ont = "ALL", keyType = 'SYMBOL')
  }

  return(res_)

}

#' Enrichment Analysis for Activated and Suppressed Genes
#'
#' This function performs separate enrichment analyses for activated and suppressed genes in a dataset.
#'
#' @param df_ A data frame containing gene data with at least `logFC` and `symbol` columns.
#' @param toGo_ A string specifying the type of analysis: `'GO'` or `'KEGG'`. Defaults to `'GO'`.
#' @param cutoff_ A numeric value specifying the FDR cutoff. Defaults to `0.05`.
#' @param organism_ A string specifying the organism. `'hsa'` for human and `'mmu'` for mouse. Defaults to `'hsa'`.
#'
#' @return A list containing two enrichment analysis results: `act` (activated) and `sup` (suppressed).
#' @examples
#'
#' df = data.frame(symbol = c("TP53", "BRCA1", "MYC"), logFC = c(1.2, -0.5, 0.8))
#' enrich_actsup = enrichActSup(df, toGo_ = "GO", cutoff_ = 0.05)
#'
#' @export
enrichActSup = function(df_, toGo_ = "GO", cutoff_ = 0.05, organism_ = 'hsa') {
  list(
    act = df_ |>
      dplyr::filter(logFC > 0) |>
      pull(symbol) |>
      enrichment_of(toGo_ = toGo_, cutoff_ = cutoff_, organism_ = organism_),
    sup = df_ |>
      dplyr::filter(logFC < 0) |>
      pull(symbol) |>
      enrichment_of(toGo_ = toGo_, cutoff_ = cutoff_, organism_ = organism_)
  )
}

## simple_vis ----

#' Enhanced Dot Plot for ActSup-Enrichment Results
#'
#' This function generates a dot plot visualization of enrichment analysis results,
#' highlighting key terms and their significance.
#'
#' @param enrichLst_ A list containing enrichment results for activated (`act`) and suppressed (`sup`) genes.
#' @param nTerms_ An integer specifying the maximum number of terms to display. Defaults to `7`.
#'
#' @return A `ggplot` object showing the dot plot.
#' @examples
#' enrich_actsup = list(
#'   act = enrichment_of(c("TP53", "BRCA1")),
#'   sup = enrichment_of(c("MYC"))
#' )
#' dotPlotPlus(enrich_actsup, nTerms_ = 5)
#'
#' @export
dotPlotPlus = function(enrichLst_, nTerms_ = 7) {

  pToLabel_ = function(vec__) {
    vec__[is.na(vec__)] = 1
    cut(
      vec__,
      breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
      labels = c("FDR<0.0001", "FDR<0.001", "FDR<0.01", "FDR<0.05", " - ")
    )
  }

  splitTerms__ = function(Des__) {

    DesVec__ = str_split(Des__, " ", simplify = F)[[1]]
    cutPoint__ = ceiling(length(DesVec__)/2) + 1

    DesVec1__ = DesVec__[1:cutPoint__]; DesVec2__ = DesVec__[cutPoint__:length(DesVec__)]

    paste0(paste(DesVec1__, collapse = " "), "\n", paste(DesVec2__, collapse = " "))

  }

  enrichLst_$act@result$type = "Activated"
  enrichLst_$sup@result$type = "Suppressed"

  evalParse_ = \(x__) eval(parse(text = x__))

  dataPlot_ = enrichLst_ |>
    lapply(
      \(each__) each__@result |>
        filter(qvalue < 0.05) |>
        mutate(GeneRatio = map_dbl(GeneRatio, evalParse_)) |>
        slice_max(-qvalue, n = nTerms_) |>
        arrange(GeneRatio) |>
        mutate(Description = capitalize(Description),
               Description = map_chr(Description, ~ ifelse(str_count(.x, " ") > 3, splitTerms__(.x), .x)),
               Description = factor(Description, levels = Description),
               qvalue = pToLabel_(qvalue))
    ) |> bind_rows()

  dataPlot_ |>
    ggplot() +
    geom_point(aes(GeneRatio, Description, color = qvalue, size = Count)) +
    scale_color_manual(values = c("FDR<0.0001" = "#4e62ab", "FDR<0.001" = "#479db4", "FDR<0.01" = "#87d0a6", "FDR<0.05" = "#cbe99d", " - " = "#f7fbae")) +
    facet_wrap(~ type, nrow = 2, scales = "free_y") +
    theme_bw() +
    ylab("")

}



