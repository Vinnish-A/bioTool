
# GEO ---------------------------------------------------------------------

## Anno ----

checkGPL = function(GPL = NULL) {

  if (length(GPL) == 0) {
    stop("please input GPL number")
  }
  GPLList = getGPLList()
  flag = (GPL %in% GPLList[, 1])

  return(flag)

}

getGPLList = function() {
  tryCatch(utils::data("gpl_list", package = "AnnoProbe"))
  GPLList = get("gpl_list")
  return(GPLList[, 1:2])
}

idmap = function(gpl = "GPL570", type = "bioc", mirror = "tercent", destdir = 'tmp') {

  gpl = toupper(gpl)
  gpl_anno = paste(gpl, c("bioc", "soft", "pipe"), sep = "_")
  if (mirror == "tercent") {
    up = "http://49.235.27.111"
  }
  if (!checkGPL(gpl)) {
    stop(
      "This platform is not in our list, please use our shinyAPP to custom annotate your probe sequences, or ask us to process and then update the R package!"
    )
  }
  else {
    tryCatch(
      utils::data("exists_anno_list", package = "AnnoProbe")
    )
    gpl_anno = gpl_anno[gpl_anno %in% exists_anno_list]
    if (T) {
      tpf = paste0(paste(gpl, type, sep = "_"), ".rda")
      down = paste0("/GEOmirror/GPL/", tpf)
      download.file(paste0(up, down), file.path(destdir, tpf), mode = "wb")
      load(file.path(destdir, tpf))
      return(get(paste(gpl, type, sep = "_")))
    }
    else {
      stop("We have that platform, but just offer other type of annotaion.")
    }
  }

}

## series ----

checkSeries = function(series) {

  name = as.character(substitute(series))
  series = get(name, envir = parent.frame())

  if (!identical(series$expr$sample, series$phen$sample)) {

    series = asSeries(series$expr, series$phen)

    assign(name, series, envir = parent.frame())

  }

  invisible()

}

as_series = function(expr, phen) {

  expr = as_tibble(expr)

  phen = expr[, 1] |>
    left_join(phen) |>
    as_tibble()

  series = list(expr = expr, phen = phen)

  class(series) = c('series', class(series))

  return(series)

}

read1series = function(series_, annotate_, table_id_) {

  gse_ = GEOquery::getGEO(filename = series_, getGPL = F)

  expr_ = gse_ |> Biobase::exprs() |> as_tibble(rownames = 'ID')
  phen_ = gse_ |> Biobase::pData() |> as_tibble() |> mutate(sample = geo_accession)

  if (annotate_) {

    if (checkGPL(gse_@annotation)) {

      table_id_ = idmap(gse_@annotation, type = 'soft', destdir = 'tmp') |>
        pull(symbol, ID)

    } else if (!is.null(table_id_)) {

      table_id_ = table_id_

    } else {

      warning('GPL Not Found With No Defaults')
      ids_ = unique(expr_$ID)
      table_id_ = setNames(ids_, ids_)

    }

  } else {

    ids_ = unique(expr_$ID)
    table_id_ = setNames(ids_, ids_)

  }

  expr_tidy_ = expr_ |>
    mutate(symbol = table_id_[ID]) |>
    filter(symbol != 'permuted_negative') |>
    distinct(symbol, .keep_all = T) |>
    column_to_rownames('symbol') |>
    select(-ID) |>
    t() |>
    as_tibble(rownames = 'sample')

  if (all(!str_detect(expr_tidy_$sample, 'GSM')) & nrow(expr_tidy_) == nrow(phen_)) {
    expr_tidy_$sample = phen_$sample
  }

  phen_tidy_ = expr_tidy_[, 1] |>
    left_join(phen_) |>
    select(sample, title, contains(':ch')) |>
    rename_with(~ .x |> str_sub(end = -5), contains(':ch')) |>
    mutate(platform = gse_@annotation)

  return(as_series(expr = expr_tidy_, phen = phen_tidy_))

}

#' Read and Process a GEO Series Matrix File
#'
#' This function reads a GEO Series Matrix file from a specified directory, processes the expression
#' and phenotype data, and returns them in tidy formats.
#'
#' @param dir_ A string specifying the directory containing the GEO Series Matrix file. The directory must contain exactly one file.
#' @param alwaysGZ_ A boolean value indicating whether to always decompress gz-compressed txt files in the folder.
#'
#' @return A list with two components:
#'   - `expr`: A tidy tibble of expression data, where rows are samples, and columns are gene symbols.
#'   - `phen`: A tidy tibble of phenotype data, including metadata for each sample.
#'
#' @details
#' - The function checks for and processes the platform annotation (GPL). If the annotation cannot be found,
#'   a warning is issued, and the gene IDs are retained as they appear in the expression data.
#' - Genes with the identifier `'permuted_negative'` are removed.
#'
#' @export
readSeries = function(dir_, annotate_ = T, table_id_ = NULL) {

  series_ = list.files(dir_, full.names = T)

  if (length(series_) == 1) {

    return(read1series(series_, annotate_, table_id_))

  } else if (length(series_) > 1) {

    lst_res_ = lapply(series_, read1series, annotate_ = annotate_, table_id_ = table_id_)

    return(
      as_series(
        expr = bind_rows(lapply(lst_res_, `[[`, 'expr')),
        phen = bind_rows(lapply(lst_res_, `[[`, 'phen'))
      )
    )

  }

}

#' Write GEO Series Data to Disk
#'
#' This function writes processed GEO Series data (expression and phenotype) to specified files in a given directory.
#'
#' @param series_ A list containing two components:
#'   - `expr`: A tibble of expression data.
#'   - `phen`: A tibble of phenotype data.
#' @param dir_ A string specifying the directory where the files will be written. If the directory does not exist, it will be created.
#' @param RData A logical value. If `TRUE`, the function saves the data as an RData file in addition to CSV files. Defaults to `TRUE`.
#' @param gse_ (Optional) A string specifying a GEO Series accession ID (e.g., "GSE12345"). If provided, it is appended to the output filenames.
#'
#' @return None. The function writes files to the specified directory.
#'
#' @details
#' - The function creates the output directory if it does not exist.
#' - Outputs include:
#'   - `expr.csv`: Expression data in CSV format.
#'   - `phen.csv`: Phenotype data in CSV format.
#'   - `all.RData`: An RData file containing the entire dataset (optional).
#'
#' @export
writeSeries = function(series_, dir_, RData = T, gse_ = NULL) {

  if (!dir.exists(dir_)) dir.create(dir_, recursive = T)

  if (!is.null(gse_)) names(series_) = paste(names(series_), gse_, sep = '_')

  cat('Writing expr...\n')
  write_csv(series_$expr, normalizePath(file.path(dir_, 'expr.csv'), '/', F))
  cat('Writing phen...\n')
  write_csv(series_$phen, normalizePath(file.path(dir_, 'phen.csv'), '/', F))

  if (RData) {
    cat('Writing RData...\n')
    save(series_, file = normalizePath(file.path(dir_, 'all.RData'), '/', F))
  }

}

## modify ----

#' filter.series
#'
#' @noRd
#'
#' @exportS3Method dplyr::filter
filter.series = function(.data, ..., .by = NULL, .preserve = FALSE) {

  .data$phen = .data$phen |>
    filter(...)

  toKeep = .data$phen$sample

  .data$expr = .data$expr |>
    filter(sample %in% toKeep)

  checkSeries(.data)

  return(.data)

}

#' mutate.series
#'
#' @noRd
#'
#' @exportS3Method dplyr::mutate
mutate.series = function(.data, ...) {

  .data$phen = .data$phen |>
    mutate(...)

  checkSeries(.data)

  return(.data)

}

#' select.series
#'
#' @noRd
#'
#' @exportS3Method dplyr::select
select.series = function(.data, ...) {

  left_join(.data$expr, .data$phen) |>
    select(...)

}

#' rename.series
#'
#' @noRd
#'
#' @exportS3Method dplyr::rename
rename.series = function(.data, ...) {

  .data$phen = .data$phen |>
    rename(...)

  checkSeries(.data)

  return(.data)

}
