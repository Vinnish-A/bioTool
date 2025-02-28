
# tidy --------------------------------------------------------------------


#' Calculate Row Sums and Add as a New Column
#'
#' This function calculates the row sums of selected numeric columns in a data frame
#' and adds the result as a new column with a specified name.
#'
#' @param df A data frame containing the data.
#' @param name A string specifying the name of the new column where the row sums will be stored.
#' @param con An expr about how to select columns.
#'            Defaults to `where(is.numeric)`, which selects all numeric columns.
#'
#' @return A modified data frame with the new column containing row sums.
#'
#' @import dplyr
#' @import rlang
#'
#' @examples
#' library(dplyr)
#' df = tibble(a = 1:3, b = 4:6, c = letters[1:3])
#' df |> tidyRowSums("row_sum")
#'
#' @export
tidy_rowSums = function(df, name, con = where(is.numeric)) {

  con = enquo(con)

  df[[name]] = rowSums(df |> select(!!con))

  return(df_)

}

#' tidy_addCol
#'
#' @noRd
#'
#' @export
tidy_addCol = function(df, name, fun, con = where(is.numeric)) {

  con = enquo(con)

  df[[name]] = fun(df |> select(!!con))

  return(df_)

}

#' Transpose a tibble
#'
#' This function transposes a data frame, converting its columns into rows,
#' and optionally renames the rownames column in the output tibble.
#'
#' @param df A data frame to be transposed.
#' @param id_new A string specifying the name of the new column that will contain the original column names (now row names).
#' @param id_old A string specifying the name of the column to be used as rownames before transposing.
#'               Defaults to the first column of the input data frame.
#'
#' @return A tibble with the transposed data, including a column with the original column names.
#'
#' @import tibble
#'
#' @examples
#' df = tibble(id = c("row1", "row2"), a = c(1, 2), b = c(3, 4))
#' df |> tidyT("new_id")
#'
#' @export
tidy_t = function(df, id_new, id_old = colnames(df_)[[1]]) {

  df |>
    column_to_rownames(id_old) |>
    t() |>
    as_tibble(rownames = id_new)

}
