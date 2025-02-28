
# metric ------------------------------------------------------------------


## regr ----

#' Metrics Regression Function
#'
#' A function to retrieve the list of regression metrics stored in the environment.
#'
#' @return A tibble containing the metadata of regression metrics.
#'
#' @export
metrics_regr = function() {

  checkReliance('mlr3')

  bind_rows(as.list(env_meta_regr))

}

#' Transform Target and Feature into Prediction Table for Regression
#'
#' This function takes a target (true values) and a feature (predicted values) and
#' converts them into a `data.table` object to be used for regression scoring.
#'
#' @param target_ A vector of true values.
#' @param feature_ A vector of predicted values.
#'
#' @return A `PredictionRegr` object containing the target and feature.
#'
#' @keywords internal
tab2pre_regr = function(target_, feature_) {

  data.frame(row_ids = seq_along(target_), truth = target_, response = feature_) |>
    mlr3::as.data.table() |>
    mlr3::as_prediction_regr()

}

#' Produce Regression Metric Function
#'
#' This function produces a metric evaluation function based on the supplied
#' metric name, which can then be used to evaluate regression models.
#'
#' @param metric_ A string containing the name of the metric to be used for evaluation in mlr3.
#'
#' @return A function that computes the score for a given target and feature based on the metric.
#'
#' @keywords internal
produce_metric_regr = function(metric_) {

  force(metric_)

  function(target, feature) {

    tab2pre_regr(target, feature)$score(msr(metric_))

  }

}

#' Add Custom Regression Metric
#'
#' This function allows you to add a custom regression metric to the environment,
#' which can be used with other functions that expect regression metrics.
#'
#' @param fun_ A function representing the custom regression metric.
#' @param key_ A string representing the key for the custom metric.
#' @param description_ A description of the custom metric (optional).
#'
#' @return NULL (invisible).
#'
#' @export
add_custom_regr = function(fun_, key_, description_ = '') {

  checkReliance('mlr3')

  meta_ = tibble(key = paste0('regr.', key_), label = description_, type = 'regr', backend = 'custom')

  assign('meta_regr_custom', bind_rows(env_meta_regr$meta_regr_custom, meta_), envir = env_meta_regr)
  assign(key_, fun_, envir = metric_regr)

  invisible()

}

### meta ----

#' Initialize Regression Meta Information
#'
#' This function initializes the metadata for regression metrics, checking whether
#' the `mlr3` package is installed. If `mlr3` is available, it loads the regression
#' metrics from `mlr3` and formats them. If not, it returns `NULL`.
#'
#' @keywords internal
init_meta_regr = function() {

  if ('mlr3' %in% rownames(installed.packages())) {

    res_ = list2env(
      lst(
        meta_regr_builtin = mlr3::msrs() |>
          mlr3::as.data.table() |>
          as_tibble() |>
          select(where(~ !is.list(.x))) |>
          filter(task_type == 'regr') |>
          select(-predict_type) |>
          rename(type = task_type) |>
          mutate(backend = 'mlr3measures'),
        meta_regr_custom = meta_regr_builtin[0, ]
      )
    )

  } else {

    res_ = NULL

  }

  return(res_)

}

#' Environment for Regression Meta Information
#'
#' This object stores the environment for regression metric metadata, initialized by
#' `init_meta_regr()`. It contains both the built-in and custom regression metrics.
#'
#' @return An environment containing `meta_regr_builtin` and `meta_regr_custom`.
#'
#' @export
env_meta_regr = init_meta_regr()

### object ----

#' Environment for Regression Meta Information
#'
#' This object stores the environment for regression metric metadata, initialized by
#' `init_meta_regr()`. It contains both the built-in and custom regression metrics.
#'
#' @keywords internal
init_metric_regr = function() {

  res_ = lapply(env_meta_regr$meta_regr_builtin$key, produce_metric_regr) |>
    set_names(sapply(strsplit(env_meta_regr$meta_regr_builtin$key, '\\.'), `[[`, 2)) |>
    list2env()

  if ('mlr3' %in% rownames(installed.packages())) {

    res_ = lapply(env_meta_regr$meta_regr_builtin$key, produce_metric_regr) |>
      set_names(sapply(strsplit(env_meta_regr$meta_regr_builtin$key, '\\.'), `[[`, 2)) |>
      list2env()

  } else {

    res_ = NULL

  }

  return(res_)

}

#' Environment for Regression Metric Functions
#'
#' This object stores the environment for regression metric functions, initialized by
#' `init_metric_regr()`. It contains the functions to compute scores for built-in regression
#' metrics from `mlr3`.
#'
#' @export
metric_regr = init_metric_regr()

### custom ----

#' Pearson Correlation Coefficient (PCC)
#'
#' This custom function calculates the Pearson Correlation Coefficient (PCC) between
#' two vectors: `target` (true values) and `feature` (predicted values). It is commonly
#' used to evaluate the linear relationship between predicted and true values in regression.
#'
#' @keywords internal
pcc = function(target, feature) {

  cor(target, feature)

}

add_custom_regr(pcc, 'prho', 'Pearson Corr Coef')
