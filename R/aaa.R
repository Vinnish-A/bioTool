
# aaa ---------------------------------------------------------------------


checkReliance = function(...) {

  if (length(unlist(list(...))) == 0) stop('No input')

  pkgs_ = unlist(list(...))
  pkgs_ = pkgs_[!pkgs_ %in% rownames(installed.packages())]
  pkgs_ = paste(pkgs_, collapse = ', ')

  if (pkgs_ != '') {

    stop(pkgs_, ' are required but not installed. To enable this, please install. ')

  }

  invisible()

}
