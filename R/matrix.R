#' to_matrix
#'
#' Convert a data.frame into a matrix, optionally setting one of the columns
#' as the row names.
#'
#' @param x data.frame to conver.
#' @param rownames column to be used as row names.
to_matrix <- function(x, rownames = NULL) {
  UseMethod("to_matrix")
}

#' @rdname to_matrix
to_matrix.data.frame <- function(x, rownames = NULL) {
  if (!is.null(rownames)) {
    rn <- x %>% pull(!!rownames)
    x <- x %>% select(- !!rownames)
  }

  y <- x %>% as.matrix()

  if (!is.null(rownames))
    rownames(y) <- rn

  y
}
