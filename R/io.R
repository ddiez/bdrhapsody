#' read_bd
#'
#' Read file in BD format.
#'
#' @param filename name of file to read.
#' @param samplename sample names to attach to cell_index.
#'
#' @export
read_bd <- function(...) {
  UseMethod("read_bd")
}

#' @rdname read_bd
#'
#' @export
read_bd.data.frame <- function(x, split = FALSE) {
  if (split) {
    y <- vector("list", nrow(x))
    for (k in seq_len(nrow(x))) {
      ytmp <- read_bd(x[["filename"]][k], samplename = x[["samplename"]][k])

      cdata <- colData(ytmp) %>% as.data.frame()
      cdata <- left_join(cdata, x, by = c("samplename", "filename"))
      rownames(cdata) <- cdata[["cell_index"]]
      colData(ytmp) <- DataFrame(cdata)

      # set rowData.
      ann <- strsplit(rownames(ytmp), "\\|")
      ann <- tibble(
        id = rownames(ytmp),
        symbol = sapply(ann, `[`, 1),
        accession = sapply(ann, `[`, 2),
        note = sapply(ann, `[`, 3)
      )

      rdata <- DataFrame(ann, row.names = rownames(ytmp))
      rowData(ytmp) <- rdata

      y[[k]] <- ytmp
    }

  } else {
    y <- read_bd(x[["filename"]], x[["samplename"]])

    # set colData.
    cdata <- colData(y) %>% as.data.frame()
    cdata <- left_join(cdata, x, by = c("samplename", "filename"))
    rownames(cdata) <- cdata[["cell_index"]]
    colData(y) <- DataFrame(cdata)

    # set rowData.
    ann <- strsplit(rownames(y), "\\|")
    ann <- tibble(
      id = rownames(y),
      symbol = sapply(ann, `[`, 1),
      accession = sapply(ann, `[`, 2),
      note = sapply(ann, `[`, 3)
    )

    rdata <- DataFrame(ann, row.names = rownames(y))
    rowData(y) <- rdata
  }

  y
}

#' @rdname read_bd
#'
#' @export
read_bd.character <- function(filename, samplename = NULL) {
  require(SingleCellExperiment)

  if (is.null(samplename))
    samplename <- paste0("S", sprintf("%05d", length(filename)))

  l <- mapply(read_bd_as_df, filename = filename, samplename = samplename, SIMPLIFY = FALSE)
  x <- bind_rows(l, .id = "filename")

  y <- to_matrix(x %>% select(-samplename, -filename), "cell_index") %>% t()

  cdata <- x %>% select(samplename, filename, cell_index)

  SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = y,
      logcounts = log10(y + 1)),
    colData = cdata
  )
}

#' @rdname read_bd
read_bd_as_df <- function(filename, samplename = NULL) {
  x <- readr::read_csv(
    filename, comment = "#",
    progress = FALSE,
    col_types = cols(.default = col_integer())
  ) %>%
    dplyr::rename(cell_index = Cell_Index)

  if (!is.null(samplename))
    x <- x %>% mutate(samplename = !!samplename, cell_index = paste0(!!samplename, "_", cell_index))

  x
}
