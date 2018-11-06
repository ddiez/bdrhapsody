#' read_bd
#'
#' Read file in BD format.
#'
#' @param x a data.frame with sample information.
#' @param split logical; whether to return a list of objects, each for a file.
#' @param filename name of file to read.
#' @param samplename sample names to attach to cell_index.
#' @param org.db and OrgDb object or a character string of such object.
#' @param ... arguments passed down to methods.
#'
#' @export
read_bd <- function(...) {
  UseMethod("read_bd")
}

#' @rdname read_bd
#'
#' @export
read_bd.data.frame <- function(x, split = FALSE, ...) {
  if (split) {
    y <- vector("list", nrow(x))
    for (k in seq_len(nrow(x))) {
      ytmp <- read_bd(x[["filename"]][k], samplename = x[["samplename"]][k])

      # add additional colData information.
      cdata <- colData(ytmp) %>% as.data.frame()
      cdata <- left_join(cdata, x, by = c("samplename", "filename"))
      rownames(cdata) <- cdata[["cell_index"]]
      colData(ytmp) <- DataFrame(cdata)

      y[[k]] <- ytmp
    }

  } else {
    y <- read_bd(x[["filename"]], x[["samplename"]])

    # add additional colData information.
    cdata <- colData(y) %>% as.data.frame()
    cdata <- left_join(cdata, x, by = c("samplename", "filename"))
    rownames(cdata) <- cdata[["cell_index"]]
    colData(y) <- DataFrame(cdata)
  }

  y
}

#' @rdname read_bd
#'
#' @export
read_bd.character <- function(filename, samplename = NULL, org.db = NULL, ...) {
  if (is.null(samplename))
    samplename <- paste0("S", sprintf("%05d", length(filename)))

  l <- mapply(read_bd_as_df, filename = filename, samplename = samplename, SIMPLIFY = FALSE)
  x <- bind_rows(l, .id = "filename")

  y <- to_matrix(x %>% select(-samplename, -filename), "cell_index") %>% t()

  # colData.
  cdata <- x %>% select_("samplename", "filename", "cell_index")

  # set rowData.
  rdata <- DataFrame(parse_bd_id(rownames(y)))

  if (! is.null(org.db)) {
    if (class(org.db) != "OrgDb")
      warning("org.db object not of class OrgDb (class: ", class(org.db), ")")
    else
      rdata[["entrezgene"]] <- mapIds(org.db, rdata[["symbol"]], column = "ENTREZID", keytype = "SYMBOL")
  }

  SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = y,
      logcounts = log10(y + 1)),
    colData = cdata,
    rowData = rdata
  )
}

#' @rdname read_bd
read_bd_as_df <- function(filename, samplename = NULL) {
  x <- readr::read_csv(
    filename, comment = "#",
    progress = FALSE,
    col_types = cols(.default = col_integer())
  ) %>%
    rename(cell_index = .data$Cell_Index)

  if (!is.null(samplename))
    x <- x %>% mutate(samplename = !!samplename, cell_index = paste0(!!samplename, "_", .data$cell_index))

  x
}


#' parse_bd_id
#'
#' @param x character vector with BD format ids.
#'
#' @export
parse_bd_id <- function(x) {
  fields <- strsplit(x, "\\|")
  data.frame(
    id = x,
    symbol = sapply(fields, `[`, 1),
    accession = sapply(fields, `[`, 2),
    note = sapply(fields, `[`, 3),
    row.names = x,
    stringsAsFactors = FALSE
  )
}
