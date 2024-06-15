#' Return Gene Wise Fit Statistics
#'
#' Convenient wrapper to interact with the fitted [NBAMSeqDataSet()] object to
#' extract the fit statistics.
#'
#' @param object An [NBAMSeqDataSet()] object that has been fitted with [NBAMSeq1()].
#' @param coef Name of the term. If NULL or mismatched, then a list of all terms will be returned.
#' @param number Number of rows to be returned after sorting. Defaults to 10.
#' @param genelist Return just the genes in this list. Expects a character vector.
#' @param adjust.method P.adj methods as detailed in [stats::p.adjust.methods]. Defaults to "BH".
#' @param sort.by Column to sort by. Defaults to "P.Value".
#' @param p.value Cut off for p-value; only rows with p-values smaller than this value are shown. Defaults to NULL.
#'
#' @return A data.frame with the statistics.
#' @export
topTable1 <- function(
    object,
    coef = NULL,
    number = 10,
    genelist = NULL,
    adjust.method = "BH",
    sort.by = "P.Value",
    p.value = NULL) {
  ## check whether NBAMSeq has been run before
  stopifnot(methods::is(object, "NBAMSeqDataSet"))
  if (!S4Vectors::metadata(object)$fitted) {
    stop("NBAMSeq function should be run before calling results function")
  }
  all_terms <- unlist(lapply(
    S4Vectors::mcols(object)$fit[[1]][c("p.table", "s.table")],
    \(x) {
      row.names(x)
    }
  ))

  # Pool terms into a single table
  if (is.null(coef) || any(!(coef %in% all_terms))) {
    stop(
      paste("Please select only one of the following coef:", paste(all_terms, collapse = ", "))
    )
  }
  stopifnot(length(coef) == 1 && is.character(coef))
  results <- as.data.frame(
    t(sapply(S4Vectors::mcols(object)$fit, function(x) {
      filter_terms(coef = coef, x)
    }))
  )

  # P Adjust
  stopifnot(length(adjust.method) == 1, is.character(adjust.method))
  names(results)[which(grepl("Pr\\(|p-value", names(results)))] <- "P.Value"
  if (!adjust.method %in% stats::p.adjust.methods) {
    stop(
      paste("Supported p.adjust:", paste(stats::p.adjust.methods, collapse = ", "))
    )
  }
  results$adj.P.Value <- stats::p.adjust(results$P.Value, method = adjust.method)

  # Filter by genelist
  if (!is.null(genelist)) {
    stopifnot(is.character(genelist), length(genelist) > 0, all(genelist %in% row.names(results)))
    results <- results[genelist, ]
  }

  # Filter by p.value
  if (!is.null(p.value)) {
    stopifnot(is.numeric(p.value))
    stopifnot(length(p.value) == 1 && p.value >= 0 && p.value <= 1)
    results <- results[results$adj.P.Value <= p.value, ]
  }

  # Sort
  if (is.null(sort.by) || any(!sort.by %in% names(results)) || !is.character(sort.by)) {
    stop(
      paste("Sort by one of the following:", paste(names(results), collapse = ", "))
    )
  }
  results <- results[do.call(order, results[sort.by]), ]

  if (S4Vectors::metadata(object)$AIC) {
    results$AIC <- vapply(S4Vectors::mcols(object)$fit, function(x) x$AIC, numeric(1))
    results$BIC <- vapply(S4Vectors::mcols(object)$fit, function(x) x$BIC, numeric(1))
  }

  # Return only head number
  stopifnot(length(number) == 1, number >= 0)
  if (!is.infinite(number)) {
    stopifnot(is.numeric(number), as.integer(number) == number)
  }
  return(utils::head(results, n = number))
}

#' Grab Term from the Smooth vs Non-Smooth Table
#'
#' This function extracts a specified term from the smooth or non-smooth table in the fit results.
#'
#' @param coef Name of the term after fitting.
#' @param fit Fit result object.
#'
#' @return A 1*n data.frame containing the term's statistics.
#' @noRd
#' @keywords internal
filter_terms <- function(coef, fit) {
  is_smooth <- grepl("s\\(", coef)
  if (is_smooth) {
    return(fit$s.table[coef, ])
  } else {
    return(fit$p.table[coef, ])
  }
  
  return(1)
}

wrap_topTable1 <- function(fit, coef) {
  fit <- topTable1(fit, coef = coef)
  p <- fit[["P.Value"]]
}