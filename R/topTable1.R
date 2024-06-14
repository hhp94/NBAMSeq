topTable1 <- function(
    object,
    coef = NULL,
    joint_results = NULL,
    number = 10,
    genelist = NULL,
    adjust.method = "BH",
    sort.by = "B",
    resort.by = NULL,
    p.value = 1) {
  ## check whether NBAMSeq has been run before
  stopifnot(is(object, "NBAMSeqDataSet"))
  if (!metadata(object)$fitted) {
    stop("NBAMSeq function should be run before calling results function")
  }
  all_terms <- unlist(lapply(
    mcols(object)$fit[[1]][c("p.table", "s.table")],
    \(x) {
      row.names(x)
    }
  ))
  if (is.null(joint_results) && (is.null(coef) || any(!(coef %in% all_terms)))) {
    message(
      paste("Please select only one of the following coef:", paste(all_terms, collapse = ", "))
    )
    return(invisible())
  }

  if (is.null(joint_results)) {
    stopifnot(length(coef) == 1 && is.character(coef))
    results <- as.data.frame(
      t(sapply(mcols(object)$fit, function(x) {
        filter_terms(coef = coef, x)
      }))
    )

    names(results)[which(grepl("Pr\\(|p-value", names(results)))] <- "P.Value"
  }

  if (!is.null(joint_results)) {
    stopifnot(
      "joint_results should be joint_results1, joint_results2, etc." =
        length(joint_results) == 1 && is.character(joint_results)
    )
    joint_columns <- grep("joint_results", names(mcols(object)), value = TRUE)
    if (!joint_results %in% joint_columns) {
      message(
        paste(
          "Please select only one joint test result:",
          paste(joint_columns, collapse = ", ")
        )
      )
      return(invisible())
    }

    results <- as.data.frame(mcols(object)[joint_results])
    names(results) <- "P.Value"
  }

  return(results)

  # stopifnot(is.numeric(p.value))
  # stopifnot(length(p.value) == 1, p.value >= 0 && p.value <= 1)
  # p.adjust.methods
  # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
  #   "fdr", "none")
}

filter_terms <- function(coef, fit) {
  is_smooth <- grepl("s\\(", coef)
  if (is_smooth) {
    fit$s.table[coef, ]
  } else {
    fit$p.table[coef, ]
  }
}
