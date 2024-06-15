#' Differential Expression Analysis Based on Negative Binomial Additive Model
#'
#' @description This function performs differential expression analysis based
#' on a negative binomial additive model.
#' @param object A [NBAMSeqDataSet()] object.
#' @param gamma A number greater than or equal to 1. Increase gamma to create
#' smoother models. Default gamma is 2.5. See [mgcv::gam()] for details.
#' @param BPPARAM An argument provided to [BiocParallel::bplapply()]. See
#' [BiocParallel::register()] for details.
#' @param use_bam Logical, whether to use [mgcv::bam()] instead of [mgcv::gam()]. `bam` can be faster for large data.
#' @param AIC Logical, whether to generate AIC/BIC for each gene.
#' @param alpha_bound Numeric, set upper bound for gene-wise dispersion parameter. Can help if there are errors. A value of 10 is sufficiently large.
#' @param save_model Logical, whether to save the final [mgcv::gam()] model for each gene. Warning: will take up a lot of memory.
#' @param verbose Logical, whether to print verbose output.
#' @param ... Additional arguments provided to [mgcv::gam()] or [mgcv::bam()].
#' @return A [NBAMSeqDataSet()] object with `metadata(object)$fitted == TRUE`.
#' @export
#' @references Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of
#' fold change and dispersion for RNA-seq data with DESeq2. Genome Biology,
#' 15:550. <https://doi.org/10.1186/s13059-014-0550-8>
#' @examples
#' gsd <- NBAMSeq::makeExample(n = 3, m = 10)
#' gsd <- NBAMSeq::NBAMSeq1(gsd)
NBAMSeq1 <- function(
    object, gamma = 2.5, BPPARAM = NULL, use_bam = FALSE, AIC = FALSE,
    alpha_bound = Inf, save_model = FALSE, verbose = FALSE,
    ...) {
  ## check input
  stopifnot(methods::is(object, "NBAMSeqDataSet"))
  stopifnot(is.numeric(gamma))
  stopifnot(length(gamma) == 1)
  stopifnot("'gamma' should be greater or equal to 1." = gamma >= 1)
  stopifnot(is.logical(use_bam))
  stopifnot(length(use_bam) == 1)
  stopifnot(is.logical(AIC))
  stopifnot(length(AIC) == 1)
  stopifnot(alpha_bound > 0, length(alpha_bound) == 1)
  stopifnot(is.logical(save_model))
  stopifnot(length(save_model) == 1)
  stopifnot(is.logical(verbose))
  stopifnot(length(verbose) == 1)

  S4Vectors::metadata(object)$AIC <- AIC
  if (length(row.names(object)) > unique(length(row.names(object)))) {
    stop("duplicated gene_id are not allowed, check row.names(object)")
  }
  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::bpparam("SerialParam")
  }
  if (use_bam) {
    fns <- mgcv::bam
  } else {
    fns <- mgcv::gam
  }
  ## construct a DESeqDataSet object
  ddsdesign <- stats::reformulate(all.vars(getDesign(object)))
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = SummarizedExperiment::assay(object),
    colData = SummarizedExperiment::colData(object),
    design = ddsdesign,
    ignoreRank = TRUE # Set ignoreRank = TRUE because we can fit mixed model with gam
  )

  if ("sizeFactors" %in% names(SummarizedExperiment::colData(dds))) {
    logsf <- log(SummarizedExperiment::colData(object)$sizeFactors)
    object$logsf <- logsf ## save logsf in object
    DESeq2::sizeFactors(dds) <- SummarizedExperiment::colData(object)$sizeFactors
  } else {
    ## estimate size factors
    dds <- DESeq2::estimateSizeFactors(dds)
    logsf <- log(DESeq2::sizeFactors(dds))
    SummarizedExperiment::colData(object)$sizeFactors <- DESeq2::sizeFactors(dds)
    object$logsf <- logsf ## save logsf in object
  }

  dat <- as.data.frame(SummarizedExperiment::colData(object))
  dat$logsf <- logsf
  formula_offset <- stats::update(getDesign(object), y ~ . + stats::offset(logsf))

  if (verbose) {
    message("\nEstimating smoothing parameters and gene-wise dispersions")
  }

  gamGeneEst <- BiocParallel::bplapply(
    seq_len(nrow(object)), function(i) {
      gamFit1(
        i = i,
        object = object,
        formula_offset = formula_offset,
        gamma = gamma,
        dat = dat,
        fns = fns,
        ...
      )
    },
    BPPARAM = BPPARAM
  )

  ## remove models that failed
  gamGeneEst <- gamGeneEst[lengths(gamGeneEst) != 0]
  stopifnot("GAM for all genes failed. Check the formula and data" = length(gamGeneEst) > 0)
  passed_genes <- vapply(gamGeneEst, function(x) x$gene, character(1))
  dds <- dds[passed_genes, ]
  object <- object[passed_genes, ]
  ##  get gam gene wise dispersion estimates and save them in S4Vectors::mcols(dds)
  S4Vectors::mcols(dds)$dispGeneEst <- (1 / vapply(gamGeneEst, function(x) x$theta, numeric(1)))

  if (!is.infinite(alpha_bound)) {
    ##  bound gene wise dispersion
    maxDisp <- pmax(alpha_bound, ncol(dds))
    S4Vectors::mcols(dds)$dispGeneEst <- pmin(S4Vectors::mcols(dds)$dispGeneEst, maxDisp)
  }

  ##  fit dispersion trend via `DESeq2::estimateDispersionsFit` function
  if (verbose) {
    message("Estimating dispersion trend")
  }
  dds <- DESeq2::estimateDispersionsFit(dds)

  ##  get gam mu estimates and save them in SummarizedExperiment::assays(dds)[["mu"]]
  muhat <- t(vapply(gamGeneEst, function(x) x$muhat, rep(1, ncol(dds))))
  colnames(muhat) <- colnames(dds)
  rownames(muhat) <- rownames(dds)
  SummarizedExperiment::assays(dds)[["mu"]] <- muhat

  ##  MAP dispersion estimates
  if (verbose) {
    message("Estimating MAP dispersion")
  }
  dds <- tryCatch(
    DESeq2::estimateDispersionsMAP(dds),
    error = function(e) {
      ## avoid possible matrix singular error in DESeq2 C++ code
      SummarizedExperiment::assays(dds)[["mu"]] <- muhat + 1e-6
      DESeq2::estimateDispersionsMAP(dds)
    }
  )

  ind <- which(is.na(S4Vectors::mcols(dds)$dispMAP))
  if (length(ind) > 0) {
    S4Vectors::mcols(dds)$dispMAP[ind] <- S4Vectors::mcols(dds)$dispGeneEst[ind]
  }

  gamDispMAP <- S4Vectors::mcols(dds)$dispMAP

  if (verbose) {
    message("Estimating model coefficients")
  }
  gamFinal <- BiocParallel::bplapply(
    seq_len(nrow(object)),
    function(i) {
      gamFit2(
        i,
        gam_fns = fns,
        dat = dat,
        object = object,
        formula_offset = formula_offset,
        gamGeneEst = gamGeneEst,
        gamDispMAP = gamDispMAP,
        AIC = AIC,
        save_model = save_model,
        ...
      )
    },
    BPPARAM = BPPARAM
  )

  ## remove models that failed and align the data by row names
  gamFinal <- gamFinal[lengths(gamFinal) != 0]
  stopifnot(
    "GAM for all genes second step failed. Check the formula and data" =
      length(gamFinal) > 0
  )
  passed_genes_final <- vapply(gamFinal, function(x) x$gene_id, character(1))
  names(gamFinal) <- passed_genes_final
  dds <- dds[passed_genes_final, ]
  object <- object[passed_genes_final, ]

  ## Add the summary tables as a list column
  S4Vectors::mcols(dds)$fit <- I(gamFinal)
  S4Vectors::mcols(object) <- S4Vectors::mcols(dds)
  S4Vectors::metadata(object)$fitted <- TRUE

  if (!save_model) {
    for (i in row.names(object)) {
      S4Vectors::mcols(object)$fit[[i]]$gamFinalFit <- NULL
    }
  }

  if (verbose) {
    message("Done!")
  }

  return(object)
}

#' Gene-wise GAM Model
#'
#' Wrapper around [mgcv::gam()] or [mgcv::bam()] that allows you to modify gamma.
#'
#' @param gam_fns Function, either [mgcv::bam()] or [mgcv::gam()].
#' @param formula The formula to fit on each gene.
#' @param gamma_value Numeric, the gamma value. See [mgcv::gam()] for details.
#' @param data Data frame, the data to fit.
#' @param ... Additional arguments passed to [mgcv::gam()] or [mgcv::bam()].
#'
#' @return A GAM fit.
#' @noRd
#' @keywords internal
fit_gam1 <- function(gam_fns, formula, gamma_value, data, ...) {
  fit <- gam_fns(
    formula = formula,
    family = mgcv::nb(link = "log"),
    method = "REML",
    gamma = gamma_value,
    data = data,
    ...
  )
  fit$gamma <- gamma_value
  return(fit)
}

#' Gene-wise GAM Model Wrapper
#'
#' Wrapper around [fit_gam1()] that removes models that failed to fit or failed to converge.
#'
#' @param i Integer, the iteration.
#' @param object A [NBAMSeqDataSet] object.
#' @param formula_offset The formula with the offset.
#' @param gamma Numeric, the gamma value.
#' @param dat Data frame, the data to fit.
#' @param fns Function, either [mgcv::bam()] or [mgcv::gam()].
#' @param ... Additional arguments passed to [mgcv::gam()] or [mgcv::bam()].
#'
#' @return A GAM fit or `NULL`.
#' @keywords internal
#' @noRd
gamFit1 <- function(i, object, formula_offset, gamma, dat, fns, ...) {
  dat$y <- SummarizedExperiment::assay(object)[i, ]

  gamfit <- tryCatch(
    expr = fit_gam1(fns, formula_offset, gamma, dat, ...),
    error = function(e) {
      tryCatch(
        fit_gam1(fns, formula_offset, 1, dat, ...),
        error = function(e) NULL,
        warning = function(w) NULL
      )
    },
    warning = function(w) {
      tryCatch(
        fit_gam1(fns, formula_offset, 1, dat, ...),
        error = function(e) NULL,
        warning = function(w) NULL
      )
    }
  )

  if (is.null(gamfit)) {
    return(NULL)
  }

  return(
    list(
      theta = gamfit$family$getTheta(TRUE),
      sp = gamfit$sp,
      coef = stats::coef(gamfit),
      muhat = gamfit$fitted.values,
      outIter = gamfit$outer.info$iter,
      gamma = gamfit$gamma,
      gene = row.names(SummarizedExperiment::assay(object))[i]
    )
  )
}

#' Fit GAM Step 2
#'
#' Fits a generalized additive model (GAM) using either [mgcv::bam()] or [mgcv::gam()].
#'
#' @param gam_fns Function, either [mgcv::bam()] or [mgcv::gam()].
#' @param formula The formula to fit.
#' @param data Data frame, the data to fit.
#' @param start List, initial values for the parameters.
#' @param sp Numeric vector, smoothing parameters.
#' @param theta Numeric, dispersion parameter for the negative binomial family.
#' @param ... Additional arguments passed to [mgcv::gam()] or [mgcv::bam()].
#'
#' @return A fitted GAM model.
#' @keywords internal
#' @noRd
fit_gam2 <- function(gam_fns, formula, data, start, sp, theta, ...) {
  fit <- gam_fns(
    formula = formula,
    family = mgcv::negbin(theta = theta, link = "log"),
    method = "REML",
    sp = sp,
    start = start,
    data = data,
    ...
  )
  return(fit)
}

#' Second Step of GAM Fit
#'
#' Fits a generalized additive model (GAM) in the second step using either [mgcv::bam()] or [mgcv::gam()].
#'
#' @param i Integer, the iteration.
#' @param gam_fns Function, either [mgcv::bam()] or [mgcv::gam()].
#' @param dat Data frame, the data to fit.
#' @param object A [NBAMSeqDataSet()] object.
#' @param formula_offset The formula with the offset.
#' @param gamGeneEst List, the GAM gene estimates.
#' @param gamDispMAP Numeric vector, the GAM dispersion MAP estimates.
#' @param AIC Logical, whether to add AIC and BIC or not.
#' @param save_model Logical, whether to save the final [mgcv::gam()] model.
#' @param ... Additional arguments passed to [mgcv::gam()] or [mgcv::bam()].
#'
#' @return A list of tables that contain p-values.
#' @keywords internal
#' @noRd
gamFit2 <- function(
    i,
    gam_fns,
    dat,
    object,
    formula_offset,
    gamGeneEst,
    gamDispMAP,
    AIC,
    save_model, ...) {
  dat$y <- SummarizedExperiment::assay(object)[i, ] # ith gene count
  gene_id <- row.names(object)[i]

  gamFinalFit <- tryCatch(
    expr = fit_gam2(
      gam_fns = gam_fns,
      formula = formula_offset,
      data = dat,
      start = gamGeneEst[[i]]$coef, # initial coefficients
      sp = gamGeneEst[[i]]$sp,
      theta = (1 / gamDispMAP[i]),
      ...
    ),
    error = function(e) NULL,
    warning = function(w) NULL
  )

  if (is.null(gamFinalFit)) {
    return(NULL)
  }

  sum_obj <- summary(gamFinalFit)

  results <- list(
    gene_id = gene_id,
    p.table = sum_obj$p.table,
    s.table = sum_obj$s.table
  )

  if (AIC) {
    results <- c(results, list(AIC = stats::AIC(gamFinalFit), BIC = stats::BIC(gamFinalFit)))
  }

  if (save_model) {
    results <- c(results, list(gamFinalFit = gamFinalFit))
  }

  return(results)
}

#' Get Terms from Formula
#'
#' Extracts the terms from a formula object.
#'
#' @param formula A formula object.
#'
#' @return A character vector of terms.
#' @keywords internal
#' @noRd
get_term.labels <- function(formula) {
  return(attr(stats::terms(formula), "term.labels"))
}

wrap_results <- function(fit, ...) {
  r <- results(fit, ...)
  r$P.Value <- r$pvalue
  as.data.frame(r)[, c("edf", "P.Value")]
}

# Version with `joint_test`
# NBAMSeq1 <- function(
#     object, gamma = 2.5, BPPARAM = NULL, use_bam = FALSE, AIC = FALSE,
#     alpha_bound = Inf, save_model = FALSE, verbose = FALSE,
#     ...) {
#   joint_test = NULL # remove `joint_test`. Doesn't work well with fixed dispersion
#   ## check input
#   stopifnot(is(object, "NBAMSeqDataSet"))
#   stopifnot(is.numeric(gamma))
#   stopifnot(length(gamma) == 1)
#   stopifnot("'gamma' should be greater or equal to 1." = gamma >= 1)
#   stopifnot(is.logical(use_bam))
#   stopifnot(length(use_bam) == 1)
#   stopifnot(is.logical(AIC))
#   stopifnot(length(AIC) == 1)
#   stopifnot(alpha_bound > 0, length(alpha_bound) == 1)
#   stopifnot(is.logical(save_model))
#   stopifnot(length(save_model) == 1)
#   stopifnot(is.logical(verbose))
#   stopifnot(length(verbose) == 1)
#   old_save_model <- save_model
#   if (!is.null(joint_test)) {
#     stopifnot(
#       "joint_test should be a list of character of terms to be joint tested" =
#         is.list(joint_test) && all(sapply(joint_test, function(x) {
#           all(is.character(x))
#         }))
#     )
#     joint_test <- unique(lapply(joint_test, function(x) {
#       get_term.labels(as.formula(paste("~", paste(x, collapse = " + "))))
#     }))
#     design_terms <- get_term.labels(getDesign(object))
#     stopifnot("some terms in joint tests are not found in design" = all(unlist(joint_test) %in% design_terms))
#     save_model <- TRUE
#   }
#   metadata(object)$AIC <- AIC
#   if (length(row.names(object)) > unique(length(row.names(object)))) {
#     stop("duplicated gene_id are not allowed, check row.names(object)")
#   }
#   if (is.null(BPPARAM)) {
#     BPPARAM <- bpparam("SerialParam")
#   }
#   if (use_bam) {
#     fns <- bam
#   } else {
#     fns <- gam
#   }
#   ## construct a DESeqDataSet object
#   ddsdesign <- reformulate(all.vars(getDesign(object)))
#   dds <- DESeqDataSetFromMatrix(
#     countData = assay(object),
#     colData = colData(object),
#     design = ddsdesign,
#     ignoreRank = TRUE # Set ignoreRank = TRUE because we can fit mixed model with gam
#   )
#
#   if ("sizeFactors" %in% names(colData(dds))) {
#     logsf <- log(colData(object)$sizeFactors)
#     object$logsf <- logsf ## save logsf in object
#     sizeFactors(dds) <- colData(object)$sizeFactors
#   } else {
#     ## estimate size factors
#     dds <- estimateSizeFactors(dds)
#     logsf <- log(sizeFactors(dds))
#     colData(object)$sizeFactors <- sizeFactors(dds)
#     object$logsf <- logsf ## save logsf in object
#   }
#
#   dat <- as.data.frame(colData(object))
#   dat$logsf <- logsf
#   formula_offset <- update(getDesign(object), y ~ . + offset(logsf))
#
#   if (verbose) {
#     message("\nEstimating smoothing parameters and gene-wise dispersions")
#   }
#
#   gamGeneEst <- bplapply(
#     seq_len(nrow(object)), function(i) {
#       gamFit1(
#         i = i,
#         object = object,
#         formula_offset = formula_offset,
#         gamma = gamma,
#         dat = dat,
#         fns = fns,
#         ...
#       )
#     },
#     BPPARAM = BPPARAM
#   )
#
#   ## remove models that failed
#   gamGeneEst <- gamGeneEst[lengths(gamGeneEst) != 0]
#   stopifnot("GAM for all genes failled. Check the formula and data" = length(gamGeneEst) > 0)
#   passed_genes <- vapply(gamGeneEst, function(x) {
#     x$gene
#   }, "character")
#   dds <- dds[passed_genes, ]
#   object <- object[passed_genes, ]
#   ##  get gam gene wise dispersion estimates and save them in mcols(dds)
#   mcols(dds)$dispGeneEst <- (1 / vapply(gamGeneEst, function(x) {
#     x$theta
#   }, 1))
#
#   if (!is.infinite(alpha_bound)) {
#     ##  bound gene wise dispersion
#     maxDisp <- pmax(alpha_bound, ncol(dds))
#     mcols(dds)$dispGeneEst <- pmin(mcols(dds)$dispGeneEst, maxDisp)
#   }
#
#   ##  fit dispersion trend via `estimateDispersionsFit` function in DESeq2
#   if (verbose) {
#     message("Estimating dispersion trend")
#   }
#   dds <- estimateDispersionsFit(dds)
#
#   ##  get gam mu estimates and save them in assays(dds)[["mu"]]
#   muhat <- t(vapply(gamGeneEst, function(x) x$muhat, rep(1, ncol(dds))))
#   colnames(muhat) <- colnames(dds)
#   rownames(muhat) <- rownames(dds)
#   assays(dds)[["mu"]] <- muhat
#
#   ##  MAP dispersion estimates
#   if (verbose) {
#     message("Estimating MAP dispersion")
#   }
#   dds <- tryCatch(
#     estimateDispersionsMAP(dds),
#     error = function(e) {
#       ## avoid possible matrix singular error in DESeq2 C++ code
#       assays(dds)[["mu"]] <- muhat + 1e-6
#       estimateDispersionsMAP(dds)
#     }
#   )
#
#   ind <- which(is.na(mcols(dds)$dispMAP))
#   if (length(ind) > 0) {
#     mcols(dds)$dispMAP[ind] <- mcols(dds)$dispGeneEst[ind]
#   }
#
#   gamDispMAP <- mcols(dds)$dispMAP
#
#   if (verbose) {
#     message("Estimating model coefficients")
#   }
#   gamFinal <- bplapply(
#     seq_len(nrow(object)),
#     function(i) {
#       gamFit2(
#         i,
#         gam_fns = fns,
#         dat = dat,
#         object = object,
#         formula_offset = formula_offset,
#         gamGeneEst = gamGeneEst,
#         gamDispMAP = gamDispMAP,
#         AIC = AIC,
#         save_model = save_model,
#         ...
#       )
#     },
#     BPPARAM = BPPARAM
#   )
#
#   ## remove models that failed and align the data by row names
#   gamFinal <- gamFinal[lengths(gamFinal) != 0]
#   stopifnot(
#     "GAM for all genes second step failed. Check the formula and data" =
#       length(gamFinal) > 0
#   )
#   passed_genes_final <- vapply(gamFinal, function(x) {
#     x$gene_id
#   }, FUN.VALUE = "character")
#   names(gamFinal) <- passed_genes_final
#   dds <- dds[passed_genes_final, ]
#   object <- object[passed_genes_final, ]
#
#   ## Add the summary tables as a list column
#   mcols(dds)$fit <- I(gamFinal)
#   mcols(object) <- mcols(dds)
#   metadata(object)$fitted <- TRUE
#
#   ## Use recursion for joint test
#   if (!is.null(joint_test)) {
#     joint_results <- lapply(joint_test, function(x) {
#       test_joint(
#         x,
#         object,
#         gamma,
#         BPPARAM,
#         use_bam,
#         AIC,
#         alpha_bound,
#         joint_test,
#         save_model,
#         verbose,
#         ...
#       )
#     })
#     joint_results <- lapply(joint_results, function(x) {
#       simplify2array(x[row.names(object)])
#     })
#
#     for(i in seq_along(joint_results)) {
#       mcols(object)[[paste0("joint_results", i)]] <- joint_results[[i]]
#     }
#   }
#
#   if (!old_save_model) {
#     for (i in row.names(object)) {
#       mcols(object)$fit[[i]]$gamFinalFit <- NULL
#     }
#   }
#
#   if (verbose) {
#     message("Done!")
#   }
#
#   return(object)
# }

# #' Joint Test Function
# #'
# #' To be called in NBAMSeq1. Basically runs NBAMSeq1 with a reduced formula and
# #' then perform glrt for each gene
# #'
# #' @return a list of genes with p.value and NA if there's an error
# #' @noRd
# #' @keywords internal
# test_joint <- function(
#     test,
#     object,
#     gamma,
#     BPPARAM,
#     use_bam,
#     AIC,
#     alpha_bound,
#     joint_test,
#     save_model,
#     verbose,
#     ...) {
#   reduced_form <- update(
#     object@design,
#     as.formula(paste("~ . -", paste(test, collapse = " - ")))
#   )
#   object@design <- reduced_form
#   reduced_fit <- mcols(NBAMSeq1(
#     object,
#     gamma = gamma,
#     BPPARAM = BPPARAM,
#     use_bam = use_bam,
#     AIC = AIC,
#     alpha_bound = alpha_bound,
#     joint_test = NULL,
#     save_model = TRUE,
#     verbose = FALSE,
#     ...
#   ))$fit
#
#   if (is.null(names(reduced_fit))) {
#     genes <- vapply(fit, function(x) {
#       x$gene_id
#     }, "character")
#     names(reduced_fit) <- genes
#   }
#   stopifnot("unexpected gene id error" = all(names(reduced_fit) %in% row.names(object)))
#   results <- bplapply(
#     row.names(object),
#     function(i) {
#       tryCatch(
#         {
#           a <- anova(
#             reduced_fit[[i]]$gamFinalFit,
#             mcols(object)$fit[[i]]$gamFinalFit,
#             test = "Chisq"
#           )
#           pr_column <- grep("Pr\\(", names(a), value = TRUE)
#           a[[pr_column]][2]
#         },
#         warning = function(w) {
#           NA
#         },
#         error = function(e) {
#           NA
#         }
#       )
#     },
#     BPPARAM = BPPARAM
#   )
#   names(results) <- row.names(object)
#   return(results)
# }
