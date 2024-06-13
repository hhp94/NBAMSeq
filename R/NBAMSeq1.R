#' @title Differential expression analysis based on negative binomial additive
#' model
#'
#' @description This function performs differential expression analysis based
#' on negative binomial additive model.
#' @param object a NBAMSeqDataSet object
#' @param gamma a number greater or equal to 1. Increase gamma to create
#' smoother models. Default gamma is 2.5. See \code{\link[mgcv]{gam}} for
#' details.
#' @param BPPARAM an argument provided to \code{\link{bplapply}}. See
#' \code{\link[BiocParallel]{register}} for details.
#' @param use_bam use bam from mgcv instead of gam. bam can be faster for large data.
#' @param AIC generate AIC/BIC for each gene or not.
#' @param alpha_bound set upper bound for gene wise dispersion parameter. Can help if there are errors. 10 is sufficiently big.
#' @param ... additional arguments provided to \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}
#' @export
#' @importFrom mgcv bam gam s nb negbin
#' @import SummarizedExperiment S4Vectors DESeq2
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom methods is new
#' @importFrom stats coef formula update AIC BIC
#' @return a NBAMSeqDataSet object
#' @references Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of
#' fold change and dispersion for RNA-seq data with DESeq2. Genome Biology,
#' 15:550. \url{https://doi.org/10.1186/s13059-014-0550-8}
#' @examples
#' gsd <- makeExample(n = 3, m = 10)
#' gsd <- NBAMSeq1(gsd)
NBAMSeq1 <- function(
    object, gamma = 2.5, BPPARAM = NULL, use_bam = FALSE, AIC = FALSE, alpha_bound = Inf,
    ...) {
  ## check input
  stopifnot(is(object, "NBAMSeqDataSet"))
  stopifnot(is.numeric(gamma))
  stopifnot(length(gamma) == 1)
  stopifnot("'gamma' should be greater or equal to 1." = gamma >= 1)
  stopifnot(is.logical(use_bam))
  stopifnot(length(use_bam) == 1)
  stopifnot(is.logical(AIC))
  stopifnot(length(AIC) == 1)
  stopifnot(alpha_bound > 0, length(alpha_bound) == 1)
  metadata(object)$AIC <- AIC

  if (is.null(BPPARAM)) {
    BPPARAM <- bpparam("SerialParam")
  }
  if (use_bam) {
    fns <- bam
  } else {
    fns <- gam
  }
  ## construct a DESeqDataSet object
  ddsdesign <- reformulate(all.vars(getDesign(object)))
  dds <- DESeqDataSetFromMatrix(
    countData = assay(object),
    colData = colData(object),
    design = ddsdesign,
    ignoreRank = TRUE # Set ignoreRank = TRUE because we can fit mixed model with gam
  )

  if ("sizeFactors" %in% names(colData(dds))) {
    logsf <- log(colData(object)$sizeFactors)
    object$logsf <- logsf ## save logsf in object
    sizeFactors(dds) <- colData(object)$sizeFactors
  } else {
    ## estimate size factors
    dds <- estimateSizeFactors(dds)
    logsf <- log(sizeFactors(dds))
    colData(object)$sizeFactors <- sizeFactors(dds)
    object$logsf <- logsf ## save logsf in object
  }

  dat <- as.data.frame(colData(object))
  dat$logsf <- logsf
  formula_offset <- update(getDesign(object), y ~ . + offset(logsf))

  message("\nEstimating smoothing parameters and gene-wise dispersions")
  gamGeneEst <- bplapply(
    seq_len(nrow(object)), \(i) {
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
  stopifnot("GAM for all genes failled. Check the formula and data" = length(gamGeneEst) > 0)
  passed_genes <- vapply(gamGeneEst, function(x) {x$gene}, "character")
  dds <- dds[passed_genes, ]
  object <- object[passed_genes, ]
  ##  get gam gene wise dispersion estimates and save them in mcols(dds)
  mcols(dds)$dispGeneEst <- (1 / vapply(gamGeneEst, function(x) {x$theta}, 1))

  if (!is.infinite(alpha_bound)) {
    ##  bound gene wise dispersion
    maxDisp <- pmax(alpha_bound, ncol(dds))
    mcols(dds)$dispGeneEst <- pmin(mcols(dds)$dispGeneEst, maxDisp)
  }

  ##  fit dispersion trend via `estimateDispersionsFit` function in DESeq2
  message("Estimating dispersion trend")
  dds <- estimateDispersionsFit(dds)

  ##  get gam mu estimates and save them in assays(dds)[["mu"]]
  muhat <- t(vapply(gamGeneEst, function(x) x$muhat, rep(1, ncol(dds))))
  colnames(muhat) <- colnames(dds)
  rownames(muhat) <- rownames(dds)
  assays(dds)[["mu"]] <- muhat

  ##  MAP dispersion estimates
  message("Estimating MAP dispersion")
  dds <- tryCatch(
    expr = {
      estimateDispersionsMAP(dds)
    },
    error = function(e) {
      ## avoid possible matrix singular error in DESeq2 C++ code
      assays(dds)[["mu"]] <- muhat + 1e-6
      estimateDispersionsMAP(dds)
    }
  )

  ind <- which(is.na(mcols(dds)$dispMAP))
  if (length(ind) > 0) {
    mcols(dds)$dispMAP[ind] <- mcols(dds)$dispGeneEst[ind]
  }

  gamDispMAP <- mcols(dds)$dispMAP

  message("Estimating model coefficients")
  gamFinal <- bplapply(
    seq_len(nrow(object)),
    \(i) {
      gamFit2(
        i,
        gam_fns = fns,
        dat = dat,
        object = object,
        formula_offset = formula_offset,
        gamGeneEst = gamGeneEst,
        gamDispMAP = gamDispMAP,
        AIC = AIC,
        ...
      )
    },
    BPPARAM = BPPARAM
  )

  ## remove models that failed and align the data by row names
  gamFinal <- gamFinal[lengths(gamFinal) != 0]
  stopifnot("GAM for all genes second step failed. Check the formula and data" = length(gamFinal) > 0)
  passed_genes_final <- vapply(gamFinal, function(x) {x$gene_id}, FUN.VALUE = "character")
  names(gamFinal) <- passed_genes_final
  dds <- dds[passed_genes_final, ]
  object <- object[passed_genes_final, ]

  ## Add the summary tables as a list column
  mcols(dds)$fit <- I(gamFinal)
  mcols(object) <- mcols(dds)
  metadata(object)$fitted <- TRUE
  message("Done!")

  return(object)
}

#' Gene-wise GAM model
#'
#' Wrapper around gam or bam that allows you to modify gamma
#'
#' @param gam_fns bam or gam
#' @param formula formula to fit on each gene
#' @param gamma_value see ?mgcv::gam
#' @param data data frame
#' @param ... argument passed to gam
#'
#' @noRd
#' @keywords internal
#' @return a gam fit
fit_gam1 <- function(gam_fns, formula, gamma_value, data, ...) {
  fit <- gam_fns(
    formula = formula,
    family = nb(link = "log"),
    method = "REML",
    gamma = gamma_value,
    data = data,
    ...
  )
  fit$gamma <- gamma_value
  return(fit)
}

#' Gene-wise GAM model wrapper
#'
#' Wrapper around fit_gam1 that remove models that failed to fit or failed to converge
#'
#' @param i iteration
#' @param object NBAMSeqDataSet object
#' @param formula_offset formula with the offset
#' @param gamma gamma value
#' @param dat data frame
#' @param fns gam or bam
#' @param ... argument passed to gam
#'
#' @noRd
#' @keywords internal
#' @return a gam fit or NULL
gamFit1 <- function(i, object, formula_offset, gamma, dat, fns, ...) {
  dat$y <- assay(object)[i, ]

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
      coef = coef(gamfit),
      muhat = gamfit$fitted.values,
      outIter = gamfit$outer.info$iter,
      gamma = gamfit$gamma,
      gene = row.names(assay(object))[i]
    )
  )
}

#' Fit Gam Step 2
#'
#' @param gam_fns bam or gam
#' @param formula form
#' @param data data
#' @param start step 1 param
#' @param sp step 1 param
#' @param theta step 1 param
#' @param ... bam or gam args
#'
#' @keywords internal
#' @noRd
fit_gam2 <- function(gam_fns, formula, data, start, sp, theta, ...) {
  fit <- gam_fns(
    formula = formula,
    family = negbin(theta = theta, link = "log"),
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
#' @param i iteration
#' @param gam_fns bam or gam
#' @param dat data
#' @param object NBAMSeqDataSet object
#' @param formula_offset formula
#' @param gamGeneEst gamGeneEst obj
#' @param gamDispMAP gamDispMAP obj
#' @param AIC add AIC and BIC or not
#' @param ... args passed to gam or bam
#'
#' @noRd
#' @keywords internal
#' @return list of tables that contain p values
gamFit2 <- function(i, gam_fns, dat, object, formula_offset, gamGeneEst, gamDispMAP, AIC, ...) {
  dat$y <- assay(object)[i, ] # ith gene count
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
    gene_id = row.names(object)[i],
    p.table = sum_obj$p.table,
    s.table = sum_obj$s.table
  )

  if (AIC) {
    results <- c(results, list(AICnonlin = AIC(gamFinalFit), BICnonlin = BIC(gamFinalFit)))
  }

  return(results)
}
