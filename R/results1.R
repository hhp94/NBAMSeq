#' @title Pulling out result
#'
#' @description This function pulls out result from NBAMSeqDataSet object
#' returned by \code{\link{NBAMSeq}}
#' @param object a NBAMSeqDataSet object returned by \code{\link{NBAMSeq}}
#' @param name the name of nonlinear variable or continuous linear variable
#' @param contrast a character of length 3. 1st element: name of factor
#' variable; 2nd element: name of numerator level; 3rd element: name of
#' denominator level. contrast = c("group", "treatment", "control") means
#' comparing treatment vs control for group variable.
#' @param indepfilter either TRUE or FALSE indicating whether independent
#' filtering should be performed. Default is TRUE.
#' @param alpha significant threshold for declaring genes as differentially
#' expressed. Default is 0.1.
#' @param pAdjustMethod pvalue adjustment method. Default is "BH". See
#' \code{\link[stats]{p.adjust}} for details.
#' @param BPPARAM an argument provided to \code{\link{bplapply}}. See
#' \code{\link[BiocParallel]{register}} for details.
#' @param ... additional arguments provided to \code{pvalueAdjustment} function
#' in DESeq2. See \code{\link[DESeq2]{results}} for details.
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom stats lowess p.adjust quantile
#' @importFrom genefilter filtered_p
#' @export
#' @return a DataFrame which contains the result
#' @references Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of
#' fold change and dispersion for RNA-seq data with DESeq2. Genome Biology,
#' 15:550. \url{https://doi.org/10.1186/s13059-014-0550-8}
#' @examples
#' gsd <- makeExample(n = 3, m = 10)
#' gsd <- NBAMSeq1(gsd)
#' res <- results1(gsd, name = "pheno")
results1 <- function(object, name, contrast, indepfilter = TRUE, alpha = 0.1,
                     pAdjustMethod = "BH", BPPARAM = NULL, ...) {
  ## check input
  stopifnot(is(object, "NBAMSeqDataSet"))
  ## check whether NBAMSeq has been run before
  if (!metadata(object)$fitted) {
    stop("NBAMSeq function should be run before calling results function")
  }
  if (is.null(BPPARAM)) {
    BPPARAM <- bpparam("SerialParam")
  }
  stopifnot(is.logical(indepfilter))
  stopifnot(length(indepfilter) == 1)

  stopifnot(is.numeric(alpha))
  stopifnot(length(alpha) == 1)
  stopifnot(alpha > 0 & alpha < 1)

  if (missing(name) & missing(contrast)) {
    stop("Either 'name' or 'contrast' argument should be provided.")
  } else if (!missing(name) & !missing(contrast)) {
    stop(" 'name' and 'contrast' arguments cannot be both provided.")
  } else if (!missing(name)) {
    stopifnot(is.character(name))
    stopifnot(length(name) == 1)
    stopifnot(name %in% all.vars(getDesign(object)))

    ## check if name is a factor
    if (is.factor(colData(object)[[name]])) {
      stop(" 'name' should be a continuous variable. For factors,
                please use 'contrast' argument.")
    }

    ## check the variable is linear or nonlinear
    flag <- grepl(
      paste0("s\\(", name, "\\)"),
      as.character(getDesign(object))[2]
    )
    if (flag) {
      clnm <- c("baseMean", paste0(c("edf_", "Chisq_", "PValue_"), name))
      resdf <- mcols(object)[clnm]
      names(resdf) <- c("baseMean", "edf", "stat", "pvalue")
    } else {
      clnm <- c("baseMean", name, paste0(c("SE_", "PValue_"), name))
      resdf <- mcols(object)[clnm]
      resdf$stat <- resdf[[name]] / resdf[[paste0("SE_", name)]]
      names(resdf) <- c("baseMean", "coef", "SE", "pvalue", "stat")
      resdf <- resdf[c("baseMean", "coef", "SE", "stat", "pvalue")]
    }
  } else {
    if (!is.character(contrast) | length(contrast) != 3) {
      stop(" 'contrast' should be a character of length 3.")
    }
    target <- contrast[1]
    firstlevel <- contrast[2]
    secondlevel <- contrast[3]
    alllevels <- levels(colData(object)[[target]])

    if (firstlevel == secondlevel) {
      stop("2nd and 3rd element in constrast should be different.")
    }
    if (!target %in% names(colData(object))) {
      stop("1st element in contrast should be a variable in colData.")
    }
    if (!is.factor(colData(object)[[target]])) {
      stop("The variable in contrast should be a factor. For continuous
                variables, please use 'name' argument.")
    }
    if (!firstlevel %in% alllevels) {
      stop("2nd element in contrast should be an appropriate level.")
    }
    if (!secondlevel %in% alllevels) {
      stop("3rd element in contrast should be an appropriate level.")
    }

    baselevel <- levels(colData(object)[[target]])[1]
    if (firstlevel != baselevel & secondlevel != baselevel) {
      ## if both first level and second level are not base level,
      ## relevel the factor and refit the model
      dat <- data.frame(colData(object))
      dat[, target] <- factor(dat[, target],
        levels = c(secondlevel, setdiff(alllevels, secondlevel))
      )
      formula_offset <- update(getDesign(object), y ~ . + offset(logsf))
      temp <- paste0("smooth_", all.vars(getDesign(object)))
      smoothcol <- temp[temp %in% names(mcols(object))]

      gamFit2 <- function(i) {
        dat$y <- assay(object)[i, ] ## ith gene count
        gamFinalFit <- gam(formula_offset,
          family = negbin(
            theta = 1 / mcols(object)[["dispMAP"]][3],
            link = "log"
          ), method = "REML",
          sp = unlist(mcols(object)[smoothcol][i, ]), data = dat
        )
        summary(gamFinalFit)$p.table[paste0(target, firstlevel), ]
      }

      reslist <- bplapply(seq_len(nrow(object)), gamFit2, BPPARAM = BPPARAM)

      resdf <- DataFrame(matrix(unlist(reslist), ncol = 4, byrow = TRUE))
      names(resdf) <- c("coef", "SE", "stat", "pvalue")
      resdf <- cbind(mcols(object)["baseMean"], resdf)
    } else {
      if (secondlevel == baselevel) {
        clnm <- c("baseMean", paste0(
          c("", "SE_", "PValue_"),
          paste0(target, "_", firstlevel, "_vs_", secondlevel)
        ))
      } else {
        clnm <- c("baseMean", paste0(
          c("", "SE_", "PValue_"),
          paste0(target, "_", secondlevel, "_vs_", firstlevel)
        ))
      }
      resdf <- mcols(object)[clnm]
      names(resdf) <- c("baseMean", "coef", "SE", "pvalue")
      resdf$stat <- resdf[["coef"]] / resdf[["SE"]]
      resdf <- resdf[c("baseMean", "coef", "SE", "stat", "pvalue")]
      if (firstlevel == baselevel) {
        resdf[["coef"]] <- -resdf[["coef"]]
        resdf[["stat"]] <- -resdf[["stat"]]
      }
    }
  }

  ## perform independent filtering
  resdf <- pvalueAdjustment(resdf,
    independentFiltering = indepfilter,
    alpha = alpha, pAdjustMethod = pAdjustMethod, ...
  )
  rownames(resdf) <- rownames(assay(object))

  if (metadata(object)$AIC) {
    AICBICnms <- c(
      "AIC", "BIC", "AICnonlin", "AIClin", "AICintercept",
      "BICnonlin", "BIClin", "BICintercept"
    )
    AICBICcolind <- which(names(mcols(object)) %in% AICBICnms)
    resdf <- cbind(resdf, DataFrame(mcols(object)[, AICBICcolind]))
  }

  return(resdf)
}
