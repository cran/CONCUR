#' @importFrom utils read.table
#' @include test.bt.R test.qt.R
.pv_func <- function(...,
                     pheno, 
                     phenoY,  
                     phenoType,  
                     X,  
                     kernel,  
                     verbose) {

  if (verbose) cat("calculating p-value\n\tverifying phenoType input\n")

  if (is.character(x = pheno)) {
    pheno <- tryCatch(expr = utils::read.table(file = pheno, 
                                               as.is = TRUE,
                                               header = TRUE),
                      error = function(e){
                                cat("unable to read", pheno,"\n")
                                stop(e$message)
                              })

  } else if (is.matrix(x = pheno)) {
    pheno <- as.data.frame(x = pheno)
  }

  outFileKernel <- NULL

  if (is.character(x = kernel)) {

    outFileKernel <- kernel

    kernel <- tryCatch(expr = utils::read.table(file = kernel, header = TRUE),
                       condition = function(e){
                                     cat("unable to read kernel\n")
                                     stop( e$message )
                                   })

    kernel <- data.matrix(frame = kernel)

    if (verbose) cat("\tread kernel from file", outFileKernel, "\n")
  }

  cnvID <- colnames(x = kernel)
  if (is.factor(x = pheno$ID)) pheno$ID <- as.character(x = levels(pheno$ID))[pheno$ID]

  # add diagonal elements for individuals with pheno data but no cnv data
  inBoth <- pheno$ID %in% cnvID

  if (!all(inBoth)) {
    mm <- match(x = cnvID, table = pheno$ID)

    if (any(is.na(x = mm))) stop("ID in cnv not present in pheno")

    minNonZero <- min(x = diag(x = kernel))
    nAdd <- sum(!inBoth)

    cnames <- colnames(x = kernel)

    kernel <- rbind(cbind(kernel, 
                          matrix(data = 0.0, 
                                 nrow = nrow(x = kernel),  
                                 ncol = nAdd)),
                    cbind(matrix(data = 0.0,
                                 nrow = nAdd,
                                 ncol = ncol(x = kernel)),
                          diag(x = minNonZero, nrow = nAdd, ncol = nAdd)))

    newnames <- c(cnames, pheno$ID[!inBoth])
    dimnames(x = kernel) <- list(newnames, newnames)

    if (verbose) {
      cat("\textended kernel to include individuals with no cnv data\n")
    }

    ok <- order(colnames(x = kernel))
    kernel <- kernel[ok,ok]
    cnvID <- colnames(x = kernel)

    if (!is.null(x = outFileKernel)) {

      tryCatch(expr = utils::write.table(x = kernel, file = outFileKernel),
               condition = function(e){
                             cat("unable to save kernel\n")
                             stop( e$message )
                           })
      if (verbose) cat("\tsaved new kernel to file", outFileKernel, "\n")

    }
  }


  if (!all(pheno$ID == cnvID)) {
    mm <- match(x = cnvID, table = pheno$ID)
    if (any(is.na(x = mm))) stop("ID in cnv not present in pheno")
    if (verbose) cat("\treorderd pheno to align with cnv\n")
    pheno <- pheno[mm,]
  }

  if (verbose) cat("\tverifying phenoY input\n")
  phenoY <- tryCatch(expr = pheno[,phenoY],
                     condition = function(x){
                                   cat("unable to identify", phenoY, "in pheno\n")
                                   stop(x$message)
                                 })

  rm(pheno)

  if (verbose) cat("\tverifying X input\n")
  if (is.character(x = X)) {
    X <- tryCatch(expr = utils::read.table(file = X, 
                                           as.is = TRUE,
                                           header = TRUE),
                  error = function(e){
                            cat("unable to read", X,"\n")
                            stop(e$message)
                          })
  } else if (is.matrix(x = X)) {
    X <- as.data.frame(x = X)
  } else if (!is.data.frame(x = X) & !is.null(x = X)) {
    stop('X must be NULL, a data.frame, or a path to the data file')
  }

  if (!is.null(x = X)) {
    if (!{"ID" %in% colnames(X)}) stop("ID must be a column of X data.frame")
  }
  if (is.factor(x = X$ID)) X$ID <- as.character(x = levels(X$ID))[X$ID]

  if (!is.null(x = X)) {

    if (!all(X$ID == cnvID)) {
      mm <- match(cnvID, X$ID)
      if (any(is.na(x = mm))) stop("ID in pheno not present in X")
      if (verbose) cat("\treorderd X to align with cnv\n")
      X <- X[mm,]
    }

    X <- X[,!{colnames(X) %in% "ID"},drop=FALSE]
  }

  if (phenoType == 'bin') {
    pv <- .test.bt(y = phenoY, 
                   K = kernel,  
                   X = X,  
                   verbose = verbose)
  } else {
    pv <- .test.qt(y = phenoY, K = kernel, X = X, verbose = verbose)
  }

  if (verbose) cat("p-value", pv, "\n")

  if (!is.null(x = outFileKernel)) kernel <- outFileKernel

  return( list("kernel" = kernel, "pValue" = pv) )
}
