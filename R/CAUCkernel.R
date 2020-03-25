# make cAUC kernel function
#' @importFrom utils write.table
#' @include deldup.R
.CAUCkernel <- function(...,
                        cnv, 
                        nCore,
                        outFileKernel,
                        verbose) {

  if (verbose) cat('calculating CONCUR kernel\n')

  K <- .deldup(cnv = cnv, nCore = nCore, verbose = verbose)

  if (length(x = outFileKernel) == 0L || is.null(x = outFileKernel)) {
    return( K )
  } else if (is.character(x = outFileKernel)) {
    tryCatch(expr = utils::write.table(x = K, file = outFileKernel),
             condition = function(e){
                           cat("unable to save kernel\n")
                           stop( e$message )
                         })
    if (verbose) cat("\tsaved kernel in file", outFileKernel, "\n")
    return( outFileKernel )
  } else {
    stop('outFileKernel must be NULL or a character')
  }

}
