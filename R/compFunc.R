#' @importFrom CompQuadForm davies
#' @importFrom CompQuadForm liu
.compFunc <- function(..., X, verbose) {

  ee <- eigen(x = X, symmetric = TRUE, only.values = TRUE)  
  lambda0 <- ee$values[abs(x = ee$values) >= 1e-10]

  if (verbose) cat("using davies method\n")
  p1 <- CompQuadForm::davies(q = 0.0, lambda = lambda0)$Qq

  if (p1 > 1.0 | p1 < 1e-8 ) {
    if (verbose) cat("result outside of expected range, using liu\n")
    p1 <- CompQuadForm::liu(q = 0.0, lambda = lambda0)
  }

  return( p1 )

}
