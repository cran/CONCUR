#' @importFrom stats lm
#' @importFrom stats residuals
#' @importFrom stats model.matrix
#' @importFrom utils read.table
#' @include compFunc.R
# Function adapted from CKAT package
.test.qt <- function (..., y, K, X = NULL, verbose) {

  if (verbose) cat("\tcontinuous phenotype\n")

  if (is.character(x = K)) {
    if (verbose) cat("\treading kernel from file", K, "\n")
    K <- utils::read.table(file = K, header = TRUE)
  }

  if (is.data.frame(x = K)) K <- data.matrix(frame = K)

  n <- length(x = y)

  if (is.null(x = X)) {
    X1 <-  matrix(data = rep(x = 1.0, times = length(x = y)), ncol = 1L)
  } else {
    X1 <- stats::model.matrix(~. , as.data.frame(X))
  }
    
  px <- ncol(x = X1)
  mod <- stats::lm(formula = y ~ X1-1)

  if (verbose) {
    cat("\tlm fit\n")
    print(mod)
  }

  res <- stats::residuals(object = mod)
  s2 <- sum(res^2) 
    
  D0 <- diag(x = length(x = y))

  P0 <- D0 - X1 %*% solve(a = crossprod(x = X1), b = t(X1))
  PKP <- P0 %*% K %*% P0
    
  q <- as.numeric(x = res %*% K %*% res / s2)

  return( .compFunc(X = PKP - q * P0, verbose = verbose) )

}

