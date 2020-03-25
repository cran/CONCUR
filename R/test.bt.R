# p-value function for binary
#' @importFrom stats glm
#' @importFrom stats model.matrix
#' @importFrom utils read.table
#' @include compFunc.R
# Function adapted from CKAT package
.test.bt <- function (..., y, K,  X = NULL, verbose, nms) {

  if (verbose) cat("\tbinary phenotype\n")

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
    
  glmfit <- stats::glm(formula = y ~ X1-1, family = "binomial")

  if (verbose) {
    cat("\tglm fit\n")
    print(glmfit)
  }
    
  mu <- glmfit$fitted.values
  res.wk <- glmfit$residuals
  res <- y - mu

  w <- mu*{1.0 - mu}
  sqrtw <- sqrt(x = w)
    
  adj <- sum({sqrtw * res.wk}^2) 
    
  DX12 <- sqrtw * X1
    
  qrX <- qr(x = DX12)
  Q <- qr.Q(qr = qrX)
  Q <- Q[, 1L:qrX$rank, drop=FALSE]

  P0 <- diag(x = length(x = y)) - tcrossprod(x = Q)
  DKD <- tcrossprod(x = sqrtw) * K

  tQK <- t(x = Q) %*% DKD
  QtQK <- Q %*% tQK 

  PKP1 <- DKD - QtQK - t(x = QtQK) + Q %*% (tQK %*% Q) %*% t(x = Q)

  q1 <- as.numeric(x = res %*% K %*% res)
  q1 <- q1 / adj

  return( .compFunc(X = PKP1 - q1 * P0, verbose = verbose) )

}

