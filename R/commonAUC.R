#' @importFrom parallel parSapply
.commonAUC <- function(..., segLength, dup, cluster) {

  n <- ncol(x = dup)
  k <- nrow(x = dup)

  dup <- segLength*dup

  .iteration <- function(i, dup) {
                  return(sapply(X = {i+1L}:n, 
                                FUN = function(j, dupi, dup) {
                                        sum(pmin(dupi, dup[,j]))
                                      }, 
                                dupi = dup[,i],
                                dup = dup))
                 }

  if (is.null(x = cluster)) {
    vec <- unlist(x = sapply(X = 1L:{n-1L},
                             FUN = .iteration, 
                             dup = dup))
  } else {
    vec <- parallel::parSapply(cl = cluster, 
                               X = 1L:{n-1L}, 
                               FUN = .iteration,
                               dup = dup)
    vec <- unlist(x = c(vec))
  }

  res <- matrix(data = 0.0, nrow = n, ncol = n)
  res[lower.tri(x = res, diag = FALSE)] <- vec
  res[upper.tri(x = res)] <- t(x = res)[upper.tri(x = res)]

  diag(x = res) <- colSums(x = dup)

  return( res )

}
