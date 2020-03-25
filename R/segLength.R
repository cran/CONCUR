#' @importFrom utils head
#' @importFrom utils tail
.segLength <- function(uniqueLoc, verbose) {

  dist <- diff(x = uniqueLoc$loc)*
          {utils::tail(x = uniqueLoc$chr, n = -1L) == 
           utils::head(x = uniqueLoc$chr, n = -1L)}

  return( c(dist, 0.0) )
}
