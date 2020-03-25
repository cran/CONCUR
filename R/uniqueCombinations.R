# unique location/chromosome combinations
#' @importFrom mgcv uniquecombs
.uniqueCombinations <- function(df, verbose) {

  if (!all(c("loc", "chr") %in% colnames(df))) stop("incorrect data provided")

  uniqueLoc <- mgcv::uniquecombs(x = df[,c('loc', 'chr')])
  uniqueLoc <- uniqueLoc[order(uniqueLoc$chr, uniqueLoc$loc),]
  r <- nrow(x = uniqueLoc)

  if (verbose) {
    cat("\tidentified", r, "unique chromosome/location combinations\n")
  }

  return( uniqueLoc )
}
