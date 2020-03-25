#' @importFrom dplyr left_join
#' @importFrom parallel makeCluster stopCluster
#' @include createDF.R uniqueCombinations.R popFunc.R segLength.R commonAUC.R
.deldup <- function(..., cnv, nCore, verbose) {

  if (verbose) cat("\tobtaining duplication/deletion matrices\n")

  m <- nrow(x = cnv)

  # {2*m x 4}
  df1 <- .createDF(cnv = cnv)

  # unique location/chromosome combinations {r x 2}
  uniqueLoc <- .uniqueCombinations(df = df1, verbose = verbose)

  # unique patient ids
  ids <- sort(x = unique(x = cnv$ID))

  n <- length(x = ids)
  r <- nrow(x = uniqueLoc)

  dup <- matrix(data = 0.0, nrow = r, ncol = n)
  del <- matrix(data = 0.0, nrow = r, ncol = n)

  len <- .segLength(uniqueLoc = uniqueLoc, verbose = verbose)

  tt <- numeric(length = r)
  cnt <- integer(length = r)

  for (i in 1L:length(x = ids)) {

    # cnv data for patient i
    cnvOfi <- df1[df1$id == ids[i],,drop = FALSE]

    # chromosomes of patient i found in others of the population
    common_chrs <- intersect(x = cnvOfi[,"chr"], 
                             y = df1[df1$id != ids[i],"chr"])

    if (length(x = common_chrs) == 0L) {
      # if no shared chromosomes, use 2.0 as type
      type <- rep(x = 2.0, times = r)
    } else {

      # transfer available type value(s) from cnvOfi to uniqueLoc
      uniqueLoci <- dplyr::left_join(x = uniqueLoc, 
                                     y = cnvOfi, 
                                     by = c('loc','chr'))

      type <- uniqueLoci$type

      tst <- uniqueLoci$chr %in% common_chrs

      # for chromosomes not shared with others in population, set to 2.0
      type[!tst] <- 2.0

      # for chromosomes shared with others in population, determine value
      type[tst] <- .popFunc(cnv = cnv[cnv$ID == ids[i],,drop = FALSE], 
                            common_chrs = common_chrs, 
                            s1 = uniqueLoci[tst,,drop = FALSE])
    }

    dif <- type - 2.0

    dup[,i] <- pmax( dif, 0.0)
    del[,i] <- pmax(-dif, 0.0)

    tst <- abs(x = dif) > 1e-8
    tt[tst] <- tt[tst] + len[tst]
    cnt <- cnt + tst

  }

  cnt[ cnt == 0L ] <- 1L

  len <- tt / cnt

  if (verbose) cat("\tobtaining kernel\n")

  if (nCore > 1L) {
    localCluster <- parallel::makeCluster(spec = nCore)
  } else {
    localCluster <- NULL
  }

  K <- .commonAUC(segLength = len, dup = dup, cluster = localCluster) +
       .commonAUC(segLength = len, dup = del, cluster = localCluster)

  if (nCore > 1L) {
    parallel::stopCluster(cl = localCluster)
  }

  dimnames(x = K) <- list(ids, ids)

  return( K )
}
