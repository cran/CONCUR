.popFunc <- function(..., cnv, common_chrs, s1) {

  type <- s1$type

  for (j in common_chrs) {

    sisj <- s1$chr == j
    cisj <- cnv$CHR == j

    need <- s1$loc[sisj & is.na(x = type)]

    for (l in need) {

      tstRange <- {cnv$BP1 <= l} & {l <= cnv$BP2} & cisj

      if (any(tstRange)) {
        val <- sum(cnv$TYPE[tstRange])
      } else {
        val <- 2.0
      }

      type[sisj & {s1$loc == l}] <- val
    }

  }

  return( type )
}
