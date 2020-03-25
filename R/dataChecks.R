.dataCheck <- function(cnv) {

  ids <- sort(x = unique(x = cnv$ID))

  msg <- NULL

  for (i in ids) {

    idMatch <- cnv$ID == i

    if (any(cnv[idMatch,"BP1"] %in% cnv[idMatch,"BP2"])) {

      idData <- cnv[idMatch,]

      chrs <- sort(x = unique(x = idData[,"CHR"]))

      for (j in chrs) {

        chrMatch <- idData[,"CHR"] == j

        chrData <- idData[chrMatch,]

        sameBP <- chrData[,"BP1"] == chrData[, "BP2"]

        if (any(sameBP)) {
          msg <- c(msg, paste("\nID", i, "CHR", j, "has same BP1 and BP2"))
          chrData <- chrData[!sameBP,,drop=FALSE]
          if (nrow(x = chrData) == 0L) next
        }

        overlap <- chrData[,"BP1"] %in% chrData[, "BP2"]

        if (any(overlap)) {
          msg <- c(msg, 
                   paste("\nID", i, "CHR", j, 
                         "CNV ends exactly where new CNV begins"))
        }

      }
    }
  }

  if (!is.null(x = msg)) {
    msg <- c(msg, "\nplease correct these data issues before proceeding")
    stop(msg, call. = FALSE)
  }

  return( NULL )
}
