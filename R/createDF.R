.createDF <- function(cnv){

  m <- nrow(x = cnv)

  # {2*m x 4}
  id <- rep(x = 1L:m, each = 2L)
  df1 <- data.frame('id' = cnv$ID, 
                    'loc' = cnv$BP1,  
                    'chr' = cnv$CHR,  
                    'type' = cnv$TYPE)[id,]
  df1[{1L:m}*2L,'loc'] <- cnv$BP2
  df1[{1L:m}*2L,'type'] <- 2.0

  rownames(x = df1) <- NULL

  return( df1 )

}
