#' Copy Number Profile Curve-Based Association Test
#'
#' Implements a kernel-based association test for CNV aggregate analysis
#'  in a certain genomic region (e.g., gene set, chromosome, or genome) that is
#'  robust to the within-locus and across-locus etiologoical heterogeneity, and 
#'  bypass the need to define a "locus" unit for CNVs.
#'
#' The CNV data must adhere to the following conditions:
#' \itemize{
#'   \item CNVs must be at least 1 unit long.
#'   \item CNVs cannot end at the exact location another begins
#' }
#' Violations of these conditions typically occur when data are rounded to
#'   a desired resolution.
#' For example
#'
#' \preformatted{
#'  ID CHR      BP1      BP2 TYPE
#'   1  13 10112087 10112414    3
#' }
#' becomes upon rounding to kilo
#' \preformatted{
#'  ID CHR   BP1   BP2 TYPE
#'   1  13 10112 10112    3     .
#' }
#' These cases should either be discarded or modified to be of length 1, e.g.,
#' \preformatted{
#'  ID CHR   BP1   BP2 TYPE
#'   1  13 10112 10113    3     .
#' }
#' As an example of condition 2 
#' \preformatted{
#'  ID CHR    BP1   BP2 TYPE
#'   1  13 100768 101100    3
#'   1  13 101100 101299    1
#' }
#' should be modified to one of
#' \preformatted{
#'  ID CHR    BP1   BP2 TYPE
#'   1  13 100768 101100    3
#'   1  13 101101 101299    1
#' }
#' or
#' \preformatted{
#'  ID CHR    BP1   BP2 TYPE
#'   1  13 100768 101099    3
#'   1  13 101100 101299    1     .
#' }
#' Additionally,
#' \preformatted{
#'  ID CHR    BP1   BP2 TYPE
#'   1  13 100768 101100    3
#'   1  13 101100 101299    3
#' }
#' should be combined as
#' \preformatted{
#'  ID CHR    BP1   BP2 TYPE
#'   1  13 100768 101299    3     .
#' }
#'
#' @param cnv A character or data.frame object. If character, the 
#'   name of the data file containing the CNV data (with a header). If 
#'   data.frame, the CNV data. The data must contain the following columns: 
#'   "ID", "CHR", "BP1", "BP2", "TYPE", where "ID" is a unique patient id,
#'   "CHR" is the CNV chromosome, "BP1" is the start location in base pairs
#'   or kilo-base pairs,
#'   "BP2" is the end location in base pairs or kilo-base pairs, and 
#'   "TYPE" is the CNV copy number.
#' @param X A character or data.frame object. If character, the 
#'   name of the data file containing the covariate data (with a header). If 
#'   data.frame, the covariate data. The data must contain a column titled 
#'   "ID" containing a unique patient id. This column must contain the
#'   same patient identifiers as contained in the CNV data specified in 
#'   input cnv. Categorical
#'   variables must be translated into design matrix format.
#' @param pheno A character or data.frame object. If character, the 
#'   name of the data file containing the phenotype data (with a header). If 
#'   data.frame, the phenotype data. The data must contain a column titled 
#'   "ID" containing a unique patient id. This column must contain the
#'   all of the patient identifiers contained in the CNV data specified in 
#'   input cnv. 
#' @param phenoY A character object. The column name in input pheno containing
#'   the phenotype of interest.
#' @param phenoType A character object. Must be one of of \{"bin", "cont"\} indicating
#'   if input phenoY (i.e., the phenotype of interest) is binary or continuous.
#' @param ... Ignored. Included to require named inputs.
#' @param nCore An integer object. If nCore > 1, package parallel is used to 
#'   calculate the kernel. Though the methods of package CompQuadForm dominate
#'   the time profile, setting nCore > 1L can improve computation times. 
#' @param outFileKernel A character object or NULL. If a character, the
#'   file in which the kernel is to be saved. If NULL, the kernel is returned
#'   by the function.
#' @param verbose A logical object. If TRUE, progress information is printed
#'   to the screen.
#'
#' @return A list containing the kernel (or its file name) and the p-value.
#'
#' @export
#'
#' @references Brucker, A., Lu, W., Marceau West, R., Yu, Q-Y., Hsiao, C. K.,
#'   Hsiao, T-H., Lin, C-H., Magnusson, P. K. E., Holloway, S. T., 
#'   Sullivan, P. F., Szatkiewicz, J. P., Lu, T-P., and
#'   Tzeng, J-Y. Association testing using Copy Number Profile Curves (CONCUR)
#'   enhances power in copy number variant analysis. <doi:10.1101/666875>.
#'
#' @include CAUCkernel.R pv_func.R dataChecks.R
#'
#' @examples
#'
#' data(cnvData)
#'
#' # binary phenoType
#' results <- concur(cnv = cnvData,
#'                   X = covData,
#'                   pheno = phenoData,
#'                   phenoY = 'PHEB',
#'                   phenoType = 'bin',
#'                   nCore = 1L,
#'                   outFileKernel = NULL,
#'                   verbose = TRUE)
#'
#' # continuous phenoType
#' results <- concur(cnv = cnvData,
#'                   X = covData,
#'                   pheno = phenoData,
#'                   phenoY = 'PHEC',
#'                   phenoType = 'cont',
#'                   nCore = 1L,
#'                   outFileKernel = NULL,
#'                   verbose = TRUE)
#'
concur <- function(cnv,
                   X,
                   pheno,
                   phenoY,
                   phenoType,
                   ...,
                   nCore = 1L,
                   outFileKernel = NULL,
                   verbose = TRUE) { 

  if (!is.character(x = phenoType)) stop("phenoType must be a character")
  phenoType <- tolower(x = phenoType)
  if (!{phenoType %in% c('bin','cont')}) {
    stop("phenoType must be one of {bin, cont}")
  }

  if (verbose) cat("\tverifying pheno input\n")

  phenoFile <- NULL
  if (is.character(x = pheno)) {
    phenoFile <- pheno

    pheno <- tryCatch(expr = utils::read.table(file = pheno, 
                                               as.is = TRUE,
                                               header = TRUE),
                      error = function(e){
                                cat("unable to read", pheno,"\n")
                                stop(e$message)
                              })
  } else if (is.matrix(x = pheno)) {
    pheno <- as.data.frame(x = pheno)
  }

  if (!is.data.frame(x = pheno)) {
    stop('pheno must be a data.frame or a path to the data file')
  }

  if (!{"ID" %in% colnames(x = pheno)}) stop("ID must be a column of pheno")

  if (!is.null(x = phenoFile)) {
    rm(pheno)
    pheno <- phenoFile
  }

  if (verbose) cat("processing copy number variant (cnv) data\n")

  # read/process cnv data  
  if (is.character(x = cnv)) {
    if (verbose) cat("\tcnv input provided as data file\n")
    cnv <- tryCatch(expr = read.table(file = cnv, 
                                      as.is = TRUE,
                                      header = TRUE),
                    error = function(e){
                              cat("unable to read", cnv,"\n")
                              stop(e$message)
                            })

  } else if (is.matrix(x = cnv)) {
    if (verbose) {
      cat("\tcnv input provided as matrix - converting to data.frame\n")
    }
    cnv <- as.data.frame(x = cnv)
  } else if (!is.data.frame(x = cnv)) {
    stop('cnv must be a data.frame or a path to the data file')
  }

  # ensure minimum data is present in cnv
  tst <- c("ID", "CHR", "BP1", "BP2", "TYPE") %in% colnames(x = cnv)
  if (any(!tst)) stop("cnv does not appear to contain the correct data")
  cnv <- cnv[,c("ID", "CHR", "BP1", "BP2", "TYPE")]

  # remove duplicate records
  dups <- duplicated(x = cnv)
  if (verbose && {sum(dups) > 0}) {
    cat("\t", sum(dups), "duplicate cnv records removed\n")
  }
  cnv <- cnv[!dups,]

  cnv <- cnv[order(cnv$ID,cnv$CHR,cnv$BP1),]

  tst <- .dataCheck(cnv = cnv)

  kernel <- .CAUCkernel(cnv = cnv,
                        nCore = nCore,
                        outFileKernel = outFileKernel,
                        verbose = verbose)

  rm(cnv)

  return( .pv_func(pheno = pheno,
                   phenoY = phenoY,
                   X = X,
                   phenoType = phenoType,
                   kernel = kernel,
                   verbose = verbose) )

}   
