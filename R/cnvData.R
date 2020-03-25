#' Pseudo Copy Number Variants Data
#' 
#' This data set includes simulated CNV data in PLINK CNV data format.
#' The data are also available from the authors through the url 
#' provided below. These data were generated following the simulation 
#' study used to illustrate the method in the original manuscript also 
#' referenced below; it has been reduced to include only 600 individuals. 
#' These data are not meaningful and are intended for demonstration purposes only.
#' 
#' @usage data(cnvData)
#'
#' @format cnvData is a data.frame containing 522 observations with 5 columns: 
#' \describe{
#' \item{ID}{character patient identifier.}
#' \item{CHR}{CNV chromosome.}
#' \item{BP1}{starting location in base pairs.} 
#' \item{BP2}{ending location in base pairs.} 
#' \item{TYPE}{copy number (0,1,2,3,or 4).}
#' }
#'
#' @references Brucker, A., Lu, W., Marceau West, R., Yu, Q-Y., Hsiao, C. K.,
#'   Hsiao, T-H., Lin, C-H., Magnusson, P. K. E., Holloway, S. T., 
#'   Sullivan, P. F., Szatkiewicz, J. P., Lu, T-P., and
#'   Tzeng, J-Y. Association testing using Copy Number Profile Curves (CONCUR)
#'   enhances power in copy number variant analysis. <doi:10.1101/666875>.
#'
#' @references \url{https://www4.stat.ncsu.edu/~jytzeng/Software/CONCUR/}
#' @keywords datasets
"cnvData"

#' Pseudo Covariate Data
#' 
#' This data set includes simulated covariate data.
#' These data were generated as draws from a Binom(1,0.5) distribution for the
#' 800 individuals in the example data provided with the package.
#' These data are not meaningful and are intended for demonstration purposes only.
#' 
#' @usage data(cnvData)
#'
#' @format covData is a data.frame containing 400 observations with 2 columns
#' \describe{
#' \item{ID}{character patient identifier.}
#' \item{SEX}{binary indicator of M/F.}
#' }
#'
#' @references Brucker, A., Lu, W., Marceau West, R., Yu, Q-Y., Hsiao, C. K.,
#'   Hsiao, T-H., Lin, C-H., Magnusson, P. K. E., Holloway, S. T., 
#'   Sullivan, P. F., Szatkiewicz, J. P., Lu, T-P., and
#'   Tzeng, J-Y. Association testing using Copy Number Profile Curves (CONCUR)
#'   enhances power in copy number variant analysis. <doi:10.1101/666875>.
#'
#' @keywords datasets
"covData"

#' Pseudo Phenotype Data
#' 
#' This data set includes simulated phenotype data.
#' These data include a binary phenotype and a normally distributed continuous
#' phenotype that are randomly generated independent of the CNV data.
#' These data are not meaningful and are intended for demonstration purposes only.
#' 
#' @usage data(cnvData)
#'
#' @format phenoData is a data.frame containing 400 observations with 3 columns
#' \describe{
#' \item{ID}{character patient identifier.}
#' \item{PHEB}{binary phenotype.}
#' \item{PHEC}{continuous phenotype.}
#' }
#'
#' @references Brucker, A., Lu, W., Marceau West, R., Yu, Q-Y., Hsiao, C. K.,
#'   Hsiao, T-H., Lin, C-H., Magnusson, P. K. E., Holloway, S. T., 
#'   Sullivan, P. F., Szatkiewicz, J. P., Lu, T-P., and
#'   Tzeng, J-Y. Association testing using Copy Number Profile Curves (CONCUR)
#'   enhances power in copy number variant analysis. <doi:10.1101/666875>.
#'
#' @keywords datasets
"phenoData"
