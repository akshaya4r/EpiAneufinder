#' Correction for the confounding effect of expression
#'
#' Correct a list of \code{\link{binned.data}} by gene expression.
#'
#' Two methods are available for GC correction:  Option
#' \code{method='quadratic'} uses the method described in the Supplementary of
#' \code{citation("AneuFinder")}. Option \code{method='loess'} uses a loess fit
#' to adjust the read count.
#' 
#' @param binned.data.list A \code{list} with \code{\link{binned.data}} objects
#'   or a list of filenames containing such objects.
#' @param GC.BSgenome A \code{BSgenome} object which contains the DNA sequence
#'   that is used for the GC correction.
#' @param same.binsize If \code{TRUE} the GC content will only be calculated
#'   once. Set this to \code{TRUE} if all \code{\link{binned.data}} objects
#'   describe the same genome at the same binsize and stepsize.
#' @param method One of \code{c('quadratic', 'loess')}. Option
#'   \code{method='quadratic'} uses the method described in the Supplementary of
#'   \code{citation("AneuFinder")}. Option \code{method='loess'} uses a loess fit
#'   to adjust the read count.
#' @param return.plot Set to \code{TRUE} if plots should be returned for visual
#'   assessment of the GC correction.
#' @param bins A \code{\link{binned.data}} object with meta-data column 'GC'.
#'   If this is specified, \code{GC.BSgenome} is ignored. Beware, no format
#'   checking is done.
#' @return A \code{list()} with \code{\link{binned.data}} objects with adjusted
#'   read counts. Alternatively a \code{list()} with
#'   \code{\link[ggplot2]{ggplot}} objects if \code{return.plot=TRUE}.
#' @author Aaron Taudt
#' @importFrom Biostrings Views alphabetFrequency
#' @importFrom stats lm predict loess
#' @importFrom reshape2 melt
#' @export
#' @examples
#' ## Get a BED file, bin it and run GC correction
#' bedfile <- system.file("extdata", "KK150311_VI_07.bam.bed.gz", package="AneuFinderData")
#' binned <- binReads(bedfile, assembly='mm10', binsize=1e6,
#'                    chromosomes=c(1:19,'X','Y'))
#' plot(binned[[1]], type=1)
#' if (require(BSgenome.Mmusculus.UCSC.mm10)) {
#'   binned.GC <- correctGC(list(binned[[1]]), GC.BSgenome=BSgenome.Mmusculus.UCSC.mm10)
#'   plot(binned.GC[[1]], type=1)
#' }
correctExpression <- function(reads) {
  print("Correcting for expression")
  UseMethod("correctExpression")
}

# check for same genome if passed a list?
# should this addGCcontent if not already present?
correctExpression.GRanges <- function(reads) {
  reads$uncorrected <- reads$counts
  reads$counts <- (reads$uncorrected/reads$expressionFactor)*width(reads[1])
  reads
}

correctExpression.list <- function(reads) {
  lapply(reads, correctExpression.GRanges)
}

correctGC.default <- function(reads) {
  stop("Do not know how to correct expression for type ", sQuote(class(reads)))
}
