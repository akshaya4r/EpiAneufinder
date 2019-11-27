calcMeanGeneExpression <- function(reads, gene.annotation=NULL) {
  UseMethod("addExpressionFactor")
}

calcMeanGeneExpression.GRanges <- function(reads, gene.annotation=NULL) {
  txdb <- getFromNamespace(gene.annotation, ns=gene.annotation)
  genes <- sort(keepStandardChromosomes(genes(txdb)))
  seqlevelsStyle(genes) <- seqlevelsStyle(reads)[1]
  subsetByOverlaps(reads,genes)
}

calcMeanGeneExpression.list <- function(reads, gene.annotation=NULL) {
  lapply(reads, addExpressionFactor, txdb)
}

calcMeanGeneExpression.default <- function(reads, gene.annotation=NULL) {
  stop("Do not know how to get expression factor for type ", sQuote(class(reads)))
}