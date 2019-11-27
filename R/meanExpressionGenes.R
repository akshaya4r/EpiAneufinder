calcMeanGeneExpression <- function(reads, gene.annotation=NULL) {
  UseMethod("calcMeanGeneExpression")
}

calcMeanGeneExpression.GRanges <- function(reads, gene.annotation=NULL) {
  #txdb <- getFromNamespace(gene.annotation, ns=gene.annotation)
  # genes <- sort(keepStandardChromosomes(genes(txdb), pruning.mode = 'coarse'))
  genes <- gene.annotation
  seqlevelsStyle(genes) <- seqlevelsStyle(reads)[1]
  genes$readcount <- countOverlaps(genes, reads, type = 'any')
  genes <- genes[which(genes$readcount>0 & genes$readcount<quantile(genes$readcount, 0.997))]
  genes <- genes[which(genes$readcount>mean(genes$readcount))]
  reads <- subsetByOverlaps(reads, genes)
  print("Keeping only genes with above mean expression")
  reads
}

calcMeanGeneExpression.list <- function(reads, transcript.db=NULL) {
  lapply(reads, addExpressionFactor, txdb)
}

calcMeanGeneExpression.default <- function(reads, transcript.db=NULL) {
  stop("Do not know how to get expression factor for type ", sQuote(class(reads)))
}
