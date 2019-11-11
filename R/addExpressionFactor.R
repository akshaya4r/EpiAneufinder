addExpressionFactor <- function(bins, gene.annotation=NULL) {
  UseMethod("addExpressionFactor")
}

addExpressionFactor.GRanges <- function(bins, gene.annotation=NULL) {
  txdb <- getFromNamespace(gene.annotation, ns=gene.annotation)
  seqlevelsStyle(txdb) <- seqlevelsStyle(bins)[1]
  trs <- reduce(transcripts(txdb))
  trs$width <- width(trs)
  trs <- mergeByOverlaps(bins,trs)
  trs <- aggregate(trs$width, by=list(as.character(trs$bins)), FUN=sum)
  colnames(trs) <- c('gr','totalTrWidth')
  trs$seqnames <- sapply(strsplit(trs$gr,":"), `[`, 1)
  trs$start <- as.numeric(sapply(strsplit(sapply(strsplit(trs$gr,":"), `[`, 2),"-"), `[`, 1))
  trs$end <- as.numeric(sapply(strsplit(trs$gr,"-"), `[`, 2))
  trs$gr <- NULL
  trs.gr <- sort(makeGRangesFromDataFrame(trs, keep.extra.columns = TRUE))
  nontrs.gr <- subsetByOverlaps(bins,trs.gr,invert=TRUE)
  nontrs.gr$GC.content <- NULL
  nontrs.gr$totalTrWidth <- 0.3
  gr <- sort(c(trs.gr,nontrs.gr))
  bins$expressionFactor <- gr$totalTrWidth
  bins
}

addExpressionFactor.list <- function(bins, gene.annotation=NULL) {
  lapply(bins, addExpressionFactor, txdb)
}

addExpressionFactor.default <- function(bins, gene.annotation=NULL) {
  stop("Do not know how to get expression factor for type ", sQuote(class(bins)))
}
