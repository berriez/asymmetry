hmap <- function(x, dominance = FALSE, ...)
{
  dendrogram <- NULL
  Rowv <- NULL
  Colv <- NULL
  trace <- NULL
  if (is.matrix(x))
    x <- .5 * (x - t(x))
  if (is(x, "skewsymmetry"))
    x <- x$A

  if (nrow(x) != ncol(x))
    stop("the same number of rows and columns are expected")
  bin <- x * 0
  bin[x > 0] <- 1
  bin[x < 0] <- -1
  rsbin <- rowSums(bin)
  x <- x[order(rsbin), order(rsbin)]
  #  keuze tussen binair en original
  if (dominance)
  {
    bin <- bin[order(rsbin), order(rsbin)]
    heatmap.2(
      bin,
      dendrogram = 'none',
      Rowv = FALSE,
      Colv = FALSE,
      trace = 'none',
      ...
    )
  }
  else {
    heatmap.2(
      x,
      dendrogram = 'none',
      Rowv = FALSE,
      Colv = FALSE,
      trace = 'none',
      ...
    )
  }
}
