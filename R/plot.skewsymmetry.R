plot.skewsymmetry <- function(x, plot.plane = 1, yplus = 0,
           xlab, ylab, ...)
  {
    #
    # add defaults
    #
    y1 <- plot.plane * 2
    x1 <- y1 - 1
    if (missing(xlab))
      xlab <- paste("Dimension", x1, sep = " ")
    else
      xlab <- xlab
    if (missing(ylab))
      ylab <- paste("Dimension", y1, sep = " ")
    else
      ylab <- ylab
    # max1 <- max(x$sv[,x1])
    # max2 <-
    max <- max(max(x$sv[, x1]), max(x$sv[, y1]))
    min <- min(min(x$sv[, x1]), min(x$sv[, y1]))


    plot(
      x$sv[, x1],
      x$sv[, y1],
      xlab = xlab,
      ylab = ylab,
      xlim = c(min, max),
      ylim = c(min, max),
      ...
    )
    if (!is.null(rownames(x$sv))) {
      text(x$sv[, x1], yplus + x$sv[, y1], rownames(x$sv))
    } else {
      text(x$sv[, x1], yplus + x$sv[, y1], c(1:x$nobj))
    }
    text(0, 0, "O")
  }
