mdsunique <- function(data, weight = NULL, ndim = 2, verbose = FALSE,
           itmax = 125, eps = 1e-12)
  {
    if (sum(data < 0) > 0)
      stop("data for this MDS model should be positive")
    if (nrow(data) != ncol(data))
      stop("the same number of rows and columns are expected")
    nobj <- nrow(data)
    if (is.null(weight))
      weight <- 0 * data + 1
    if (sum(weight < 0) > 0)
      stop("weights for this model should be positive")
    if (nrow(weight) != ncol(weight))
      stop("the same number of rows and columns are expected in the weight matrix")

    we <-  rbind(cbind(weight * 0, weight), cbind(t(weight), weight * 0))
    v <- diag(rowSums(we)) - we
    da <-  rbind(cbind(data * 0, data), cbind(t(data), data * 0))
    ztil <- matrix(runif(nobj * 2 * ndim), 2 * nobj, ndim)
    dist <- as.matrix(dist(ztil))
    ind <- dist < .000000001
    dist[ind] <- 1
    BMAT <- we * da / dist
    BMAT[ind] <- 0
    BMAT <- diag(rowSums(BMAT)) - BMAT #VX
    z <- rbind(diag(nobj), diag(nobj))
    zvz <-  t(z) %*% v %*% z
    nsp <- c(rep(1, nobj)) %*% t(c(rep(1, nobj))) / nobj
    fdim <- ndim + 2 * nobj
    ztil <- cbind(ztil, diag(2 * nobj))
    zvz_inv <- solve(zvz + nsp) - nsp
    dist <- as.matrix(dist(ztil))
    stressold <- sum(weight * (data - dist[1:nobj, (nobj + 1):(2 * nobj)]) ^
                       2)
    for (i in 1:itmax) {
      ind <- dist < .000000001
      dist[ind] <- 1
      BMAT <- we * da / dist
      BMAT[ind] <- 0
      BMAT <- diag(rowSums(BMAT)) - BMAT #VX
      ztil <- BMAT %*% ztil
      ztil[, 1:ndim] <- z %*% zvz_inv %*% t(z) %*% ztil[, 1:ndim]
      uniq <- diag(ztil[, (ndim + 1):fdim]) / diag(v)
      ztil <- cbind(ztil[, 1:ndim], diag(uniq))
      dist <- as.matrix(dist(ztil))
      if (verbose)
        print(sum(weight * (data - dist[1:nobj, (nobj + 1):(2 * nobj)]) ^
                    2))
      stress <- sum(weight * (data - dist[1:nobj, (nobj + 1):(2 * nobj)]) ^
                      2)
      if (stressold - stress < eps)
        break
      stressold <- stress
    }
    mat <- ztil
    colnames(mat) <- paste("D", 1:(dim(mat)[2]), sep = "")

    pred <- as.matrix(dist(ztil))

    pred <- pred[1:nobj, (nobj + 1):(2 * nobj)] # for the todo list
    pred <- dist[1:nobj, (nobj + 1):(2 * nobj)]
    fulldim <- ndim + 2 * nobj
    resid <- data - pred
    confi <- mat[1:nobj, 1:ndim]
    uni <- diag(ztil[, (ndim + 1):ncol(ztil)])
    row <- uni[1:nobj]
    col <- uni[(nobj + 1):(2 * nobj)]
    unique <- cbind(row, col)
    result <-
      list(
        ndim = ndim,
        stress = stress,
        X = ztil,
        confi = confi,
        resmat = resid,
        fulldim = fulldim,
        niter = i,
        nobj = nobj,
        model = "MDS with unique dimensions",
        unique = unique,
        row = row,
        col = col
      )
    class(result) <- "mdsunique"
    return(result)
  }

plot.mdsunique <-
  function(x,
           plot.dim = c(1, 2),
           yplus = 0,
           xlab,
           ylab,
           ...) {
    #
    # add defaults
    #
    x1 <- plot.dim[1]
    y1 <- plot.dim[2]
    if (missing(xlab))
      xlab <- paste("Dimension", x1, sep = " ")
    else
      xlab <- xlab
    if (missing(ylab))
      ylab <- paste("Dimension", y1, sep = " ")
    else
      ylab <- ylab
    plot(x$confi[, x1], x$confi[, y1], xlab = xlab, ylab = ylab, ...)
    if (!is.null(rownames(x$confi))) {
      text(x$confi[, x1], yplus + x$confi[, y1], rownames(x$confi))
    } else {
      text(x$confi[, x1], yplus + x$confi[, y1], c(1:x$nobs))
    }
  }

print.mdsunique <-
  function(x, ...)
  {
    cat("Dimensions:              ")
    cat(x$ndim)
    cat("\n")
    cat("Number of objects:       ")
    cat(x$nobj)
    cat("\n")
    cat("Number of iterations:    ")
    cat(x$niter)
    cat("\n")
    cat("Stress:                   ")
    cat(x$stress)
    cat("\n")
  }
summary.mdsunique <-
  function(object, ...)
  {
    cat("\n")
    cat("Configurations:\n")
    xda <- cbind(object$confi, object$row, object$col)
    nam <- colnames(object$confi)
    names <- c(nam, "Row", "Col")
    colnames(xda) <- names
    print(round(xda, 4))
  }
