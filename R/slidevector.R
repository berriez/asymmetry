slidevector <- function(data, weight = NULL, ndim = 2, verbose = FALSE,
                        itmax = 125, eps = 1e-12, rotate = TRUE)
{
    if (sum(data < 0) > 0)
      stop("data for the slide-vector model should be non-negative")
    if (nrow(data) != ncol(data))
      stop("the same number of rows and columns are expected in the data matrix")
    nobj <- nrow(data)
    if (is.null(weight))
      weight <- 0 * data + 1
    if (sum(weight < 0) > 0)
      stop("weight for the slide-vector model should be positive")
    if (nrow(weight) != ncol(weight))
      stop("the same number of rows and columns are expected in the weight matrix")
    we <- rbind(cbind(weight * 0, weight), cbind(t(weight),  weight * 0))
    v <- diag(rowSums(we)) - we
    da <- rbind(cbind(data * 0, data), cbind(t(data), data * 0))
    ztil <- matrix(runif(nobj * 2 * ndim), 2 * nobj, ndim)
    dist <- as.matrix(dist(ztil))
    ind <- dist < .000000001
    dist[ind] <- 1
    BMAT <- we * da / dist
    BMAT[ind] <- 0
    BMAT <- diag(rowSums(BMAT)) - BMAT #VX
    zsym <- rbind(diag(nobj), diag(nobj))
    z <- cbind(zsym, c(rep(1, nobj), rep(0, nobj)))
    zvz <- t(z) %*% v %*% z
    nsp <- c(rep(1, nobj), 0) %*% t(c(rep(1, nobj), 0)) / nobj
    zvz_inv <- solve(zvz + nsp) - nsp
    zvz_inv_s <-
      solve(zvz[1:nobj, 1:nobj] + nsp[1:nobj, 1:nobj]) - nsp[1:nobj, 1:nobj]
    #ztil <- z%*%zvz_inv%*%t(z)%*%ztil #hier ook startoplossing toevoegen
    if (rotate) {
      ztil[, 1] <- z %*% zvz_inv %*% t(z) %*% ztil[, 1]
      ztil[, 2:ndim] <- zsym %*% zvz_inv_s %*% t(zsym) %*% ztil[, 2:ndim]
    } else {
      ztil <- z %*% zvz_inv %*% t(z) %*% ztil
    }
    # model satisfies constraints
    dist <- as.matrix(dist(ztil))
    stressold <- sum(weight * (data - dist[1:nobj, (nobj + 1):(2 * nobj)]) ^
                       2)
    #print(stressold)
    for (i in 1:itmax) {
      ind <- dist < .000000001
      dist[ind] <- 1
      BMAT <- we * da / dist
      BMAT[ind] <- 0
      BMAT <- diag(rowSums(BMAT)) - BMAT #VX
      zup <- BMAT %*% ztil
      if (rotate) {
        ztil[, 1] <- z %*% zvz_inv %*% t(z) %*% zup[, 1]
        ztil[, 2:ndim] <- zsym %*% zvz_inv_s %*% t(zsym) %*% zup[, 2:ndim]

      } else {
        ztil <- z %*% zvz_inv %*% t(z) %*% zup
      }
      dist <- as.matrix(dist(ztil))
      stress <- sum(weight * (data - dist[1:nobj, (nobj + 1):(2 * nobj)]) ^
                      2)
      if (verbose)
        print(stress)
      fu <- stressold - stress
      if (stressold - stress < eps)
        break
      stressold <- stress
    }
    niter <- i
    ## if loop
    ## second
    if (rotate) {
      X <- matrix(0, nobj + 1, 2)
      X[, 1] <- zvz_inv %*% t(z) %*% zup[, 1]
      X[1:nobj, 2:ndim] <- zvz_inv_s %*% t(zsym) %*% zup[, 2:ndim]
    } else {
      X <- zvz_inv %*% t(z) %*% zup #moet v tussen
    }
    ## end
    mat <- X[1:nobj,]
    Z <- X[nobj + 1,]
    if (!is.null(rownames(data)))
      rownames(mat) <- rownames(data)
    colnames(mat) <- paste("D", 1:(dim(mat)[2]), sep = "")
    pred <- dist[1:nobj, (nobj + 1):(2 * nobj)]
    resid <- data - pred
    result <-
      list(
        ndim = ndim,
        stress = stress,
        confi = mat,
        slvec = Z,
        resid = resid,
        niter = niter,
        nobj = nobj,
        pred = pred,
        model = "Slide-vector model",
        ztil = ztil
      )
    class(result) <- "slidevector"
    return(result)
}
print.slidevector <-
  function(x, ...)
  {
    cat("Dimensions:              ")
    cat(x$ndim)
    cat("\n")
    cat("Number of objects:       ")
    cat(x$nobs)
    cat("\n")
    cat("Number of iterations:    ")
    cat(x$niter)
    cat("\n")
    cat("Stress:                   ")
    cat(x$stress)
    cat("\n")
  }
summary.slidevector <-
  function(object, ...)
  {
    cat("\n")
    cat("Configurations:\n")
    print(round(object$confi, 4))
    cat("\n")
    cat("Slide-vector:\n")
    cat(round(object$slvec, 4))
  }
