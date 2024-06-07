skewsymmetry <- function(x)
{
  if (nrow(x) != ncol(x))
    stop("the same number of rows and columns are expected")
  S <- .5 * (t(x) + x)
  A <- .5 * (x - t(x))
  lin <- rowMeans(A)
  nobj <- nrow(x)
  Svectors <- matrix(0.0, nrow(x), nrow(x))
  s <- svd(A)
  for (i in seq(1, (nrow(x) - 1), by = 2)) {
    pred1 <-
      s$u[, i:(i + 1)] %*% (s$d[i:(i + 1)] * diag(2)) %*% t(s$v[, i:(i + 1)]) #good
    pred2 <-
      (s$u[, i] %*% t(s$u[, (i + 1)]) - s$u[, (i + 1)] %*% t(s$u[, i])) * s$d[i]
    if (sum(abs(pred1 - pred2) > 1e-10) == 0)
    {
      Svectors[, i] <- s$u[, i]
      Svectors[, i + 1] <- s$u[, i + 1]
    } else {
      Svectors[, i + 1] <- s$u[, i]
      Svectors[, i] <- s$u[, i + 1]
    }
  }
  if (!is.null(rownames(x)))
    rownames(Svectors) <- rownames(x)

  lk <-
    list(
      S = S,
      A = A,
      linear = lin,
      sv = Svectors,
      sval = s$d,
      nobj = nobj
    )
  class(lk) <- "skewsymmetry"
  lk
}
