summary.skewsymmetry <- function(object, ...)
{
  tssq <- sum(as.vector(object$S + object$A) ^ 2)
  tab <-
    rbind(c(sum(as.vector(object$S) ^ 2), 100 * sum(as.vector(object$S) ^ 2) /
              tssq),
          c(sum(as.vector(object$A) ^ 2), 100 * sum(as.vector(object$A) ^
                                                      2) / tssq),
          c(tssq, 100 * tssq / tssq))
  colnames(tab) <- c("SSQ", "Percent")
  rownames(tab) <- c("Symmetry", "Skew-symmetry", "Total")
  tab
}
