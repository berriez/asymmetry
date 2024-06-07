\name{summary.skewsymmetry}
\alias{summary.skewsymmetry}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary method of the decomposition
}
\description{
Prints a decomposition of the sum of squares of an asymmetric matrix. The first column gives the sum of squares, and the second column gives the percentages of the two components. This decomposition can be applied to data, but also to a matrix of residuals obtained from a fitted model.}
\usage{
\method{summary}{skewsymmetry}(object, ...)
}
\arguments{
  \item{object}{
An object of class \code{decomposition}
}
  \item{\dots}{
Further parameters
}
}
\examples{
data(Englishtowns)
q <- skewsymmetry(Englishtowns)
summary(q)
}
\keyword{Methods}
