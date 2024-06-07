\name{hmap}
\alias{hmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Heatmap for skew-symmetric data
}
\description{
This heatmap displays the values of a skew-symmetric matrix by colors. The option \code{dominance} orders the rows and columns of the matrix in such a way that the values in the uppertriangle are positive and the values in the lower triangle are negative. The order is calculated from the row-sums of the signs obtained from the skew-symmetric matrix.
}
\usage{
hmap(x, dominance = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A square matrix, either skew-symmetric or asymmetric, or an object of class \code{decomposition}. If an asymmetric matrix is given, the skew-symmetric part is computed.
}
  \item{dominance}{
If true the signs of the skew-symmetric matrix are shown in the heatmap, if set to false the values in this matrix are shown.
}
  \item{\dots}{Further plot arguments: see \code{\link[gplots]{heatmap.2}} for detailed information.}
}
\examples{
data(studentmigration)
hmap(studentmigration, dominance = TRUE, col = c("red", "white", "blue"))
}
\keyword{Methods}

