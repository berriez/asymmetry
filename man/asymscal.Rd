\name{asymscal}
\alias{asymscal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weighted Euclidean Model for Asymmetric Matrices
}
\description{
This function fits a weighted multidimensional scaling model that is known as the asymscal model. This model is an extension of the symmetric Euclidean distance model proposed by Young (1975). The model is fitted in a stress majorization framework called SMACOF, whereas Young fitted this model using a least squares algorithm. Asymmetry is modelled by differential weighting of the dimensions of a multidimensional scaling configuration. When a subject compares object \emph{i} to \emph{j} he or she may use different weights when comparing object \emph{j} to \emph{i} In addition to these weights, the locations of the objects are jointly estimated from the data.
\deqn{d_{ij}(X)=\sqrt{\sum_{s=1}^pv_{is}(x_{is}-x_{js})^2}.}

}
\usage{
asymscal(data, ndim = 2, start = NULL, verbose = FALSE, itmax = 10000, eps = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Asymmetric dissimilarity matrix
}
  \item{ndim}{
Number of dimensions
}
  \item{start}{
Optional configuration with starting values, the default is a random start configuration
}
  \item{verbose}{
If TRUE, stress values during the iterations are printed
}

  \item{itmax}{
Maximum number of iterations
}
  \item{eps}{
Convergence criterion for Stress
}
}
\details{
This function exploits a connection between the INDSCAL model and the asymscal model. This method inherits the methods for plotting an printing from the \code{smacofIndDiff} in the smacof package. Basically, the asymscal takes two steps. First, this function sets up the appropriate dissimilarity and missing data structure for a three-way multidimensional scaling model, then a call to the method \code{smacofIndDiff} in the imported package smacof is made. After correcting for the normalization applied to the data by \code{smacofIndDiff}, the results can be displayed and plotted by the methods in the package \code{smacof}.
The original algorithm for fitting the asymscal model fits squared distances. This function is based on majorization, and fits distances and not squared distances. The configuration matrix is normalized, the sum of squares of the columns of this matrix are equal to one.
}
\value{
  \item{delta}{Observed dissimilarities}
  \item{obsdiss}{List of observed dissimilarities, normalized}
  \item{gspace}{Joint configurations aka group stimulus space}
  \item{cweights}{Configuration weights}
  \item{stress}{Stress-1 value}
  \item{resmat}{Matrix with residuals}
  \item{rss}{Residual sum-of-squares}
  \item{spp}{Stress per point}
  \item{ndim}{Number of dimensions}
  \item{model}{Type of the asymmetric scaling model}
  \item{niter}{Number of iterations}
  \item{nobj}{Number of objects}
}
\references{
Young, F. W. (1975). An asymmetric Euclidean model for multi-process asymmetric data. Paper presented at the U.S.-Japan Seminar on Multidimensional scaling, San Diego, U.S.A.
}
\examples{
\dontrun{
data("asymscalexample")
t<-asymscal(asymscalexample, ndim = 2, itmax = 10000, eps = 1e-10)
t$cweights
round(t$cweights, 3)
plot(t, plot.type = "confplot")
plot(t, plot.type = "bubbleplot")
plot(t, plot.type = "stressplot")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MDS}
