\name{mdsunique}
\alias{mdsunique}
\alias{print.mdsunique}
\alias{summary.mdsunique}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MDS Model with Unique Dimensions
}
\description{
This asymmetric MDS model proposed by Holman (1979) and is related to a constrained scaling model developed by Bentler & Weeks (1982).The model has two sets of dimensions, shared or common dimensions and the other set are unique dimensions. There are common dimensions that apply to all objects in the analysis, and unique dimensions that apply to one object and not the other objects. A unique dimension has a non zero value for only one object, the coordinates for the other objects are zero. There are as many unique dimensions as there are objects. An asymmetric version of this model has two sets of unique dimensions: one for the rows and one for the columns. The distance in this model is defined as:
\deqn{d_{ij}(X)=\sqrt{\sum_{s=1}^p (x_{is}-x_{js})^2 + r_{i}^{2}+c_{j}^{2}}.}
}
\usage{
mdsunique(data, weight = NULL, ndim = 2, verbose = FALSE, itmax = 125, eps = 1e-12)

}
%- maybe also 'usage' for other objects documented here.



\arguments{
  \item{data}{ Asymmetric dissimilarity matrix

%%     ~~Describe \code{data} here~~
}
  \item{weight}{Optional non-negative matrix with weights, if no weights are given all weights are set equal to one}
  \item{ndim}{Number of dimensions
%%     ~~Describe \code{dim} here~~
}
  \item{verbose}{ If true, prints the iteration history to screen
%%     ~~Describe \code{verbose} here~~
}
  \item{itmax}{Maximum number of iterations
%%     ~~Describe \code{itmax} here~~
}
  \item{eps}{  Convergence criterion for Stress
%%     ~~Describe \code{eps} here~~
}
}
\value{
  \item{ndim }{Number of dimensions of the configuration}
  \item{fulldim }{Number of dimensions of the full model, this equals \code{ndim} + the number of rows + the number of columns}
  \item{stress }{The raw stress for this model}
  \item{confi}{Returns the configuration matrix of shared dimensions of this multidimensional scaling model}
  \item{X}{Returns the configuration matrix of the full model consisting of shared and unique dimensions}
  \item{niter}{The number of iterations for the algorithm to converge}
  \item{nobj}{The number of objects in this model}
  \item{resid}{A matrix with raw residuals}
  \item{model}{Name of this asymmetric multidimensional scaling model}
  \item{row}{The unique dimensions for the rows}
  \item{col}{The unique dimensions for the columns}
  \item{unique}{The unique dimensions}
}
\examples{
\dontrun{
data("studentmigration")
mm<-studentmigration
mm[mm==0]<-.5          # replace zeroes by a small number
mm <- -log(mm/sum(mm)) # convert similarities to dissimilarities
v<-mdsunique(mm, ndim = 2, itmax = 2100, verbose=FALSE, eps = .0000000001)
plot(v, yplus = .3)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MDS}% use one of  RShowDoc("KEYWORDS")
