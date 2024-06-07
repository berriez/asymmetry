\name{asymmetry-package}
\alias{asymmetry-package}
\alias{asymmetry}
\docType{package}
\title{
\packageTitle{asymmetry}
}
\description{
\packageDescription{asymmetry}
}

\author{
\packageAuthor{asymmetry}

Maintainer: \packageMaintainer{asymmetry}
}
\details{
Asymmetry in general, and in proximity relations in particular means that the relation
from $i$ to $j$  is not equal to the relation in the opposite direction, that is, from $j$ to $i$.
This package offers functions for the analysis of asymmetry. For instance, to obtain a heatmap of the skew symmetric part of the data, we use the hmap function. Other functions available are the slidevector model, the asymscal model, a multidimensional scaling model with unique dimensions and, of course. A cornerstone of this package is the decomposition of an asymmetric matrix into a symmetric part and a skew symmetric part. This is a well-known mathematical decomposition and is used extensively in this package.

The analysis of asymmetry was developed as an extension to a symmetric method such as multidimensional scaling. We start with the definition of an asymmetric matrix. An asymmetric data matrix has the same number of rows and columns, and these rows and columns refer to the same set of objects. At least some elements in the upper-triangle are different from the corresponding elements in the lower triangle for a matrix to be asymmetric.

Usually this decomposition is applied to data to study the two components separately. But it can also be applied to model parameters. Here, we use a decomposition to residuals , that is, to the deviations from a model to the data.

}
\references{
Zielman, B., and Heiser, W. J. (1993), The analysis of asymmetry by a slide-vector, Psychometrika, 58, 101-114.
}


