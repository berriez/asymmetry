\name{plot.skewsymmetry}
\alias{plot.skewsymmetry}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting method for the skew-symmetric part of an asymmetric matrix
}
\description{
This plotting method provides a multidimensional representation of skew-symmetry based on the singular value decomposition (SVD). The properties of the SVD of a skew-symmetric matrix were given by Gower (1977) where also the guidelines for the interpretation of diagrams obtained by plotting pairs of singular vectors is described. The singular vectors of a skew-symmetric matrix come in pairs with equal singular values. The diagrams are not interpreted by comparing distances between point as is usual in multidimensional scaling, but by comparing areas formed by two points and the origin. The singular vectors span a plane, and the area of the triangle between two points and the origin represents skew-symmetry.
The sign of the skew-symmetry between two points is modelled by a direction in the plane. Going clockwise the area between two points and the origin is negative, goint counter clockwise the area is positive.
}
\usage{
\method{plot}{skewsymmetry}(x, plot.plane = 1, yplus = 0, xlab, ylab, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class skewsymmetry
}
  \item{plot.plane}{
Integer indicating which plane to plot
}
  \item{yplus}{
Offset for the labels above the object points
}
  \item{xlab}{
Label for the x-axis
}
  \item{ylab}{
Label for the y-axis
}
  \item{\dots}{
Further plot arguments
}
}
\references{
Gower, J.C. (1977) The analysis of asymmetry and orthogonality. In: Recent Developments in Statistics ( J. Barra, F. Brodeau, G. Romier & B. van Cutsem, Eds.), 109-123. North Holland, Amsterdam.
}
\keyword{Methods}

