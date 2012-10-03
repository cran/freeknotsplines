\name{coef.freekt}

\Rdversion{1.1}

\alias{coef.freekt}

\title{
Compute Coefficients of B-Splines For Free-Knot Splines
}

\description{
This function computes the coefficients of the B-splines for free-knot splines, given the amount of the penalty (if applicable) and the locations of the knots.
}

\usage{
\method{coef}{freekt}(object, \dots)
}

\arguments{
  \item{object}{
  An object of class "\code{freekt}" obtained by using one of the fitting algorithms.
}
  \item{\dots}{
  Additional arguments to be passed to coef.freekt.  Currently ignored.
}
}

\value{
A vector containing the coefficients of the B-splines.
}

\author{
Steven Spiriti
}

\seealso{
\code{\link{fitted.freekt}} to compute the fitted values and \code{\link{residuals.freekt}} to compute the residuals.
}

\examples{
x <- 0:30/30
truey <- x*sin(10*x)
set.seed(10556)
y <- truey + rnorm(31, 0, 0.2)
xy.freekt <- freelsgen(x, y, degree = 2, numknot = 2, 555)
coef(xy.freekt)
}

\keyword{
nonparametric
regression
smooth
}