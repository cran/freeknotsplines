\name{fitted.freekt}

\Rdversion{1.1}

\alias{fitted.freekt}

\title{
Compute Fitted Values For Free-Knot Spline
}

\description{
This function computes the fitted values, given the amount of the penalty (if applicable) and the locations of the knots.
}

\usage{
\method{fitted}{freekt}(object, xfit = object@x, \dots)
}

\arguments{
  \item{object}{
  An object of class "\code{freekt}" obtained by using one of the fitting algorithms.
}
  
  \item{xfit}{
  A vector of x values at which to compute the fitted values.  Defaults to the x values of the data.
}

  \item{\dots}{
  Additional arguments to be passed to fitted.freekt.  Currently ignored.
}
}

\value{
A vector containing the fitted values.
}

\author{
Steven Spiriti
}

\seealso{
\code{\link{residuals.freekt}} for the residuals.
}

\examples{
x <- 0:30/30
truey <- x*sin(10*x)
set.seed(10556)
y <- truey + rnorm(31, 0, 0.2)
xy.freekt <- freelsgen(x, y, degree = 2, numknot = 2, 555)
fitted(xy.freekt)
}

\keyword{
nonparametric
regression
smooth
}

