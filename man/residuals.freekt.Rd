\name{residuals.freekt}

\Rdversion{1.1}

\alias{residuals.freekt}

\title{
Compute Residuals For Free-Knot Spline
}

\description{
This function computes the residuals, given the optimal values for lambda and the locations of the knots.
}

\usage{
\method{residuals}{freekt}(object, \dots)
}

\arguments{
  \item{object}{
  An object of class "\code{freekt}" obtained by using one of the fitting algorithms.
}
  \item{\dots}{
  Additional arguments to be passed to residuals.freekt.  Currently ignored.
}
}

\value{
A vector containing the residuals.
}

\author{
Steven Spiriti
}

\seealso{
\code{\link{fitted.freekt}} to compute the fitted values.
}

\examples{
x <- 0:30/30
truey <- x*sin(10*x)
set.seed(10556)
y <- truey + rnorm(31, 0, 0.2)
xy.freekt <- freelsgen(x, y, degree = 2, numknot = 2, 555)
plot(x, residuals(xy.freekt))
}

\keyword{
nonparametric
regression
smooth
}

