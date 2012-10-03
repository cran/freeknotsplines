\name{summary.freekt}

\Rdversion{1.1}

\alias{summary.freekt}

\title{
Summarize Free-Knot Spline Fit
}

\description{
This function displays a summary of the fit obtained using a free-knot spline.
}

\usage{
\method{summary}{freekt}(object, \dots)
}

\arguments{
  \item{object}{
   An object of class "\code{freekt}" obtained by using one of the fitting algorithms.
}
  \item{\dots}{
  Additional arguments to be passed to summary.freekt.  Currently ignored.} 
}

\value{
A table containing the values of the optimal amount of penalty (when applicable), the optimal knots, sum of squared errors (SSE), and generalized cross-validation (GCV).
}

\author{
Steven Spiriti
}

\seealso{
\code{\link{freeknotfit}} for the fitting algorithms.
}

\examples{
x <- 0:30/30
truey <- x*sin(10*x)
set.seed(10556)
y <- truey + rnorm(31, 0, 0.2)
xy.freekt <- freelsgen(x, y, degree = 2, numknot = 2, 555)
summary(xy.freekt)
}

\keyword{
nonparametric
regression
smooth
}


