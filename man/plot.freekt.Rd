\name{plot.freekt}

\Rdversion{1.1}

\alias{plot.freekt}

\title{
Plot Fitted Values for Free-Knot Spline
}

\description{
This function plots the fit obtained using a free-knot spline.
}

\usage{
\method{plot}{freekt}(x, xfit = x@x, linecolor = "blue", lwd = 1, lty = 1, \dots)
}

\arguments{
  \item{x}{
  An object of class "\code{freekt}" obtained by using one of the fitting algorithms.
}
  
  \item{xfit}{
  A vector of x values at which to plot the fitted values.  Defaults to the x values of the data.
}

  \item{linecolor}{
  The color of the line.  Defaults to blue.
}

  \item{lwd}{
  The line width.  It is passed to the \code{lines} function.  Defaults to 1.
}

  \item{lty}{
  The line type.  It is passed to the \code{lines} function.  Defaults to 1.
}

  \item{\dots}{
  Additional arguments to be passed to the \code{plot} function.
}
}

\value{
A plot of the data, together with the spline estimator.
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
plot(xy.freekt, xfit = 0:1000/1000)
}

\keyword{
nonparametric
regression
smooth
}
