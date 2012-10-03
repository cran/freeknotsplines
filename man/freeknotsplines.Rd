\name{freeknotspline-package}
\Rdversion{1.1}
\alias{freeknotspline-package}
\docType{package}

\title{
Free-Knot Splines
}

\description{
This package is for fitting free-knot splines for data with one independent variable and one dependent variable.
}

\details{
Four fitting methods are included for the case where the number of knots is known in advance; for details, see \code{\link{freeknotfit}}.  In addition, methods are available to compute the fitted values, the residuals, and the coefficients of the splines, and to plot the results, along with a method to summarize the results.  Finally, a function (see \code{\link{fit.search.numknots}}) is provided to optimize the number of knots using a specified fit criterion.  Several fit criteria are provided (see \code{\link{fitcriteria}}).
}

\author{
Steven Spiriti, Philip Smith, and Pierre Lecuyer
}

\references{
Spiriti, S., Eubank, R., Smith, P., Young, D., "Knot Selection for Least-Squares and Penalized Splines," \emph{Journal of Statistical Computation and Simulation}, in press.
}

\keyword{
package
nonparametric
regression
smooth
}
