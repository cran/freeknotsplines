\name{fitcriteria}

\Rdversion{1.1}

\alias{fitcriteria}
\alias{AIC.freekt}
\alias{AICc.freekt}
\alias{BIC.freekt}
\alias{adjAIC.freekt}
\alias{adjGCV.freekt}

\title{
Fit Criteria for Free-Knot Splines
}

\description{
These functions compute various criteria for determining the fit of a free-knot spline.  \code{AIC.freekt} computes the Akaike Information Criterion, with \code{k} determining the amount of the penalty.  \code{AICc.freekt} computes the corrected Akaike Information Criterion.  \code{BIC.freekt} computes the Bayesian Information Criterion, also known as Schwarz Information Criterion.  \code{adjAIC.freekt} computes an adjusted Akaike Information Criterion with the penalty increased to account for the greater flexibility of free knots.  \code{adjGCV.freekt} computes an adjusted GCV with the degrees of freedom increased to account for the greater flexibility of free knots.}

\usage{
\method{AIC}{freekt}(object, \dots, k = 2)
AICc.freekt(object)
\method{BIC}{freekt}(object, \dots)
adjAIC.freekt(object)
adjGCV.freekt(object, d = 3)
}

\arguments{
  \item{object}{
  An object of class "\code{freekt}" obtained by using one of the fitting algorithms.
}
  \item{k}{
  The amount of the penalty.  Used only for \code{AIC.freekt}.
}
  \item{d}{
  The amount of the penalty.  Used only for \code{adjGCV.freekt}.
} 
  \item{\dots}{
  Additional arguments to be passed to the \code{AIC.freekt} and \code{BIC.freekt} functions.
}
}

\value{
Returns the value of the specified fit criterion.
}

\references{
Spiriti, S., Eubank, R., Smith, P., Young, D., "Knot Selection for Least-Squares and Penalized Splines," \emph{Journal of Statistical Computation and Simulation}, in press.
}

\author{
Steven Spiriti
}

\seealso{
\code{\link{fit.search.numknots}}, which uses these fit criteria to determine the number of knots.
}

\keyword{
nonparametric
regression
smooth
}
