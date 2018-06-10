\name{freeknotfit}

\Rdversion{1.1}

\alias{freeknotfit}
\alias{freelsgen}
\alias{freelsgold}
\alias{freepsgen}
\alias{freepsgold}

\title{
Fit Free-Knot Splines To Data
}

\description{
These functions fit free-knot splines to data with one independent variable and one dependent variable.  It is assumed that the number of knots is known in advance.  \code{freelsgen} and \code{freelsgold} fit least-squares splines with no penalty, while \code{freepsgen} and \code{freepsgold} fit penalized splines.  \code{freelsgen} and \code{freepsgen} use a genetic algorithm, while \code{freelsgold} and \code{freepsgold} use a blind search augmented with a golden section algorithm.
}

\usage{
freelsgen(x, y, degree, numknot, seed = 5, stream = 0)
freelsgold(x, y, degree, numknot, seed = 5, stream = 0)
freepsgen(x, y, degree, numknot, seed = 5, stream = 0)
freepsgold(x, y, degree, numknot, seed = 5, stream = 0)
}

\arguments{
  \item{x}{
  A vector containing the values of the independent variable.
}
  \item{y}{
  A vector containing the values of the dependent variable.
}
  \item{degree}{
  The degree of the spline fit.
}
  \item{numknot}{
  The number of knots.
}
  \item{seed}{
  The value of the initial seed.  Defaults to 5.
}
  \item{stream}{
  The value of the initial stream to be used for parallel programming.
  Defaults to 0.
}
}

\value{
An object of class "\code{freekt}" containing the following components: 

\item{x}{A vector containing the x values.}
\item{y}{A vector containing the y values.}
\item{degree}{The degree of the spline fit.}
\item{seed}{The value of the initial seed.}
\item{stream}{The value of the stream.}
\item{lambda}{The optimum amount of penalty.  This is automatically equal to 0 for \code{freelsgen} and \code{freelsgold}.}
\item{optknot}{A vector containing the optimal knots.}
\item{tracehat}{The trace of the hat matrix for the optimal fit.}
\item{GCV}{The value of generalized cross validation (GCV) for the optimal fit.}
\item{GSJS}{The GSJS estimator, an estimator of the variance of the data.}
\item{call}{The function call.}
}

\references{
Eubank, R. (1999), \emph{Nonparametric Regression and Spline Smoothing}, New York: Marcel Dekker, Inc., Second ed.\cr
\cr
Spiriti, S., Eubank, R., Smith, P., Young, D., "Knot Selection for Least-Squares and Penalized Splines," \emph{Journal of Statistical Computation and Simulation}, in press.
}

\author{
Steven Spiriti, Philip Smith, and Pierre Lecuyer
}

\seealso{
\code{\link{fit.search.numknots}} for the case where the number of knots is not specified in advance.
}

\examples{
x <- 0:30/30
truey <- x*sin(10*x)
set.seed(10556)
y <- truey + rnorm(31, 0, 0.2)
xy.freekt <- freelsgen(x, y, degree = 2, numknot = 2, 555)
plot.freekt(xy.freekt, xfit = 0:1000/1000)
}

\keyword{
nonparametric
regression
smooth
}
