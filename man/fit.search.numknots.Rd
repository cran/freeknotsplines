\name{fit.search.numknots}

\Rdversion{1.1}

\alias{fit.search.numknots}

\title{
Perform a Search on the Number of Knots and Fit Free-Knot Splines To Data Using the Optimal Number of Knots
}

\description{
This function fits free-knot splines to data using every value for the number of knots between \code{minknot} and \code{maxknot}.  The number of knots is then chosen to optimize a fit criterion.  The free-knot spline with the optimum number of knots is returned.
}

\usage{
fit.search.numknots(x, y, degree, minknot = 1, maxknot = 5, 
                   alg = "LS", search = "genetic",
                   knotnumcrit = "adjGCV", k = 2, d = 3, seed = 5,
                   stream = 0)
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
  \item{minknot}{
  The minimum number of knots to search.  Defaults to 1.
}
  \item{maxknot}{
  The maximum number of knots to search.  Defaults to 5.
}
  \item{alg}{
  The spline-fitting algorithm.  Choices are "LS" for least-squares and "PS" for P-splines.  Defaults to "LS."
}
  \item{search}{
  The random search algorithm.  Choices are "genetic" for a genetic algorithm and "golden" for a blind random search with golden section adjustment.  Defaults to "genetic."
}
  \item{knotnumcrit}{
  The criterion to be used for determining the number of knots.  Choices are "GCV" for generalized cross-validation, "AIC" for the Akaike information criterion, "AICc" for corrected Akaike information criterion, "BIC" for Bayesian information criterion, "adjAIC" for an adjusted version of the Akaiki information criterion, and "adjGCV" for an adjusted version of generalized cross-validation.  Defaults to "adjGCV."  
}
  \item{k}{The amount of penalty when AIC is used.  Has no effect with criteria other than AIC.  Defaults to 2.
}
  \item{d}{The amount of penalty when adjGCV is used.  Has no effect with criteria other than adjGCV.  Defaults to 3.
}
  \item{seed}{
  The value of the initial seed.  Defaults to 5.
}
  \item{stream}{
  The value of the initial stream to be used for parallel programming.  Defaults to 0.
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
\item{GCV}{The value of generalized cross validation (GCV) for the optimal fit.)}
\item{GSJS}{The GSJS estimator, an estimator of the variance of the data.}
\item{call}{The function call.}
}

\references{
Eubank, R. (1999), \emph{Nonparametric Regression and Spline Smoothing}, New York: Marcel Dekker, Inc., Second ed.\cr
\cr
Spiriti, S., Eubank, R., Smith, P., Young, D., "Knot Selection for Least-Squares and Penalized Splines," \emph{Journal of Statistical Computation and Simulation}, in press.
}

\author{
Steven Spiriti
}

\seealso{
\code{\link{fitcriteria}} for the fit criteria, \code{\link{freeknotfit}} for the free-knot spline algorithms.
}

\examples{
x <- 0:30/30
truey <- x*sin(10*x)
set.seed(10556)
y <- truey + rnorm(31, 0, 0.2)
xy.freekt <- fit.search.numknots(x, y, degree = 2, minknot = 1, maxknot = 3, seed = 555)
plot.freekt(xy.freekt, xfit = 0:1000/1000)
}

\keyword{
nonparametric
regression
smooth
}
