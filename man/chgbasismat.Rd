\name{chgbasismat}

\Rdversion{1.1}

\alias{chgbasismat}

\title{
Compute Change-of-Basis Matrix to Convert Truncated Power Basis to B-Spline Basis
}

\description{
This function computes a change-of-basis matrix that converts from the truncated power basis to the B-spline basis.  It is used by \code{coef.freekt} and should not be called directly by the user. 
}

\usage{
chgbasismat(knot, ord)
}

\arguments{
    \item{knot}{The vector of knots.  The first and last knots are repeated \code{ord} times to make the length equal to the dimension of the spline space.  
   }

    \item{ord}{The order of the spline, which is one more than the degree.
   }
}

\value{
The change-of-basis matrix that converts from the truncated power basis to the B-spline basis.
}

\references{
Smith, P. (1982), "Hypothesis Testing B-Spline Regression," \emph{Communications in Statistics, Series B}, 11, 143-157.\cr
\cr
Spiriti, S., Eubank, R., Smith, P., Young, D., "Knot Selection for Least-Squares and Penalized Splines," \emph{Journal of Statistical Computation and Simulation}, in press.
}

\author{
Steven Spiriti
}
