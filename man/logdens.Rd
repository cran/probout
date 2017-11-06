\name{logdens}
\alias{logdens}

\title{Log Density for Gaussian Mixture Model}

\description{
   Computes the log density for observations in a univariate or multivariate
   Gaussian mixture model with spherical or diagonal (co)variance that 
   varies across components.
}
\usage{
logdens( x, simData, shrink = 1)
}
\arguments{
  \item{x}{
     A numeric vector or matrix for which the log density is to be computed.
  }
 \item{simData}{
     Observations from a call to \code{simData}, which includes
     the partition centroids and variance information for the underlying
     simulation model.
  }
 \item{shrink}{
     Shrinkage parameter for the mixture model variance. To be
     consistent with the shrinkage as described in \code{partProb},
     the variance is scaled by the \emph{square} of \code{shrink}.
     The default value is \code{shrink = 1}, so that no shrinkage is 
     applied to the variance.
 }
}
\value{
   A vector giving the log density of \code{x} in the model as
   specified by \code{simData}, with optional shrinkage applied to the
   variance.
}
\details{
If either \code{radius = 0}, or \code{simData} returns only centroids 
(\code{nsim = 0}), then no density estimate is attempted.
}
\section{References}{
 G. Celeux and G. Govaert, Gaussian Parsimonious Mixture Models,
 \emph{Pattern Recognition}, 1995.

 G. J. McLachlan and D. Peel, \emph{Finite Mixture Models}, Wiley, 2000.
 
 C. Fraley and A. E. Raftery, Model-based clustering, discriminant analysis
 and density estimation, 
 \emph{Journal of the American Statistical Association}, 2002.
}
\seealso{
  \code{\link{partProb}}
}
\examples{

 lead <- leader(faithful)
 sim <- simData( lead)

 logdens( faithful, sim)

}

\keyword{models}
