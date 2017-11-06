\name{OutlierStatistic}
\alias{OutlierStatistic}

\title{Nonparametric Outlier Statistic}

\description{
   Robust nonparametric outlier statistic for univariate or multivariate data.
}
\usage{
OutlierStatistic( x, nproj=1000, prior=NULL, seed=NULL)
}
\arguments{
  \item{x}{
     A numeric vector or matrix for which the outlier statistic is to
     be determined.
  }
  \item{nproj}{
     If \code{x} is multivariate, the number of random projections to be used
     in computing the statistic.
  }
  \item{prior}{
     If \code{x} is multivariate, a prior estimate of the statistic for each 
     observation in \code{x}, to be used as a base line for maximization 
     relative to new random projections.
  }
 \item{seed}{
    An optional integer argument to \code{set.seed} for reproducible
    simulations. By default the current seed will be used.
    Reproducibility can also be achieved by calling \code{set.seed}
    before calling \code{OutlierStatistic}.
  }
}
\value{
   A vector giving the maximum value of the outlier statistic for each 
   observation over all projections.
}
\note{
Note that partition probabilities are computed from an exponential distribution
fit to the outlier statistic, rather than from the empirical distribution
of the outlier statistic.
}
\section{References}{
 W. A. Stahel, \emph{Breakdown of Covariance Estimators}, PhD thesis,
 Fachgruppe Fur Statistik, Eidgenossische Technische Hochshule (ETH), 1981.

 D. L. Donoho, \emph{Breakdown Properties of Multivariate Location Estimators},
 PhD thesis, Department of Statistics, Harvard University, 1982.
}
\seealso{
  \code{\link{partProb}}
}
\examples{

 stat <- OutlierStatistic(faithful)
 q.99 <- quantile(stat,.99)
 out <- stat > q.99

 plot( faithful[,1], faithful[,2], 
       main="red : .99 quantile for outlier statistic", cex=.5)
 points( faithful[out,1], faithful[out,2], 
         pch = 4, col = "red", lwd = 1, cex = .5)

 require(mvtnorm)

 set.seed(0)
 Sigma <- crossprod(matrix(rnorm(2*2),2,2))
 x <- rmvt( 10000, sigma = Sigma, df = 2) 

 stat <- OutlierStatistic(x)
 q.95 <- quantile(stat,.95)

 hist(x, main = "gray : .95 quantile for outlier statistic", col = "black")
 abline( v = x[stat > q.95], col = "gray")
 hist(x, col = "black", add = TRUE)
}

\keyword{nonparametric}
