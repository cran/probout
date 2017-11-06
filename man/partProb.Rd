\name{partProb}
\alias{partProb}

\title{Partition outlier probabilities}

\description{
   Assigns outlier probabilities to the partitions by fitting an exponential 
   distribution to a nonparametric outlier statistic for simulated
   data or partition centroids.
}
\usage{
partProb( simData, method = c("intrinsic","distance","logdensity","distdens",
          "density"), shrink = 1, nproj = 1000, seed = NULL)
}
\arguments{
  \item{simData}{
     Observations from a call to \code{simData}, which includes
     the partition centroids and (optionally) simulated data as well.
  }
 \item{method}{One of the following options:
\tabular{lcl}{
  \code{"intrinsic"} \tab : \tab outlier statistic applied to simulation data 
                     (centroids if no simulation) \cr
  \code{"distance"} \tab : \tab outlier statistic applied to distances between
                      NN partitions \cr
  \code{"logdensity"} \tab : \tab outlier statistic applied to
             differences in log density between NN partitions\cr
  \code{"distdens"} \tab : \tab outlier statistic applied to a matrix
        consisting of the \code{"distance"} and
             \code{"logdensity"} values \cr
  \code{"density"} \tab : \tab outlier statistic applied to
             smallest/largest ratios of density between NN partitions
   }
The default is to use the \code{"intrinsic"} method.
 }
  \item{shrink}{
     Shrinkage parameter for outlier detection data. The offsets from
     \code{simData} are scaled by this factor before adding them to the
     partition centroids as data for outlier detection. The default
     value is \code{shrink = 1}, so that no shrinkage is applied to
     simulation offsets.
 }
 \item{nproj}{
    If the data is multivariate or \code{method = "distdens"}, 
    the number of random projections to be used
    to obtain the outlier statistic.
 }
 \item{seed}{
    An optional integer argument to \code{set.seed} for reproducible
    outlier statistics. By default the current seed will be used.
    Reproducibility can also be achieved by calling \code{set.seed}
    before calling \code{partProb}.
 }
}
\value{
 A vector of probabilities for each partition, obtained by fitting an
 exponential distribution to the outlier statistic.
}
\details{
\code{"logdensity"} is generally prefered over \code{"density"}, because
negative values that are large in magniude
of the logarithm of the density will not be 
numerically distinguishable as density values.
}
\section{References}{
  C. Fraley, Estimating Outlier Probabilities for Large Datasets, 2017.
  }
\seealso{
  \code{\link{simData}},
  \code{\link{OutlierStatistic}},
  \code{\link{allProb}}
}
\examples{

 set.seed(0)

 lead <- leader(faithful)
 nlead <- length(lead[[1]]$partitions)

# repeat multiple times to account for randomness
 ntimes <- 100
 probs <- matrix( NA, nlead, ntimes)
 for (i in 1:ntimes) {
    probs[,i] <- partProb( simData(lead[[1]]), method = "distance")
 }

# median probability for each partition
 partprobs <- apply( probs, 1, median)

 quantile(probs)

# plot leaders with outlier probability > .95
 plot( faithful[,1], faithful[,2], pch = 16, cex = .5,
       main = "red : leaders with outlier probability > .95")
 out <- partprobs > .95
 l <- lead[[1]]$leaders
 points( faithful[l[out],1], faithful[l[out],2], pch = 8, cex = 1, col = "red")

}

\keyword{nonparametric}
