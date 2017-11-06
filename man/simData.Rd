\name{simData}
\alias{simData}

\title{Simulates observations for outlier determination.}

\description{
  Simulates observations from a mixture model based on information on
  partitions from the \code{leader} function.
}
\usage{
simData( leaderInstance, nsim=NULL, model=c("diagonal","spherical"), seed=NULL)
}
\arguments{
 \item{leaderInstance}{
    A	single component from a call to \code{leader}, giving Leader Algorithm
    results for one value of the partitioning radius.
  }
  \item{nsim}{
    The number of observations to be simulated. Only the radius and centroids
    are returned of \code{nsim = 0} or \code{leaderInstance$radius == 0})
    --- no observations are simulated. \cr
    Default: \code{min(# observations,max(# partitions, 1000))}.
  }
  \item{model}{
    For multivariate data, a vector of character strings indicating the type of
    Gaussian mixture model covariance to be used in generating the simulated 
    observations (see \code{details}). \cr
    For univariate data, the observations are generated from a model
    in which the variances may vary across components.
  }
 \item{seed}{
    An optional integer argument to \code{set.seed} for reproducible
    simulations. By default the current seed will be used.
    Reproducibility can also be achieved by calling \code{set.seed}
    before calling \code{simData}.
  }
}
\value{
 A list with the following components:
  \item{radius}{The value of the radius associated with \code{leaderInstance}.}
  \item{location}{
     The vector or matrix of centroids of the partitions. If a matrix,
     rows correspond to the partitions and columns to the variables.
  }
  \item{index}{
    A vector of integer values giving the index of the partition associated 
    with each simulated observation.
  }
  \item{offset}{
    A vector of numeric values giving offset for the simulated observations 
    from their associated centroids.
  }
  \item{weight}{
    A vector of numeric values between 0 and 1 giving the proportion of data 
    observations in each partition.
  }
  \item{scale}{
     The scale (variance) of the mixture components in a univariate or
     spherical model. Set to 1 for each component in the diagonal model.
  }
  \item{shape}{
     A matrix giving the variances of the mixture component in a diagonal
     model. The rows correspond to the dimensions of the data, while the
     columns correspond to the mixture components (partitions).
  }
}
\details{
The following models are available for multivariate data:
 \tabular{lcl}{
    \code{"spherical"} \tab : \tab spherical, varying volume \cr
    \code{"diagonal"} \tab : \tab diagonal, varying volume and shape \cr
   }
An ellipsoidal model is also possible, but has not yet been implemented. \cr
If \code{nsim = 0} or \code{leaderInstance$radius == 0}, no observations are 
simulated, and only the radius and partition centroids are returned.
}
\section{References}{
  C. Fraley, Estimating Outlier Probabilities for Large Datasets, 2017.
  }
\seealso{
  \code{\link{leader}},
  \code{\link{partProb}}
}
\examples{

 radius.default <- LWradius(nrow(faithful),ncol(faithful))
 lead <- leader(faithful, radius = c(0,radius.default))

# (simulated) data for outlier statistic (no simulation for radius = 0)
 sim <- lapply( lead, simData)

# components of simData output
 lapply( sim, names)

}

\keyword{datagen}
