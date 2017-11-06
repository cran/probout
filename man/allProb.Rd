\name{allProb}
\alias{allProb}

\title{Outlier probabilities for all observations}

\description{
   Outlier probabilities for all of  the data, obtained by assigning
   to each observation the probabilty of the its associated leader partition.
}
\usage{
allProb( leaderInstance, partprob)
}
\arguments{
  \item{leaderInstance}{
    A single component from a call to \code{leader}, giving leader algorithm
    results for one value of the partitioning radius.
  }
  \item{partprob}{
     A vector of probabilities for each partition in \code{leaderInstance}.
  }
}
\value{
  A vector of probabilities for each observation in the data underlying
  \code{leaderInstance}. Each observation inherits the probability of its
  associated partition.
}
\seealso{
  \code{\link{leader}},
  \code{\link{partProb}}
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

 quantile(partprobs)

# plot leaders with outlier probability > .95
 plot( faithful[,1], faithful[,2], pch = 16, cex = .5,
       main = "red : instances with outlier probability > .95")
 allprobs <- allProb( lead[[1]], partprobs)
 out <- allprobs > .95
 points( faithful[out,1], faithful[out,2], pch = 8, cex = 1, col = "red")

}

\keyword{misc}
