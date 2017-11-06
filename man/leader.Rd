\name{leader}
\alias{leader}

\title{Leader Algorithm for Data Partitioning}

\description{
  Partitions the data according to Hartigan's leader algorithm, and 
  provides ranges, centroids, and variances for the partitions.
}
\usage{
leader(data, radius = NULL, scale = T)
}
\arguments{
  \item{data}{
    A numeric vector or matrix of observations. If a matrix, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{radius}{
     A vector of values for the partitioning radius. Wilkinson's default
     radius is used if \code{radius} is left unspecified (see function 
\code{LWradius}).
  }
  \item{scale}{
     A logical variable indicating whether or not the data should be mapped
     to the unit hypercube. The default is to scale the data. Values of the
     radius will not be scaled; they should be specifed relative to the unit 
     hypercube unless \code{scale = F}.
   }
}
\value{
 A list with one component for each value of \code{radius}, each having the
 following sub-components:

  \item{radius}{The value of the radius associated with the partitioning.}
  \item{partitions}{
     A list with one component for each partition, giving the indexes (as
     observations in the data) of the members of the partition. The first 
     index is that of the associated \emph{leader} (sometimes called 
     \emph{exemplar}).
  }
  \item{leaders}{
    The indexes of the leaders for each partition.
  }
  \item{centroids}{
    The centroids for each partition, as a matrix with rows corresponding to
    the partitions and columns corresponding to variables if multidimensional.
    These will be the data if \code{radius == 0}.
  }
  \item{variances}{
    The variances for each partition, as a matrix with rows corresponding to
    the partitions and columns corresponding to variables if multidimensional.
  }
  \item{ranges}{
    A list with two components: \code{min} and \code{max} giving the
    minimum and maximum values for each variable for each partition.
    These \code{range} components are given as a matrix with rows 
    corresponding to the partitions and columns corresponding to variables 
    if multidimensional.
  }
  \item{maxdist}{
    A vector with one value for each partition, giving the largest distance
    from each leader to any member of its partition.
  }
}
\details{
   Given a partitioning radius \code{r}, the leader algorithm makes one pass 
   through the data, designating an observation as a new leader if it is not 
   within \code{r} of an existing leader, and otherwise assigning it to the
   partition associated with the nearest existing leader. The set of leaders 
   typically depends on the order of the data observations. \cr
   If \code{radius = 0}, then all of the data observations are leaders, and
   only \code{radius} and \code{leaders} are returned as output components. \cr
   This implementation does a completely new nearest-neighbor search for
   each observation and for each radius. A more efficient approach would be to 
   maintain, for each radius, a data structure (such as a kd-tree) allowing 
   fast nearest-neighbor search. These data structures could then be updated 
   to account for new observations. Currently, there doesn't seem to be a way 
   to do this in R.
}
\section{References}{
  J. A. Hartigan, \emph{Clustering Algorithms}, Wiley, 1975.
  
  L. Wilkinson, Visualizing Outliers, Technical Report, University of
  Illinois at Chicago, 2016.
  \code{https://www.cs.uic.edu/~wilkinson/Publications/outliers.pdf}
}
\seealso{
  \code{\link{LWradius}}
}
\examples{

 radius.default <- LWradius(nrow(faithful),ncol(faithful))
 lead <- leader(faithful, radius = c(0,radius.default))

# number of partitions for each radius
 sapply(lead, function(x) length(x$partitions))

# plot the leaders for the non-zero radius
 plot( faithful[,1], faithful[,2], 
       main = "blue indicates leaders (default radius)", 
       pch = 16, cex = .5)
 ldrs <- lead[[2]]$leaders
 points( faithful[ldrs,1], faithful[ldrs,2], 
         pch = 8, col = "dodgerblue", cex = .5)

}

\keyword{cluster}
