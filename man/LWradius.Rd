\name{LWradius}
\alias{LWradius}

\title{Wilkinson's default leader-partitioning radius}

\description{
   Wilkinson's default leader-partitioning radius.
}
\usage{
LWradius( n, p)
}
\arguments{
  \item{n}{
     The number of observations (rows) in the data.
  }
  \item{p}{
     The number of variables (columns) in the data; \code{p = 1} if
     univariate.
  }
}
\value{
   Wilkinson's default leader partitioning radius \code{0.1/(log(n)^(1/p))}.
}
\section{References}{
  L. Wilkinson (2016), Visualizing Outliers, Technical Report, University of
  Illinois at Chicago, 
  \code{https://www.cs.uic.edu/~wilkinson/Publications/outliers.pdf}.
}
\seealso{
  \code{\link{leader}}
}
\examples{

 x1 <- rnorm(10000)
 LWradius(length(x1),1)

 LWradius(nrow(faithful),ncol(faithful))
}

\keyword{misc}
