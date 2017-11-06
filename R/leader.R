leader <-
function(data, radius = NULL, scale = T) 
{
#  require(FNN)

unitize <-
function(z) {
        zrange <- range(z)
        if (!(dif <- diff(zrange))) 
            return(rep(0, length(z)))
        (z - zrange[1])/dif
    }

  data <- if (scale) apply(as.matrix(data),2,unitize) else as.matrix(data)

  n <- nrow(data)
  p <- ncol(data)

  LWradius <- function(n,p) 0.1/(log(n)^(1/p))
  if (is.null(radius)) radius <- LWradius(n,p)

# radius <- sort(unique(pmax(radius,0)))
  radius <- sort(pmax(radius,0))

  if (length(radius) == 1 && radius == 0) {
    return(list(list(radius = 0, leaders = data)))
  }

  rad0 <- radius[1] == 0

  nrad <- length(radius)
  LIST <- rep( list(list( radius = NA, 
                     count = 1,
                     partitions = list(1),
                     leaders = 1, 
                     centroids = list(data[1,]),
                     variances = list(rep(0,p)),
                     ranges = list(min = list(data[1,]), max = list(data[1,])),
                     maxdist = 0)), length(radius))



  if (radius[1] == 0) {

   LIST[[1]]$radius <- 0
   LIST[[1]]$partitions <- 1:n
   LIST[[1]]$leaders <- data
   LIST[[1]]$count <- LIST[[1]]$maxdist <- LIST[[1]]$centroids <- NULL
   LIST[[1]]$variances <- LIST[[1]]$ranges <- NULL
   irad <- 2

  }
  else irad <- 1

  for (i in 2:n) {
# one pass through data
     for (r in irad:nrad) {
        KNN <- get.knnx(data = data[LIST[[r]]$leaders, , drop = F], 
                        query = data[i, , drop = F], k = 1)
        m <- KNN$nn.index[1, 1]
        d <- KNN$nn.dist[1, 1]
        if (d < radius[r]) {
# new member
          l <- LIST[[r]]$leaders[m]
          j <- LIST[[r]]$count[m]
          LIST[[r]]$ranges$min[[m]] <- pmin(data[i,],LIST[[r]]$ranges$min[[m]])
          LIST[[r]]$ranges$max[[m]] <- pmax(data[i,],LIST[[r]]$ranges$max[[m]])
          LIST[[r]]$count[m] <- LIST[[r]]$count[m] + 1
          e <- (data[i,] - LIST[[r]]$centroids[[m]])
          h <- e/(j+1)
          LIST[[r]]$centroids[[m]] <- LIST[[r]]$centroids[[m]] + h
          LIST[[r]]$variances[[m]] <- LIST[[r]]$variances[[m]] + (j/(j+1))*(e*e)
          LIST[[r]]$partitions[[m]] <- c(LIST[[r]]$partitions[[m]], i)
          LIST[[r]]$maxdist[m] <- max(LIST[[r]]$maxdist[m],d)
          next
        }
# new exemplar
      LIST[[r]]$radius <- radius[r]
      LIST[[r]]$count <- c(LIST[[r]]$count,1)
      LIST[[r]]$maxdist <- c(LIST[[r]]$maxdist,0)
      LIST[[r]]$ranges$min <- c(LIST[[r]]$ranges$min, list(data[i,]))
      LIST[[r]]$ranges$max <- c(LIST[[r]]$ranges$max, list(data[i,]))
      LIST[[r]]$centroids <- c(LIST[[r]]$centroids, list(data[i,]))
      LIST[[r]]$variances <- c(LIST[[r]]$variances, list(rep(0,p)))
      LIST[[r]]$leaders <- c(LIST[[r]]$leaders, i)
      LIST[[r]]$partitions <- c(LIST[[r]]$partitions, list(i))
   }}


   LIST[irad:nrad] <- lapply( LIST[irad:nrad], function(x) {
       x$ranges$min <- drop(do.call( "rbind", x$ranges$min))
       x$ranges$max <- drop(do.call( "rbind", x$ranges$max))
       x$centroids <- drop(do.call( "rbind", x$centroids))
       for (k in 1:length(x$count)) if (x$count[k] > 1) x$variances[[k]] <- x$variances[[k]]/(x$count[k] - 1)
       x$variances <- drop(do.call("rbind", x$variances))
       x <- x[names(x) != "count"] 
       x
      })
  
   LIST  
}
