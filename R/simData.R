simData <-
function( leaderInstance, nsim=NULL,  model=c("diagonal","spherical"), seed=NULL)
{
# assumes Euclidean distance

    if (is.null(leaderInstance$radius)) {
      if (length(leaderInstance) == 1) {
        leaderInstance <- leaderInstance[[1]]
        if (is.null(leaderInstance$radius)) stop("leaderInstance improperly specified")
      }
      else {
        stop("leaderInstance improperly specified")
      }
    }

    if (leaderInstance$radius == 0 || (!is.null(nsim) && nsim <= 0)) {
      return(list( radius = leaderInstance$radius, location = leaderInstance$centroids))
    }

    nex <- length(leaderInstance$partitions) # number of partitions


    n <- length(unlist(leaderInstance$partitions))

    if (is.null(nsim)) nsim <- 10*nex
    nsim <- min(n,max(nex,nsim,1000))

    p <- ncol(as.matrix(leaderInstance$centroids))

    weightFun <- function(x) {
           xmax <- max(abs(x))
           x <- x/xmax
           x/sum(x)
          }

    weight <- weightFun(sapply(leaderInstance$partitions, length))

    if (p > 1) {
      switch(model[1],
            "spherical" = {
                           scale <- apply(leaderInstance$variances, 1, mean)
                           scale <- pmax(scale,.Machine$double.eps)
                           shape <- NULL
                         },
            "diagonal" = {
                          scale <- NULL
                          shape <- apply( leaderInstance$variances, 1, 
                              function(x) pmax( x, .Machine$double.eps))
                         },
            stop("unrecognized model")
            )
      }
    else {
      scale <- pmax(leaderInstance$variances,.Machine$double.eps)
      shape <- NULL
    }

    if (!is.null(seed)) set.seed(seed)

#     require(mclust)

      if (p > 1) {

        X <- switch(model[1],
                   "spherical" = {
                                  simVII(list( pro = weight, 
                                         mean = t(leaderInstance$centroids),
                                   variance = list(sigmasq = scale)), nsim)
                                 },
                   "diagonal" = {
                                 simVVI(list( pro = weight, 
                                        mean = t(leaderInstance$centroids),
                      variance = list(scale = rep(1,length(weight)),
                                       shape = shape)), nsim)
                           },
               stop("unrecognized model")
          )

         E <- X[,-1,drop=F] - leaderInstance$centroids[X[,1],,drop = F]    

      }
      else  {

        X  <- simV(list( pro = weight, mean = drop(leaderInstance$centroids),
                      variance = list(sigmasq = scale)), nsim)

        E <- X[,-1] - drop(leaderInstance$centroids)[X[,1]]

      }

      list( radius = leaderInstance$radius,
            location = leaderInstance$centroids, index = X[,1], offset = E,
            weight = weight, scale = scale, shape = shape)

 }
