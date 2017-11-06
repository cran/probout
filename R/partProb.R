partProb <-
function( simData, method = c("intrinsic","distance","logdensity","distdens",
          "density"), shrink = 1, nproj = 1000, seed = NULL)
{

pseudoData <- function( simData, shrink = 1)
{
 if (is.null(dim(simData$location)) || ncol(simData$location) == 1) {
   c(drop(simData$location), 
     drop(simData$location)[simData$index]+shrink*drop(simData$offset))
 }
 else {
   rbind(simData$location, 
   simData$location[simData$index,,drop=F]+shrink*as.matrix(simData$offset))
 }
}

 nex <- if (is.null(dim(simData$location))) {
          length(simData$location) 
        }
        else {
          nrow(simData$location)
        }
 

stat <- 
switch(method[1],
       "intrinsic" = {
        if (simData$radius != 0 && !is.null(simData$index)) {
           x <- pseudoData(simData,shrink=shrink)
           OutlierStatistic(x,nproj=nproj,seed=seed)
         }
         else {
           OutlierStatistic(simData$location,nproj=nproj,seed=seed)
         }
        },
       "distance" = {
        if (simData$radius != 0 && !is.null(simData$index)) {
           x <- pseudoData(simData,shrink=shrink)
           d <- get.knnx(simData$location,query=x,k=2)$nn.dist[,2]
           OutlierStatistic(d,seed=seed)
         }
         else {
           dex <- knn.dist(simData$location,1)
           OutlierStatistic(dex,nproj=nproj,seed=seed)
         }
        },
       "density" = {
        if (simData$radius != 0 && !is.null(simData$index)) {
           x <- pseudoData(simData,shrink=shrink)
           den <- exp(logdens(x, simData))
           KNN2 <- get.knnx(data = simData$location, query = x, k = 2)
           z <- numeric(length(den))
           for (j in seq(along = den)) {
              r <- range(c(den[j],den[KNN2$nn.index[j,2]]))
              z[j] <- min(r)/max(r)
           }
           OutlierStatistic(z,seed=seed)
         }
         else {
           stop("don't have a density estimate for this case")
           OutlierStatistic( simData$location, seed=seed)
         }
        },
       "logdensity" = {
        if (simData$radius != 0 && !is.null(simData$index)) {
           x <- pseudoData(simData,shrink=shrink)
           logden <- logdens(x, simData)
           KNN2 <- get.knnx(data = simData$location, query = x, k = 2)
           z <- numeric(length(logden))
           for (j in seq(along = logden)) {
              z[j] <- abs(logden[j] - logden[KNN2$nn.index[j,2]])
           }
           OutlierStatistic( z, seed=seed)
         }
         else {
           stop("don't have a density estimate for this case")
           OutlierStatistic( simData$location, seed=seed)
         }
        },
       "distdens" = {
        if (simData$radius != 0 && !is.null(simData$index)) {
           x <- pseudoData(simData,shrink=shrink)
           logden <- logdens(x, simData)
           KNN2 <- get.knnx(data = simData$location, query = x, k = 2)
           z <- numeric(length(logden))
           for (j in seq(along = logden)) {
              z[j] <- abs(logden[j] - logden[KNN2$nn.index[j,2]])
           }
           d <- get.knnx(simData$location,query=x,k=2)$nn.dist[,2]
           OutlierStatistic( cbind(d,z), nproj, seed=seed)
         }
         else {
           stop("don't have a density estimate for this case")
           OutlierStatistic( simData$location, seed=seed)
         }
        },
         stop("method not recognized")
       )


#emcdf <- function (x) 
#{
# x <- sort(x)
# sapply( x, function(z) sum(x <= z))/length(x)
#}
#Pstat <- empcdf(stat)
#ord <- order(stat)

# exponential distribution of outlier statistic
# require(MASS)
 fit <- fitdistr( stat, "exponential")
 pex <- pexp(stat[1:nex],rate=fit[[1]])

 pex
}
