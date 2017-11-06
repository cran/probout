logdens <-
function( x, simData, shrink = 1) 
{

logdensV <- 
function( data, prob, centroids, variances) 
{
    n <- length(data)
    k <- length(prob)
    out <- .Fortran("lgd1v", as.double(data), as.double(prob), 
        as.double(centroids), as.double(variances), 
        as.integer(n), as.integer(k), 
        z = double(n*k), hood = double(max(n,k)), PACKAGE = "probout")

#   matrix(out$z,n,k)
    out$hood[1:n]
}

logdensVII <-
function (data, prob, centroids, variances) 
{
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    k <- length(prob)
    out <- .Fortran("lgdvii", as.double(data), as.double(prob), 
        as.double(centroids), as.double(variances), 
        as.integer(n), as.integer(p), as.integer(k), 
        z = double(n*k), hood = double(max(n,k)), PACKAGE = "probout")

#    matrix(out$z,n,k)
     out$hood[1:n]
}

logdensVVI <- 
function (data, prob, centroids, scale, shape) 
{
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    k <- length(prob)
    out <- .Fortran("lgdvvi", as.double(data), as.double(prob), 
        as.double(centroids), as.double(scale), as.double(shape),
        as.integer(n), as.integer(p), as.integer(k), 
        z = double(n*k), hood = double(max(n,k)), PACKAGE = "probout")

#    matrix(out$z,n,k)
     out$hood[1:n]
}

   if (simData$radius == 0 || is.null(simData$index)) {
     stop("radius = 0 or nsim = 0: no log density for simData")
   }

   p <- if (is.null(dim(x))) 1 else ncol(x)

   shrink <- shrink*shrink

   if (p == 1) {
     return(logdensV( x, prob = simData$weight,
                      centroids = drop(simData$location),
                      variances = shrink*simData$scale))
   }

   if (is.null(simData$shape)) {
     logdens <- logdensVII( x, prob = simData$weight,
                            centroids = t(simData$location),
                            variances = shrink*simData$scale)
   }
   else {
     logdens <- logdensVVI( x, prob = simData$weight,
                            centroids = t(simData$location),
                            scale = rep(1,length(simData$weight)), 
                            shape = shrink*simData$shape)
   }

 logdens
}
