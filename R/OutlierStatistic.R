OutlierStatistic <-
function( x, nproj=1000, prior=NULL, seed=NULL) {

   x <- as.matrix(x)
   n <- nrow(x)
   p <- ncol(x)

# don't reset the seed unless multivariate ...

   if (p == 1) nproj <- 1 else if (!is.null(seed)) set.seed(seed)

   if (is.null(prior)) prior <- rep(-Inf,n)

   stat <- apply( apply( x %*% matrix( rnorm(nproj*p), p, nproj), 2, 
       function(vTx) {
                      MED <- median(vTx)
                      MAD <- median(abs(vTx - MED))
                      if (!MAD) return(rep(Inf,n)) else abs(vTx - MED)/MAD
               }), 1, max)

  pmax(stat,prior)
}
