allProb <-
function( leaderInstance, partprob) {

 if (is.null(leaderInstance$radius)) {
      if (length(leaderInstance) == 1) {
        leaderInstance <- leaderInstance[[1]]
        if (is.null(leaderInstance$radius)) stop("leaderInstance improperly specified")
      }
      else {
        stop("leaderInstance improperly specified")
      }
    }

mapx <- function(leaderInstance)
{
# maps data to partitions
  l <- length(unlist(leaderInstance$partitions))
  z <- numeric(length(unlist(leaderInstance$partitions)))
  for (i in 1:length(leaderInstance$partitions)) {
     z[leaderInstance$partitions[[i]]] <- i
  }
  z
}

  if (!is.null(dim(partprob))) partprob <- apply(partprob,1,median)

  partprob[mapx(leaderInstance)]
}
