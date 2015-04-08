
getpeaks <- function(qtlprofiles, cutoff = 3.0){
  cat("Starting peak detection above",cutoff,"\n")
  peaks <- vector("list",nrow(qtlprofiles))
  for(x in 1:nrow(qtlprofiles)){
    peak <- FALSE
    curmax <- 0
    curmaxindex <- 1
    marker <- 1
    maximums <- NULL
    for(ab in (qtlprofiles[x,]>cutoff | qtlprofiles[x,]<(-cutoff))){
      if(ab){
        peak <- TRUE
        if(qtlprofiles[x,marker]/abs(qtlprofiles[x,marker]) > 0){
          if(qtlprofiles[x,marker] > curmax){
            curmax <- qtlprofiles[x,marker]
            curmaxindex <- marker
          }
        }else{
          if(qtlprofiles[x,marker] < (-curmax)){
            curmax <- qtlprofiles[x,marker]
            curmaxindex <- -marker
          }
        }
        if(ncol(qtlprofiles)==marker){
          if(curmax!=0) maximums <- c(maximums,curmaxindex)
        }
      }else{
        if(curmax!=0) maximums <- c(maximums,curmaxindex)
        peak <- FALSE
        curmax <- 0
      }
      marker <- marker+1
    }
    peaks[[x]] <- maximums
  }
  peaks
}