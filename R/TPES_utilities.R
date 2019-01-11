################################################################################
getSegmentsPos <- function(chr, initialPos, finalPos, segments){

  goodSegs <- which(
    ( as.numeric(gsub("chr","",segments[,1])) == as.numeric(gsub("chr","",chr)) |
        gsub("chr","",segments[,1]) == gsub("chr","",chr)  ) & (
          (as.numeric(segments[,2]) >= (as.numeric(initialPos)) & as.numeric(segments[,2]) <= (as.numeric(finalPos)) ) |
            (as.numeric(segments[,3]) >= (as.numeric(initialPos)) & as.numeric(segments[,3]) <= (as.numeric(finalPos)) ) |
            (as.numeric(segments[,2]) <= (as.numeric(initialPos)) & as.numeric(segments[,3]) >= (as.numeric(finalPos)) )))

  return(goodSegs)}

################################################################################
nr.modes <- function (y){
  d1 <- diff(y)
  signs <- diff(d1/abs(d1))
  length(signs[signs == -2])}

################################################################################
findPeaks <- function (AF, bw) {

  dens <- stats::density(AF, bw = bw, na.rm = T)
  ydens <- dens$y
  xdens <- dens$x
  d1 <- diff(ydens)
  signs <- diff(d1/abs(d1))
  signs2 <- which(signs == -2)

  if(length(signs2) > 0){
    peaks <- data.frame(t(sapply(signs2, function(x, xdens, ydens){
      list(xID=xdens[x], yID=ydens[x])}, xdens, ydens)))
  } else {
    cat(paste0("[",Sys.time() ,"] No peak found\n"))
    peaks <- data.frame()
  }

  return(peaks)}

################################################################################
findMin <- function (AF, bw) {

  dens <- stats::density(AF, bw = bw, na.rm=T)
  ydens <- dens$y
  xdens <- dens$x
  d1 <- diff(ydens)
  signs <- diff(d1/abs(d1))
  signs2 <- which(signs == 2)

  minima <- data.frame(t(sapply(signs2, function(x, xdens){xdens[x]}, xdens)))

  return(minima)}

################################################################################
