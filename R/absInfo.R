absInfo <- function(x, ...){
  UseMethod("absInfo")
}
absInfo.rmx <- function(x, ...){
  IFx <- x$rmxIF$IFun(x$x)
  if(ncol(IFx) == 1){
    info <- IFx^2
  }else{
    info <- rowSums(IFx^2)
  }
  cbind(x = x$x, absInfo = info)
}
