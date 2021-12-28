absInfo <- function(object, ...){
  UseMethod("absInfo")
}
absInfo.rmx <- function(object, ...){
  IFx <- object$rmxIF$IFun(object$x)
  if(ncol(IFx) == 1){
    info <- IFx^2
  }else{
    info <- rowSums(IFx^2)
  }
  cbind(x = object$x, absInfo = info)
}
