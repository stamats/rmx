relInfo <- function(object, ...){
  UseMethod("relInfo")
}
relInfo.rmx <- function(object, ...){
  IFx <- object$rmxIF$IFun(object$x)
  if(ncol(IFx) == 1){
    info <- IFx^2
  }else{
    info <- rowSums(IFx^2)
  }
  res <- cbind(object$x, IFx^2/info)
  colnames(res) <- c("x", paste0("relInfo ", names(object$rmxIF$parameter)))
  res
}
