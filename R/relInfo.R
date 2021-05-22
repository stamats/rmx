relInfo <- function(x, ...){
  UseMethod("relInfo")
}
relInfo.rmx <- function(x, ...){
  IFx <- x$rmxIF$IFun(x$x)
  if(ncol(IFx) == 1){
    info <- IFx^2
  }else{
    info <- rowSums(IFx^2)
  }
  res <- cbind(x$x, IFx^2/info)
  colnames(res) <- c("x", paste0("relInfo ", names(x$rmxIF$parameter)))
  res
}
