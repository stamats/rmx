###############################################################################
## normal location and scale
###############################################################################
.getA1.norm <- function(r){
    approx(x = .radius.gitter.norm, y = .A1.norm, xout = r, yleft = 1)$y
}
.getA2.norm <- function(r){
    approx(x = .radius.gitter.norm, y = .A2.norm, xout = r, yleft = 0.5)$y
}
.geta.norm <- function(r){
    approx(x = .radius.gitter.norm, y = .a.norm, xout = r, yleft = 0)$y
}
.getb.norm <- function(r){
    approx(x = .radius.gitter.norm, y = .b.norm, xout = r, yleft = Inf)$y
}

