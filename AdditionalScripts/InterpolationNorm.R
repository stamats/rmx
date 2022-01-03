###############################################################################
## Interpolated functions to speed up computation of Lagrange Multipliers
###############################################################################

## source optIFNorm.R
radius <- c(1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, seq(1e-4, 0.01, by = 0.001),
            seq(0.02, 5, by = 0.01), seq(5.05, 10, by = 0.05))

fun <- function(radius){
  print(radius)
  optIF.norm(radius, delta = 1e-9)
}
locationScale <- sapply(radius, fun)

## location and scale
n <- length(radius)

.A1.norm <- unlist(locationScale[4,])[seq(1, 4*n-3, by = 4)]
.A2.norm <- unlist(locationScale[4,])[seq(4, 4*n, by = 4)]
.a.norm <- unlist(locationScale[5,])[seq(2, 2*n, by = 2)]
.b.norm <- unlist(locationScale[6,])
.radius.gitter.norm <- radius
.asVar.mean.norm <- sapply(locationScale[10,], function(x) x[1,1])
.asVar.sd.norm <- sapply(locationScale[10,], function(x) x[2,2])


plot(radius, .A1.norm, type = "l")
plot(radius, .A2.norm, type = "l")
plot(radius, .a.norm, type = "l")
plot(radius, .b.norm, type = "l")
plot(radius, .asVar.mean.norm, type = "l")
plot(radius, .asVar.sd.norm, type = "l")

radius[radius > 1.62 & radius < 1.8]
print(.b.norm[radius > 1.62 & radius < 1.8], digits = 10)

## Saving the results in sysdata.rda
#load("sysdata.rda")
save(.radius.gitter.norm, .A1.norm, .A2.norm, .a.norm, .b.norm, 
     .asVar.mean.norm, .asVar.sd.norm, file = "sysdata.rda")

## and proceed with FiniteSampleCorrectionFactorNorm.R
