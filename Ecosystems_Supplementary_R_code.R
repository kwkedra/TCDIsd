
#Supplementary Information
#TITLE: Spectrally Guided Delineation of Structurally 
#Homogeneous Vegetation Patches using Sentinel-2 Imagery 
#in Heterogeneous Post-Industrial Landscapes
#
#Kamil Kedra* and Andrzej M. Jagodzinski
#Institute of Dendrology, Polish Academy of Sciences, Kórnik, 62-035, Poland
#*k.w.kedra@gmail.com


#######################
# 1. R packages       #
#######################

if(!require("raster")) install.packages("raster")
if(!require("stars")) install.packages("stars")
if(!require("terra")) install.packages("terra")
if(!require("sf")) install.packages("sf")
if(!require("lidR")) install.packages("lidR")
if(!require("car")) install.packages("car")

#######################
# 2. Load functions   #
#######################

##calc.TCDIsd
##Description: Calculate a TCDIsd raster
##Arguments:
##H - Sentinel-2 raster::stack object with 12 reflectance bands
calc.TCDIsd <- function(H) {
#citation: 
#Crist EP. 1985. A TM Tasseled Cap equivalent transformation 
#  for reflectance factor data. Remote Sensing of Environment 
#  17(3): 301–306. https://doi.org/10.1016/0034-4257(85)90102-6
#Healey S, Cohen W, Zhiqiang Y, Krankina O. 2005. Comparison of 
#  Tasseled Cap-based Landsat data structures for use in forest 
#  disturbance detection. Remote Sensing of Environment 
#  97(3): 301–310. https://doi.org/10.1016/j.rse.2005.05.009
TCB <- 0.2043*H[[2]] + 0.4158*H[[3]] + 0.5524*H[[4]] + 0.5741*H[[8]] + 0.3124*H[[11]] + 0.2303*H[[12]]
TCG <- -0.1603*H[[2]] - 0.2819*H[[3]] - 0.4934*H[[4]] + 0.7940*H[[8]] - 0.0002*H[[11]] - 0.1446*H[[12]]
TCW <- 0.0315*H[[2]] + 0.2021*H[[3]] + 0.3102*H[[4]] + 0.1594*H[[8]] - 0.6806*H[[11]] - 0.6109*H[[12]]
TCDI <- TCB - (TCG + TCW)
m <- matrix(1, nc=5, nr=5)
b <- raster::focal(x=TCDI, w=m, fun=sd)
return(b)
}

##calc.SeLImin
##Description: Calculate a SeLImin raster
##Arguments:
##H - Sentinel-2 raster::stack object with 12 reflectance bands
calc.SeLImin <- function(H) {
#citation: 
#Pasqualotto N, Delegido J, Van Wittenberghe S, Rinaldi M, 
#  Moreno J. 2019. Multi-Crop Green LAI Estimation with a New 
#  Simple Sentinel-2 LAI Index (SeLI). Sensors, 
#  19(4): 904. https://doi.org/10.3390/s19040904
SeLI <- (H[[9]] - H[[5]]) / (H[[9]] + H[[5]])
m <- matrix(1, nc=5, nr=5)
d <- raster::focal(x=SeLI, w=m, fun=min)
return(d)
}

##getPolygons
##Description: Compute a set of polygons for a desired level 
##of homogeneity
##Arguments:
##rast1 - TCDIsd raster
##rast2 - SeLImin raster
##thres1 - TCDIsd threshold (upper limit)
##thres2 - SeLImin threshold (lower limit)
##col1, col2, lwd., lty. - graphical parameters
##plot. - logical, should the polygons be added to a plot?
##temp.wd - working directory for storing temporary files (required)
getPolygons <- function(rast1, rast2, thres1, thres2, col1=1, col2=0, lwd.=1, lty.=1, plot.=TRUE, temp.wd) {
  TCDIsdVal <- raster::getValues(rast1)
  SeLIminVal <- raster::getValues(rast2)
  idx <- which( TCDIsdVal < thres1 & SeLIminVal > thres2 )
if(length(idx) > 0) {
  rast.new <- rast1
  idxna <- which(is.na(raster::values(rast.new)))
  raster::values(rast.new) <- 0
    #plot(rast.new)
  raster::values(rast.new)[idx] <- 1
  raster::values(rast.new)[idxna] <- NA
    #raster::plot(rast.new)
raster::writeRaster(rast.new, paste0(temp.wd, "H_temp.tif"), overwrite=TRUE )
  rast.new2 <- stars::read_stars( paste0(temp.wd, "H_temp.tif") )
  rast.new2p = sf::st_as_sf(rast.new2, as_points=FALSE, merge=TRUE)
  rast.new2p
  idx_1  <- which(rast.new2p[[1]] == 1)
  rast.new2p_1  <- rast.new2p$geometry[idx_1]
  rast.new2p_1
area_m2 <- as.numeric(sf::st_area(rast.new2p_1))
area_m2
  if( length(which(area_m2==0)) > 0 ) {
  idx.null <- which(area_m2==0)
  rast.new2p_1 <- rast.new2p_1[-idx.null,]
  area_m2 <- as.numeric(sf::st_area(rast.new2p_1))
  }
if( length(rast.new2p_1) > 0 ) {
if(plot.) { plot(rast.new2p_1, add=TRUE, border=col1, col=col2, lwd=lwd., lty=lty. ) }
return(rast.new2p_1)
}else{ return(NULL) }
}else{ return(NULL) }
}

##check.in
##Description: Check if a point falls within homogeneity class polygons (logical output)
##Arguments:
##core - a set of polygons (sf object)
##pt.sel - a selected point (sf object)
check.in <- function(core, pt.sel) {
options(warn=-1)
if(!is.null(core)) {
core_int <- sf::st_intersects(pt.sel, core)
core_int_val <- length(core_int[[1]]) > 0
}else{ core_int_val <- FALSE }
return(core_int_val)
options(warn=0)
}

##LiDAR_var.test
##Description: Perform the Levene's test for homogeneity of variances, 
##between four groups of LiDAR height standard deviations
##Arguments:
##Ldr_data - a 40 by 40 m portion of a LiDAR point cloud (lidR object)
##i - index (number) of a selected point (or a field plot)
##b.col, f.col - graphical parameters (border color and font color, respectively) to mark a LiDAR plot on a map
##wd - working directory for saving tif files: maps and boxplots for LiDAR height standard deviations
##Output: Levene's test p-value
LiDAR_var.test <- function(Ldr_data, i, b.col=0, f.col=0, wd=NA) {
polygon( c(min(Lid$X),min(Lid$X),max(Lid$X),max(Lid$X)) , c(min(Lid$Y),max(Lid$Y),max(Lid$Y),min(Lid$Y)) ,lty=1, lwd=2, border=b.col  ) 
text(mean(Lid$X), mean(Lid$Y), labels=i, col=f.col)
  idx.ground <- which( Ldr_data$Classification == 2 & Ldr_data$Z < 0.5 )
  if(length(idx.ground) > 0) { Ldr_data <- Ldr_data[-idx.ground,] }
  f.sd  <- ~list(sdz = sd(Z))
Ldr_data$X <- Ldr_data$X - min(Ldr_data$X)
Ldr_data$Y <- Ldr_data$Y - min(Ldr_data$Y)
Ldr_data$Z <- Ldr_data$Z - min(Ldr_data$Z)
  idx.p1 <- which(Ldr_data$X < 20  & Ldr_data$Y >= 20)
  idx.p2 <- which(Ldr_data$X >= 20 & Ldr_data$Y >= 20)
  idx.p3 <- which(Ldr_data$X < 20  & Ldr_data$Y < 20)
  idx.p4 <- which(Ldr_data$X >= 20 & Ldr_data$Y < 20)
p1 <- Ldr_data[idx.p1,]
p2 <- Ldr_data[idx.p2,]
p3 <- Ldr_data[idx.p3,]
p4 <- Ldr_data[idx.p4,]
  p1r <- raster::raster(resolution=4, ext=raster::extent(c(0, 20,20,40)), crs=lidR::crs(p1) )
  p2r <- raster::raster(resolution=4, ext=raster::extent(c(20,40,20,40)), crs=lidR::crs(p2) )
  p3r <- raster::raster(resolution=4, ext=raster::extent(c(0, 20,0, 20)), crs=lidR::crs(p3) )
  p4r <- raster::raster(resolution=4, ext=raster::extent(c(20,40,0, 20)), crs=lidR::crs(p4) )
p1_pm.sd <- raster::rasterize( x=cbind(p1$X,p1$Y), field=p1$Z , y=p1r, fun=function(Z, ...) c(length(Z),mean(Z),sd(Z)) )
p2_pm.sd <- raster::rasterize( x=cbind(p2$X,p2$Y), field=p2$Z , y=p2r, fun=function(Z, ...) c(length(Z),mean(Z),sd(Z)) )
p3_pm.sd <- raster::rasterize( x=cbind(p3$X,p3$Y), field=p3$Z , y=p3r, fun=function(Z, ...) c(length(Z),mean(Z),sd(Z)) )
p4_pm.sd <- raster::rasterize( x=cbind(p4$X,p4$Y), field=p4$Z , y=p4r, fun=function(Z, ...) c(length(Z),mean(Z),sd(Z)) )
    if(!is.na(wd)) {
    devcur <- dev.cur()
    main. <- paste0("point_", i ) 
tiff(filename= paste0(wd, "raster_", main., ".tif"), width= 665, height= 700, units= "px")
par(mfrow=c(2,2), mar=c(5,2,4.5,6)+1); raster::plot(p1_pm.sd[[3]], col=terrain.colors(255), main="p1"); raster::plot(p2_pm.sd[[3]], col=terrain.colors(255), main="p2"); raster::plot(p3_pm.sd[[3]], col=terrain.colors(255), main="p3"); raster::plot(p4_pm.sd[[3]], col=terrain.colors(255), main="p4")
    dev.off()
    dev.set(devcur)
    }
  p1_zsd <- raster::values(p1_pm.sd[[3]])
  p2_zsd <- raster::values(p2_pm.sd[[3]])
  p3_zsd <- raster::values(p3_pm.sd[[3]])
  p4_zsd <- raster::values(p4_pm.sd[[3]])
iq1 <- which(p1_zsd > quantile(p1_zsd,0.025,na.rm=TRUE) & p1_zsd < quantile(p1_zsd,0.975,na.rm=TRUE) )
iq2 <- which(p2_zsd > quantile(p2_zsd,0.025,na.rm=TRUE) & p2_zsd < quantile(p2_zsd,0.975,na.rm=TRUE) )
iq3 <- which(p3_zsd > quantile(p3_zsd,0.025,na.rm=TRUE) & p3_zsd < quantile(p3_zsd,0.975,na.rm=TRUE) )
iq4 <- which(p4_zsd > quantile(p4_zsd,0.025,na.rm=TRUE) & p4_zsd < quantile(p4_zsd,0.975,na.rm=TRUE) )
  npts <- c( length(p1_zsd), length(p2_zsd), length(p3_zsd), length(p4_zsd) )
  npts.iq <- c( length(iq1), length(iq2), length(iq3), length(iq4) )
  var. <- c( var(p1_zsd[iq1]), var(p2_zsd[iq2]), var(p3_zsd[iq3]), var(p4_zsd[iq4]) )
  sd.  <- c( sd(p1_zsd[iq1]), sd(p2_zsd[iq2]), sd(p3_zsd[iq3]), sd(p4_zsd[iq4]) )
  mea. <- c( mean(p1_zsd[iq1]), mean(p2_zsd[iq2]), mean(p3_zsd[iq3]), mean(p4_zsd[iq4]) )
print( round(rbind(npts, npts.iq, var., sd., mea.),2) ); flush.console()
if( min(npts.iq) > 15 ) {
  #hist(p1_zsd[iq1])
  #hist(p2_zsd[iq2])
  #hist(p3_zsd[iq3])
  #hist(p4_zsd[iq4])
dat.z <- rbind(
cbind(z=p1_zsd[iq1], gr="idx.p1"),
cbind(z=p2_zsd[iq2], gr="idx.p2"),
cbind(z=p3_zsd[iq3], gr="idx.p3"),
cbind(z=p4_zsd[iq4], gr="idx.p4")
)
dat.z <- as.data.frame(dat.z)
dat.z$z <- as.numeric(dat.z$z) 
dat.z$gr <- as.factor(dat.z$gr) 
head(dat.z)
  LT <- car::leveneTest(z ~ gr, data=dat.z)
  LTdf <- as.data.frame(LT)[1,]
  colnames(LTdf) <- c("Df", "F.value", "p.value")
  rownames(LTdf) <- 1
  LTp <- round(LT$Pr[1], 4)
    if(!is.na(wd)) {
    devcur <- dev.cur()
    main. <- paste0("point_", i ) 
    tiff(filename= paste0(wd, "boxplot_", main.,".tif"), width= 300, height= 500, units= "px")
    boxplot(p1_zsd[iq1], p2_zsd[iq2], p3_zsd[iq3], p4_zsd[iq4], col=c(8,2,3,4), names=c("p1","p2","p3","p4") )
    title(main=main., sub=paste0("p-value = ", LTp) )
    dev.off()
    dev.set(devcur)
    }
return( data.frame(LTdf) )
}else{ return( data.frame(Df=NA, F.value=NA, p.value=NA) ) }
}


#######################
# 3. Usage            #
#######################

##Load a Sentinel-2 image, with 12 reflectance bands (unified10-m resolution),
##in a projected coordinate reference system (units in meters), 
##cropped and masked to a desired region of interest with a sufficient buffer
H <- raster::stack("C:/.../...bsq")

##Calculate the TCDIsd raster
TCDIsd <- calc.TCDIsd(H)

##Plot the TCDIsd raster, other layers will be added to this plot
raster::plot(TCDIsd, col=terrain.colors(255))

##Calculate the SeLImin raster
SeLImin <- calc.SeLImin(H)
#raster::plot(SeLImin)

##Define a working directory folder (it must be created on a disk first)
temp. <- "C:/.../"

##Compute polygons for a desired homogeneity level; input rasters without a buffer
core010 <- getPolygons(rast1=TCDIsd, rast2=SeLImin, thres1=0.010, thres2=0.02,
           col1="darkblue", col2=rgb(0,0,1,.1), lwd.=3, lty.=1, plot.=TRUE, 
           temp.wd=temp.)

##Calculate areas of the polygons (in m2)
a <- sf::st_area(core010)
a

##Define an index to select one polygon (for demonstration)
i <- 1
##Create a centroid of the selected polygon
pt.sel <- sf::st_centroid(core010[i])

##or iterate over a set of points
#  pts.sel <- sf::st_centroid(core010)
#for(i in 1:length(pts.sel)) {
#pt.sel <- pts.sel[i]

##Check if the selected point falls within the selected polygons
check.in010 <- check.in(core010, pt.sel=pt.sel)
print(check.in010); flush.console()

##Get the coordinates of a point
xyi <- sf::st_coordinates(pt.sel)

##Define a filter to load only the required portion of a LiDAR point cloud, define redundant classes to skip
fltr <- paste("-keep_tile", round(xyi[1],3)-20, round(xyi[2],3)-20, 40, "-drop_class 1 6")

##Set working directory to the folder with LiDAR data
setwd("C:/.../")
list.files()

##Load the LiDAR data (in the same crs as the Sentinel-2 image) and normalize height
Lid <- lidR::readLAS("….las", filter=fltr)
Lidn <- lidR::normalize_height(Lid, lidR::tin())
#lidR::plot(Lidn)

##Perform the Levene's test (check if a plot is homogeneous in terms of LiDAR data)
LVT <- LiDAR_var.test(Lidn, i=i, b.col=2, f.col=0, wd=temp.)
LVT
#}

#dev.off()