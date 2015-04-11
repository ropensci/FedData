# FedData Tester
devtools::install_github("bocinsky/FedData")
library(FedData)
library(raster)
library(png)
library(RColorBrewer)

setwd("~/Desktop/FedData Test")

# Get a random contiguous USA county for testing
wgetDownload("http://dds.cr.usgs.gov/pub/data/nationalatlas/countyp010g.shp_nt00934.tar.gz",destdir=getwd())
untar("./countyp010g.shp_nt00934.tar.gz")
county <- rgdal::readOGR(".","countyp010g")
county <- county[!(county$STATE %in% c("AK","VI","PR","HI")),]
county <- county[sample(1:length(county),1),]
# county <- county[which(county$NAME=='Napa'),]

# Get the NED (USA ONLY)
# Returns a raster
NED <- getNED(template=county, label=paste(county$STATE,'_',county$NAME, sep=''), raw.dir="/Users/Bocinsky/Desktop/FedData Test/RAW/NED/",extraction.dir="/Users/Bocinsky/Desktop/FedData Test/EXTRACTIONS/NED/", res='1')

# Get the daily GHCN data (GLOBAL)
# Returns a list: the first element is the spatial locations of stations,
# and the second is a list of the stations and their daily data
GHCN.prcp <- getGHCNDaily(template=county, label=paste(county$STATE,'_',county$NAME, sep=''), elements=c('prcp'), standardize=F)
GHCN.temp <- getGHCNDaily(template=county, label=paste(county$STATE,'_',county$NAME, sep=''), elements=c('tmin','tmax'), standardize=T)

# Get the NHD (USA ONLY)
NHD <- getNHD(template=county, label=paste(county$STATE,'_',county$NAME, sep=''))

# Get the NRCS SSURGO data (USA ONLY)
SSURGO <- getSSURGO(template=county, label=paste(county$STATE,'_',county$NAME, sep=''))

slope <- terrain(NED, opt='slope')
aspect <- terrain(NED, opt='aspect')
NED.hill <- hillShade(slope, aspect, 40, 230)

NED <- raster::mask(NED,county)
NED.hill <- raster::mask(NED.hill,county)

plot(NED)
plot(county, add=T)
plot(NHD$NHDFlowline, col='gray50', border='gray50', add=T)
plot(NHD$NHDWaterbody, col='gray50', border='gray50', add=T)
plot(NHD$NHDArea, col='gray50', border='gray50', add=T)
plot(county, add=T)
plot(GHCN.prcp[[1]], pch=17, add=T)
plot(GHCN.temp[[1]], pch=19, add=T)

# plot(SSURGO[['spatial']])

wine <- merge(SSURGO[['spatial']],SSURGO[['tabular']]$mucropyld[SSURGO[['tabular']]$mucropyld$cropname=='Wine grapes',],by.x='MUKEY',by.y="mukey",all=T)[,c("nonirryield.r","irryield.r")]
wine$yield.max <- pmax(wine$irryield.r,wine$nonirryield.r,na.rm=T)




# Plot NED
colors <- paste0(colorRampPalette(brewer.pal(9, "Greens"),bias=2)(1000),"CC")
plot(county)
plot(NED.hill, col=grey(30:100/100), maxpixels=2000000, useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
plot(NED, col=colors, maxpixels=2000000, useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
plot(county, add=T)

# Plot NHD
plot(NHD$NHDFlowline)
sp::plot(NHD$NHDWaterbody, col='black', add=T)
sp::plot(NHD$NHDArea, col='black', add=T)

# Plot average potential wine grape yield
colors <- paste0(colorRampPalette(brewer.pal(9, "Reds"))(1000),"CC")
wine$colors <- colors[wine$yield.max*100]
plot(county)
plot(NED.hill, col=grey(70:100/100), maxpixels=2000000, useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
sp::plot(wine, col=wine$colors,border=NA, add=T)
plot(county, add=T)





itrdb.data <- getITRDB(recon.years=1:2000, 
                       calib.years=1924:1983, 
                       species=NULL, 
                       measurement.type="Ring Width", 
                       chronology.type="Standard", 
                       makeSpatial=T)







fig.ratio <- 0.64
plot.height <- 4.3
plot.width <- plot.height*fig.ratio
legend.height <- 0.5
edge.margin <- 0.1
horizontal.between <- (10-(plot.width*3)-(edge.margin*2))/2

fig.width <- (plot.width*3)+(edge.margin*2)+(horizontal.between*2)
fig.height <- (plot.width*(1/fig.ratio))+legend.height+(edge.margin*2)
fig.ratio <- fig.height/fig.width



plot.map <- function(the.rasters, background, label, colors, color.breaks, range.z, ...){
  plot(1, type='n', xlab="", ylab="", xlim=c(xmin(the.rasters[[1]]),xmax(the.rasters[[1]])),ylim=c(ymin(the.rasters[[1]]),ymax(the.rasters[[1]])), xaxs="i", yaxs="i", axes=FALSE, main='')
  
  quartz(file='../FIGURES/temp.png', width=diff(c(xmin(the.rasters[[1]]),xmax(the.rasters[[1]]))), height=diff(c(ymin(the.rasters[[1]]),ymax(the.rasters[[1]]))), antialias=FALSE, bg="white", type='png', family="Gulim", pointsize=8, dpi=1200)
  par(mai=c(0,0,0,0))
  plot(1, type='n', xlab="", ylab="", xlim=c(xmin(the.rasters[[1]]),xmax(the.rasters[[1]])),ylim=c(ymin(the.rasters[[1]]),ymax(the.rasters[[1]])), xaxs="i", yaxs="i", axes=FALSE, main='')
  plot(background, col=grey(30:100/100), maxpixels=500000, useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
  plot(the.rasters[[1]], zlim=range.z, breaks=color.breaks, col=colors[[1]], useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
  if(length(the.rasters)>1){
    for(i in 2:length(the.rasters)){
      plot(the.rasters[[i]], zlim=range.z, breaks=color.breaks, col=colors[[i]], useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
    }
  }
  dev.off()
  img <- readPNG('../FIGURES/temp.png')
  #   system("rm ../FIGURES/temp.png")
  rasterImage(img, xmin(the.rasters[[1]]),ymin(the.rasters[[1]]),xmax(the.rasters[[1]]),ymax(the.rasters[[1]]), interpolate=F)
  
  #   points(y=35.700275, x=-105.908613, pch=17, cex=2)
  #   text(y=35.700275, x=-105.908613, labels="SFI", font=2, cex=3, adj=c(1.1,-0.1))
  
  inch <- (extent(the.rasters[[1]])@xmax-extent(the.rasters[[1]])@xmin)/(fig.width-par('mai')[2]-par('mai')[4])
  scalebar.new(d=10, cex=2, font=2, side='left',lab.side='left', height=0.075*inch, label="", line.offset=c(-0.05*inch,0.05*inch), xy=c(xmax(the.rasters[[1]]),ymin(the.rasters[[1]])), lwd=4, lend=1)
  text(x=xmin(the.rasters[[1]])+(0.075*inch), y=ymin(the.rasters[[1]])+(0.075*inch), labels=label, col="black", font=2, cex=3, adj=c(0,0), xpd=T)  
}





cairo_pdf(filename=paste('./NAPA.pdf',sep=''), width=fig.width, height=fig.height, bg="white", pointsize=8)

## NED

par(mai=c(edge.margin+(plot.width*fig.ratio*0)+(plot.width*fig.ratio*0)+(legend.height*1),
          edge.margin,
          edge.margin,
          edge.margin+(plot.width*2)+(horizontal.between*2)), 
    oma=c(0,0,0,0), 
    lend=2, 
    ljoin=1, 
    xpd=F)

colors <- paste0(colorRampPalette(brewer.pal(9, "Greens"),bias=2)(1000),"CC")
plot(1, type='n', xlab="", ylab="", xlim=c(extent(county)@xmin,extent(county)@xmax),ylim=c(extent(county)@ymin,extent(county)@ymax), xaxs="i", yaxs="i", axes=FALSE, main='')
quartz(file='./temp.png', width=diff(c(extent(county)@xmin,extent(county)@xmax)), height=diff(c(extent(county)@ymin,extent(county)@ymax)), antialias=FALSE, bg="white", type='png', family="Gulim", pointsize=8, dpi=1200)
par(mai=c(0,0,0,0))
plot(1, type='n', xlab="", ylab="", xlim=c(extent(county)@xmin,extent(county)@xmax),ylim=c(extent(county)@ymin,extent(county)@ymax), xaxs="i", yaxs="i", axes=FALSE, main='')
plot(NED.hill, col=grey(30:100/100), maxpixels=2000000, useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
plot(NED, col=colors, maxpixels=2000000, useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
dev.off()
img <- readPNG('./temp.png')
#   system("rm ../FIGURES/temp.png")
rasterImage(img, extent(county)@xmin,extent(county)@ymin,extent(county)@xmax,extent(county)@ymax, interpolate=F)
plot(county, ljoin=0, lend=0, add=T)

inch <- (extent(county)@xmax-extent(county)@xmin)/(fig.width-par('mai')[2]-par('mai')[4])
scalebar.new(d=10, cex=1.25, font=2, side='right',lab.side='right', height=0.05*inch, label="10 km", line.offset=c(0*inch,0*inch), xy=c(extent(county)@xmin,extent(county)@ymin), lwd=4, lend=1)
text(x=extent(county)@xmin, y=extent(county)@ymax, labels='a', col="black", font=2, cex=3, adj=c(0,1), xpd=T)  

# Legend
par(mai=c(edge.margin+(plot.width*fig.ratio*0),
          edge.margin,
          edge.margin+(plot.height*1),
          edge.margin+(horizontal.between*2)+(plot.width*2)), 
    oma=c(0,0,0,0), 
    lend=2, 
    ljoin=1, 
    xpd=F, 
    new=T)
plot(1, type='n', xlab="", ylab="", xlim=c(0,1000),ylim=c(0,legend.height), font=2, xaxs="i", yaxs="i", axes=FALSE, main='')

quartz(file='./temp.png', width=plot.width, height=1, antialias=FALSE, bg="white", type='png', family="Gulim", pointsize=8, dpi=1200)
par(mai=c(0,0,0,0))
plot(1, type='n', xlab="", ylab="", xlim=c(0,1000),ylim=c(0,legend.height), xaxs="i", yaxs="i", axes=FALSE, main='')
rect(col=colors, border=NA, xleft=0:999, xright=1:1000, ybottom=0, ytop=0.2, xpd=T)
# rect(xleft=0, xright=1, ybottom=c(0,cumsum(VEPIIS.grow.GDD.retro.union.hist.matrix.mean)[1:999]), ytop=c(cumsum(VEPIIS.grow.GDD.retro.union.hist.matrix.mean),1), col=rev(GDD.colors), border=NA)
dev.off()
img <- readPNG('./temp.png')
# system("rm ../FIGURES/temp.png")
rasterImage(img, 0,0,1000,legend.height, interpolate=F)

# rect(col=NA,border=NA, xleft=H2O_YEAR.brk[1:length(H2O_YEAR.colors)], xright=H2O_YEAR.brk[2:(length(H2O_YEAR.colors)+1)], ybottom=0.6, ytop=0.85, xpd=T)
# segments(x0=30,x1=30,y0=0.6, y1=0.85, lwd=2, lend=1)
text(x=0, y=0.35, labels=0, cex=1.5, font=2, adj=c(0,0.5), xpd=T)  
text(x=1000, y=0.35, labels=1270, cex=1.5, font=2, adj=c(1,0.5), xpd=T)
text(x=500, y=0.35, labels="elevation (masl)", font=2, cex=1.5, adj=c(0.5,0.5), xpd=T)  


## NHD
par(mai=c(edge.margin+(plot.height*0)+(legend.height*1),
          edge.margin+(plot.width*1)+(horizontal.between*1),
          edge.margin,
          edge.margin+(plot.width*1)+(horizontal.between*1)), 
    oma=c(0,0,0,0), 
    lend=2, 
    ljoin=1, 
    xpd=F, 
    new=T)

plot(1, type='n', xlab="", ylab="", xlim=c(extent(county)@xmin,extent(county)@xmax),ylim=c(extent(county)@ymin,extent(county)@ymax), xaxs="i", yaxs="i", axes=FALSE, main='')
plot(NHD$NHDFlowline, ljoin=0, lend=0, add=T)
sp::plot(NHD$NHDWaterbody, col='black', ljoin=0, lend=0, add=T)
sp::plot(NHD$NHDArea, col='black', ljoin=0, lend=0, add=T)
# plot(county, add=T)

inch <- (extent(county)@xmax-extent(county)@xmin)/(fig.width-par('mai')[2]-par('mai')[4])
scalebar.new(d=10, cex=1.25, font=2, side='right',lab.side='right', height=0.05*inch, label="10 km", line.offset=c(0*inch,0*inch), xy=c(extent(county)@xmin,extent(county)@ymin), lwd=4, lend=1)
text(x=extent(county)@xmin, y=extent(county)@ymax, labels='b', col="black", font=2, cex=3, adj=c(0,1), xpd=T)  

# Legend
par(mai=c(edge.margin,
          edge.margin+(horizontal.between*1)+(plot.width*1),
          edge.margin+(plot.height*1),
          edge.margin+(horizontal.between*1)+(plot.width*1)), 
    oma=c(0,0,0,0), 
    lend=2, 
    ljoin=1, 
    xpd=F, 
    new=T)
plot(1, type='n', xlab="", ylab="", xlim=c(0,1000),ylim=c(0,legend.height), font=2, xaxs="i", yaxs="i", axes=FALSE, main='')
legend(x=500,y=0.25,legend="flowlines/areas",col='black',lty=1,lwd=2,xjust=0.5,yjust=0.5,cex=1.5,bty='n',text.font=2, xpd=T)
# 
# 
# # rect(col=NA,border=NA, xleft=H2O_YEAR.brk[1:length(H2O_YEAR.colors)], xright=H2O_YEAR.brk[2:(length(H2O_YEAR.colors)+1)], ybottom=0.6, ytop=0.85, xpd=T)
# # segments(x0=30,x1=30,y0=0.6, y1=0.85, lwd=2, lend=1)
# text(x=0, y=0.35, labels=0, cex=1.6, font=2, adj=c(0,0.5), xpd=T)  
# text(x=1000, y=0.35, labels=1270, cex=1.6, font=2, adj=c(1,0.5), xpd=T)
# text(x=500, y=0.35, labels="Elevation (masl)", font=2, cex=1.6, adj=c(0.5,0.5), xpd=T)  


## SSURGO
par(mai=c(edge.margin+(plot.height*0)+(legend.height*1),
          edge.margin+(plot.width*2)+(horizontal.between*2),
          edge.margin,
          edge.margin+(plot.width*0)+(horizontal.between*0)), 
    oma=c(0,0,0,0), 
    lend=2, 
    ljoin=1, 
    xpd=F, 
    new=T)
plot(1, type='n', xlab="", ylab="", xlim=c(extent(county)@xmin,extent(county)@xmax),ylim=c(extent(county)@ymin,extent(county)@ymax), xaxs="i", yaxs="i", axes=FALSE, main='')
colors <- paste0(colorRampPalette(brewer.pal(9, "Reds"))(1000),"CC")
wine$colors <- colors[wine$yield.max*100]
plot(NED.hill, col=grey(70:100/100), maxpixels=2000000, useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
sp::plot(wine, col=wine$colors,border=NA, ljoin=0, lend=0, add=T)
plot(county, ljoin=0, lend=0, add=T)

inch <- (extent(county)@xmax-extent(county)@xmin)/(fig.width-par('mai')[2]-par('mai')[4])
scalebar.new(d=10, cex=1.25, font=2, side='right',lab.side='right', height=0.05*inch, label="10 km", line.offset=c(0*inch,0*inch), xy=c(extent(county)@xmin,extent(county)@ymin), lwd=4, lend=1)
text(x=extent(county)@xmin, y=extent(county)@ymax, labels='c', col="black", font=2, cex=3, adj=c(0,1), xpd=T)  

# Legend
par(mai=c(edge.margin,
          edge.margin+(horizontal.between*2)+(plot.width*2),
          edge.margin+(plot.height*1),
          edge.margin), 
    oma=c(0,0,0,0), 
    lend=2, 
    ljoin=1, 
    xpd=F, 
    new=T)
plot(1, type='n', xlab="", ylab="", xlim=c(0,1000),ylim=c(0,legend.height), font=2, xaxs="i", yaxs="i", axes=FALSE, main='')

quartz(file='./temp.png', width=plot.width, height=1, antialias=FALSE, bg="white", type='png', family="Gulim", pointsize=8, dpi=1200)
par(mai=c(0,0,0,0))
plot(1, type='n', xlab="", ylab="", xlim=c(0,1000),ylim=c(0,legend.height), xaxs="i", yaxs="i", axes=FALSE, main='')
rect(col=colors, border=NA, xleft=0:999, xright=1:1000, ybottom=0, ytop=0.2, xpd=T)
# rect(xleft=0, xright=1, ybottom=c(0,cumsum(VEPIIS.grow.GDD.retro.union.hist.matrix.mean)[1:999]), ytop=c(cumsum(VEPIIS.grow.GDD.retro.union.hist.matrix.mean),1), col=rev(GDD.colors), border=NA)
dev.off()
img <- readPNG('./temp.png')
# system("rm ../FIGURES/temp.png")
rasterImage(img, 0,0,1000,legend.height, interpolate=F)

# rect(col=NA,border=NA, xleft=H2O_YEAR.brk[1:length(H2O_YEAR.colors)], xright=H2O_YEAR.brk[2:(length(H2O_YEAR.colors)+1)], ybottom=0.6, ytop=0.85, xpd=T)
# segments(x0=30,x1=30,y0=0.6, y1=0.85, lwd=2, lend=1)
text(x=0, y=0.35, labels=1.5, cex=1.5, font=2, adj=c(0,0.5), xpd=T)  
text(x=1000, y=0.35, labels=8, cex=1.5, font=2, adj=c(1,0.5), xpd=T)
text(x=500, y=0.35, labels="wine grape yield (tons/acre)", font=2, cex=1.5, adj=c(0.5,0.5), xpd=T)  

dev.off()

system("gs -sDEVICE=pdfwrite -dNOPAUSE -dQUIET -dBATCH -dCompatibilityLevel=1.9 -dEmbedAllFonts=true -sOutputFile=./NAPA_FINAL.pdf ./NAPA.pdf")
system("rm ./NAPA.pdf")
system("mv ./NAPA_FINAL.pdf ./NAPA.pdf")






scalebar.new <- function (d, xy = NULL, height = NULL, line.offset=c(0,0), side="right", lab.side='top', lonlat = NULL, label, adj = c(0.5, -0.5), lwd = 2, ...){
  pr <- par()
  if (is.null(lonlat)) {
    if (pr$usr[1] > -181 & pr$usr[2] < 181 & pr$yaxp[1] > 
          -200 & pr$yaxp[2] < 200) {
      lonlat <- TRUE
    }
    else {
      lonlat <- FALSE
    }
  }
  
  if (lonlat) {
    lat <- mean(pr$yaxp[1:2])
    if (missing(d)) {
      dx <- (pr$usr[2] - pr$usr[1])/10
      d <- pointDistance(cbind(0, lat), cbind(dx, lat), 
                         TRUE)
      d <- signif(d/1000, 2)
      label <- NULL
    }
    p <- cbind(0, lat)
    dd <- raster:::.destPoint(p, d * 1000)
    dd <- dd[1, 1]
  } else {
    if (missing(d)) {
      d <- round(10 * (pr$usr[2] - pr$usr[1])/10)/10
      label <- NULL
    }
    dd <- d
  }
  
  if (is.null(xy)) {
    padding = c(5, 5)/100
    parrange <- c(pr$usr[2] - pr$usr[1], pr$usr[4] - pr$usr[3])
    xy <- c(pr$usr[1] + (padding[1] * parrange[1]), pr$usr[3] + 
              (padding[2] * parrange[2]))
  }
  
  xy <- xy + line.offset
  
  if(side=='right'){
    xstart = xy[1]
    xend = xy[1] + dd
  }else{
    xstart = xy[1] - dd
    xend = xy[1]
  }
  
  if(is.null(height)){
    height <- dd * 0.1
  }
  
  rect(xleft=xstart, ybottom=xy[2], xright=xend, ytop=xy[2]+height, col='black',border=NA, lend=1)
  
  #   lines(matrix(c(xstart, xy[2], xend, xy[2]), byrow = T, nrow = 2), lend=1, lwd = lwd, ...)
  
  if (missing(label)) {
    label <- paste(d)
  }
  if (is.null(label)) {
    label <- paste(d)
  }
  if (missing(adj)) {
    adj <- c(0.5, -0.2 - lwd/20)
  }
  
  if(lab.side=='top'){
    text(mean(c(xstart,xend)), xy[2], labels = label, adj = c(0.5,-0.5), 
         ...)
  }else if(lab.side=='right'){
    text(xend, xy[2], labels = label, adj = c(-0.1,0), 
         ...)
  }else if(lab.side=='left'){
    text(xstart, xy[2], labels = label, adj = c(1.1,0), 
         ...)
  }else if(lab.side=='bottom'){
    text(mean(c(xstart,xend)), xy[2], labels = label, adj = c(0.5,1.5), 
         ...)
  }
}

