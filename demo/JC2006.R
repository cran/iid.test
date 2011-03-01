# R.E Benesta Oslo, Norway
# R-script: see URL http://cran.r-project.org
#
# Figures for Benestad, R.E (2006) Can we expect more extreme precipitation on the monthly time scale?  J.Clim  Vol. 19, No. 4, pages 630-637

#rm(list=ls())
library(clim.pact)
library(ncdf)

record.analysis <- function(N.test=1000,plot.every=FALSE,
                            season=NULL,bw=TRUE,force.read=FALSE,
                            dir="GCMs/",scenario="A1b",nx=71,ny=27,nt=99) {

months <- c("J","F","M","A","M","J","J","A","S","O","N","D")
if (!is.null(season)) {
  seas <- "";
  for (i in 1:length(season)) seas <- paste(seas,months[season[i]],sep="")
} else seas <- ""

# Functions called inside record.analysis

make.data.file <- function(dir="GCMs/",scenario="A1b",nx=71,ny=27,nt=99) {

print(" === Record-statistics ===")
prints("Takes GCM results gridded onto a common grid")

files.all <- list.files(pattern="regriddedPrecip",path=dir)
files.all <- files.all[grep(".nc",files.all)]

gcms.a1b <- files.all[grep("SRESA1B",files.all)]
gcms.b1  <- files.all[grep("SRESB1",files.all)]
gcms.a2  <- files.all[grep("SRESA2",files.all)]

gcms <- switch(lower.case(scenario),"a1b"=,gcms.a1b,
               "b1"=,gcms.b1,"a2"=,gcms.a2)
               
#gcms <- gcms.b1
print(gcms)

igcm <- 0
nrecs <- rep(0,nx*ny); dim(nrecs) <- c(nx,ny)
nrecs.all <- nrecs 
JC2006 <- list(nrecs.1=nrecs)
mm <- rep(1:12,nt)
ar1 <- nrecs 

for (gcm in gcms) {
  igcm <- igcm + 1
  print(gcm)
  filename <- paste(dir,gcm,sep="")
  ncid <- open.ncdf(filename)
  data <- get.var.ncdf(ncid,"PRECIP")
  lon <- get.var.ncdf(ncid,"XNS")
  lat <- get.var.ncdf(ncid,"YNS")
  close.ncdf(ncid)
  
  for (i in 1:nx) {
    for (j in 1:ny) {
       nrecs[i,j] <- nrecords(data[i,j,],mm,season)
    } 
  }
  if (plot.every) {
    image(lon,lat,nrecs,main=gcm)
    addland()
    contour(lon,lat,nrecs,add=TRUE)
    if (is.null(season)) fname <- paste("temp/JC2004_",igcm,sep="") else
                         fname <- paste("temp/JC2004_",igcm,seas,sep="")
    dev.print(file=paste(fname,"ps",sep=""),horizontal=FALSE)
    dev2bitmap(file=paste(fname,".jpg",sep=""),type="jpeg"); dev.off()
  }
  
  
  eval(parse(text=paste("JC2006$nrecs.",igcm," <- nrecs",sep="")))
  eval(parse(text=paste("attr(JC2006$nrecs.",igcm,",'GCM') <- gcm",sep="")))
  JC2006$lon <- lon - 360
  JC2006$lat <- lat
  nrecs.all <- nrecs.all + nrecs

  AR1 <-autocor(data,mm,season,seas)
  eval(parse(text=paste("JC2006$a.",igcm," <- AR1",sep="")))
  #print(c(dim(ar1),NA,dim(AR1)))
  ar1 <- ar1 + AR1
}

JC2006$nrecs.all <- nrecs.all/length(gcms)
JC2006$AR1.all <- ar1/length(gcms)

MC.1 <- rep(0,N.test)
MC.all <- rep(0,N.test)
if (is.null(season)) season <- 1:12

#print("Monte-Carlo test")
ii <- is.element(mm,season)
for (i in 1:N.test) {
     vec <- rnorm(nt*length(season))
     MC.1[i] <- nrecords(vec,mm[ii],season)
     for (ie in 1:length(gcms)) {
       vec <- rnorm(nt*length(season))
       MC.all[i] <- MC.all[i] + nrecords(vec,mm[ii],season)
     }
     MC.all[i] <- MC.all[i]/length(gcms)
}

JC2006$MC.1 <- MC.1
JC2006$MC.all <- MC.all
AR1 <- JC2006$AR1.all 
attr(JC2006$MC.1,"method") <-
  "Monte-Carlo integration: 99 x 12 months using stochastic iid record (rnorm)"
attr(JC2006$MC.all,"method") <-
  "Monte-Carlo integration: 10 x 99 x 12 months using stochastic iid record (rnorm)"

save(file="JC2006.rda",JC2006)
}

nrecords <- function(x,mm,season=NULL) {
   print(paste("nrecords: season=",season))
   nrecs <- 0              
   good <- is.finite(x) & is.finite(mm)
   x <- x[good]; mm <- mm[good]
   if (!is.null(season)) {
     good <- is.element(mm,season)
     x <- x[good]; mm <- mm[good]
   } else season <- 1:12
   for (im in season) {
      vec <- x[is.element(mm,im)]
      for (it in 2:length(vec)) {
        if (vec[it] > max(vec[1:(it-1)],na.rm=TRUE)) nrecs <- nrecs + 1
      }
   }
  recs <- nrecs/length(season) + 1    # The first value is always a record by definition!
  invisible(recs)
}

autocor <- function(x,mm,season=NULL,seas="") {
  #print("autocor")
  dims <- dim(x)
  #print(dims)
  nx <- dims[1]; ny <- dims[2]; nt <- dims[3]
  a <- rep(NA,nx*ny); dim(a) <- c(nx,ny)
  dim(x) <- c(nx*ny,nt)
  clim <- rep(NA,nx*ny*12); dim(clim) <- c(nx*ny,12)
  #print(summary(c(x)))
  if (!is.null(season)) months <- season else months <- 1:12
  #print(months)
  for (im in months) {
    ii <- is.element(mm,im)
    clim[,im] <- rowMeans(x[,ii],na.rm=TRUE)
    x[,ii] <- x[,ii] - clim[,im]  # subtract temporal mean for given month
                                  # i.e. remove annual cycle
  }
  good <- is.element(mm,months)
  x <- x[,good]
  #print(summary(c(x)))
  dim(clim) <- c(nx,ny,12)
  dim(x) <- c(nx,ny,sum(good))
  #print(dim(x))
  for (i in 1:nx) {
    for (j in 1:ny) {
       a[i,j] <- acf(c(x[i,j,]),plot=FALSE, na.action =na.pass)$acf[2]
    } 
  }
  invisible(a)
}




if (force.read) make.data.file(dir,scenario,nx,ny,nt) else {
  data(JC2006)
  AR1 <- JC2006$AR1.all
  gcms <- names(JC2006); gcms <- gcms[grep("nrecs.",gcms)]
}




#print(" ------- Saved record-statistics data for paper29/JC2004 -------- ") 

if (bw) {
  my.col <- rgb( c(seq(0.3,0.7,by=0.4)),
                 c(seq(0.3,0.7,by=0.4)),
                 c(seq(0.3,0.7,by=0.4)) )
} else {
  my.col <- rgb( c(rep(1.0,1), rep(1.0,1), rep(0.5,1)),
                 c(rep(0.5,1), rep(1.0,1), rep(0.5,1)),
                 c(rep(0.5,1), rep(1.0,1), rep(1.0,1)) )
}

nrecs.all <- JC2006$nrecs.all
MC.all <- JC2006$MC.all
lon <- JC2006$lon; lat <- JC2006$lat; nx <- length(lon); ny <- length(lat)
nrecs <- nrecs.all
nrecs[(nrecs > quantile(MC.all,0.025)) & (nrecs < quantile(MC.all,0.975))] <- NA

x11()
image(lon,lat,nrecs,main="Number of records",
      sub=paste(length(gcms),"scenarios",seas),col=my.col)
addland(col="grey30",lwd=2)
contour(lon,lat,nrecs.all,add=TRUE,cex=1.6)
if (is.null(season)) fname <- paste("temp/JC2004_all.",sep="") else
                     fname <- paste("temp/JC2004_all_",seas,".",sep="")
dev.print(file=paste(fname,"ps",sep=""),horizontal=FALSE)
dev.copy2eps(file=paste(fname,"eps",sep=""))
dev2bitmap(file=paste(fname,"jpg",sep=""),type="jpeg"); dev.off()

h <- hist(MC.all,breaks=seq(0,10,by=0.01))
edf <- cumsum(h$density)/sum(h$density)
for (i in 1:nx) {
  for (j in 1:ny) {
    nrecs[i,j] <- 1 - min(edf[h$mids >= nrecs.all[i,j]])
  }
} 

x11()
contour(lon,lat,nrecs,main="Number of records: p-values",levels=seq(0,1,by=0.1),
        sub=paste(length(gcms),"scenarios",seas),col="grey")
contour(lon,lat,nrecs,levels=seq(0,0.1,by=0.01),add=TRUE)
addland(col="grey30",lwd=3)
if (is.null(season)) fname <- paste("temp/JC2004_all_p-values.",sep="") else
                     fname <- paste("temp/JC2004_all_p-values_",seas,".",sep="")
dev.print(file=paste(fname,"ps",sep=""),horizontal=FALSE)
dev.copy2eps(file=paste(fname,"eps",sep=""))
dev2bitmap(file=paste(fname,"jpg",sep=""),type="jpeg"); dev.off()

x11()
plot(h$mids,edf,type="l",lwd=2,main="edf: Monte-Carlo results for 10 GCMs",
     xlab="numer of records",ylab="edf")
grid()
dev.print(file="temp/JC2004_edf.ps",horizontal=FALSE)
dev.copy2eps(file="temp/JC2004_edf.eps")

print(paste("Monte-Carlo=",min(h$mids[edf > 0.5]),"=",mean(MC.all),"= theoretical=",sum(1/(1:99)), 
             "q_0.025=",quantile(MC.all,0.025),"q_0.975=",quantile(MC.all,0.975)))

my.col <- rgb( c(seq(0.2,0.9,length=8),rep(1,5),seq(0.9,0.2,length=8)),
               c(seq(0.2,0.9,length=8),rep(1,5),seq(0.9,0.2,length=8)),
               c(seq(0.2,0.9,length=8),rep(1,5),seq(0.9,0.2,length=8)) )

if (is.null(season)) fname <- paste("temp/JC2004_all_ar1.",sep="") else
                     fname <- paste("temp/JC2004_all_ar1_",seas,".",sep="")
filled.contour(lon,lat,AR1,main="Autocorrelation",col=my.col,levels=seq(-1,1,length=21),
      sub=paste(length(gcms),"scenarios",seas))
las <- 1

    # From filled.contour in base
    mar.orig <- (par.orig <- par(c("mar","las","mfrow")))$mar
    on.exit(par(par.orig))

    w <- (3 + mar.orig[2]) * par('csi') * 2.54
    layout(matrix(c(2, 1), nc=2), widths=c(1, lcm(w)))

    par(las = las)
    mar <- mar.orig
    mar[4] <- 1
    par(mar=mar)

contour(lon,lat,AR1,add=TRUE,cex=1.6,levels=seq(-1,1,length=21))
addland(col="black",lwd=3)
dev.print(file=paste(fname,"ps",sep=""),horizontal=FALSE)
dev.copy2eps(file=paste(fname,"eps",sep=""))
dev2bitmap(file=paste(fname,"jpg",sep=""),type="jpeg"); dev.off()
}

record.analysis(bw=FALSE)

#For other seasons, first set force.read=TRUE:
#record.analysis(season=c(12,1,2))
#record.analysis(season=3:5)
#record.analysis(season=6:8,bw=FALSE)
#record.analysis(season=9:11)
