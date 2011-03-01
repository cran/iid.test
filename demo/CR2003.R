detrend <- function(st.obj,
                    trendmod="lm(y ~ x + I(x^2) + I(x^3),data=X)") {
  nt <- length(st.obj$val[,1])
  x <- seq(0,1,length=nt)
  trend <- matrix(rep(NA,12*nt),nt,12)
  for (im in 1:12) {
    y <- st.obj$val[,im]
    X <- data.frame(y=y,x=x)
    trnd <- predict(eval(parse(text=trendmod)),newdata=X)
    st.obj$val[,im] <- y - trnd
    trend[,im] <- trnd
  }
  st.obj$trend <- trend
  rm(y,x,nt,X)
  invisible(st.obj)
}


#-----------------------------------


data(CR2003)
attach(CR2003)
N.test <- 1000
ar <- 0.5
ldt=TRUE
nt <- 200

#-----------------------------------

if (ar == 0) mc.descr <- "white noise" else
             mc.descr <- paste("red noise ( ar=",ar,")")

AR <- matrix(rep(20,12),20,12)
for (i in 1:12) {
  AR[,i] <- acf(oslo$val[,i],plot=FALSE)$acf
}

acf(oslo$val[,1],main="Oslo-Blindern auto-correlation (Jan--Dec)",lwd=3)
grid()
for (i in 1:12) {
  lines(0:20+(i-1)/12,c(AR[,i],AR[20,i]),type="S",lwd=3)
}


if (ldt) {
  x11()
  plot(yy,oslo$val[,1],
       main="Remove trend")
  grid()
  oslo <- detrend(oslo)
  lines(yy,oslo$trend[,1],lwd=3,col="red")
  lines(yy,oslo$trend[,1]+oslo$val[,1],lty=3,col="red")
  sodankylae <- detrend(sodankylae)
  bergen <- detrend(bergen)
  bodoe <- detrend(bodoe)
  haparanda <- detrend(haparanda)
  jokkmokk <- detrend(jokkmokk)
  koebenhavn <- detrend(koebenhavn)
  turku <- detrend(turku)
  helsinki <- detrend(helsinki)
  tampere <- detrend(tampere)
  jyvaeskylae <- detrend(jyvaeskylae)
  lappeenranta <- detrend(lappeenranta)
  torshavn <- detrend(torshavn)
  mariehamn <- detrend(mariehamn)
  stockholm <- detrend(stockholm)
  goeteborg <- detrend(goeteborg)
  stensele <- detrend(stensele)
  karesuando <- detrend(karesuando)
  falun <- detrend(falun)
  karlstad <- detrend(karlstad)
  sveg <- detrend(sveg)
  oestersund <- detrend(oestersund)
  vestervig <- detrend(vestervig)
  varnes <- detrend(varnes)
  haernoesand <- detrend(haernoesand)
}
print(dim(oslo$val))


ny <- length(yy)
n.series <- 25*12
recs <- rep(-999,n.series)
recs.2 <- rep(-999,n.series)
y <- rep(0,ny)
y.2 <- rep(0,ny)
for (i in 1:ny) {
  ii <- ny-i+1
  vec <- c(oslo$val[i,],sodankylae$val[i,],bergen$val[i,],
           bodoe$val[i,],haparanda$val[i,],jokkmokk$val[i,],
           koebenhavn$val[i,],turku$val[i,],helsinki$val[i,],
           tampere$val[i,],jyvaeskylae$val[i,],lappeenranta$val[i,],
           torshavn$val[i,],mariehamn$val[i,],
           stockholm$val[i,],goeteborg$val[i,],stensele$val[i,],
           karesuando$val[i,],falun$val[i,],karlstad$val[i,],
           sveg$val[i,],oestersund$val[i,],vestervig$val[i,],
           varnes$val[i,],haernoesand$val[i,])
#  vec.2 <- c(oslo$val[ii,],sodankylae$val[ii,],bergen$val[ii,],
#             bodoe$val[ii,],haparanda$val[ii,],jokkmokk$val[ii,],
#             koebenhavn$val[ii,],turku$val[ii,],helsinki$val[ii,],
#             tampere$val[ii,],jyvaeskylae$val[ii,],lappeenranta$val[ii,],
#             torshavn$val[ii,],mariehamn$val[ii,],
#             stockholm$val[ii,],goeteborg$val[ii,],stensele$val[ii,],
#             karesuando$val[ii,],falun$val[ii,],karlstad$val[ii,],
#             sveg$val[ii,],oestersund$val[ii,],vestervig$val[ii,],
#             varnes$val[ii,],haernoesand$val[ii,])
  vec.2 <- vec +
    c(oslo$trend[i,],sodankylae$trend[i,],bergen$trend[i,],
      bodoe$trend[i,],haparanda$trend[i,],jokkmokk$trend[i,],
      koebenhavn$trend[i,],turku$trend[i,],helsinki$trend[i,],
      tampere$trend[i,],jyvaeskylae$trend[i,],lappeenranta$trend[i,],
      torshavn$trend[i,],mariehamn$trend[i,],
      stockholm$trend[i,],goeteborg$trend[i,],stensele$trend[i,],
      karesuando$trend[i,],falun$trend[i,],karlstad$trend[i,],
      sveg$trend[i,],oestersund$trend[i,],vestervig$trend[i,],
      varnes$trend[i,],haernoesand$trend[i,])
  
  i.rec <- (vec > recs) & is.finite(vec)
  i.rec.2 <- (vec.2 > recs.2)  & is.finite(vec.2)
  if (sum(i.rec,na.rm=TRUE)>0) {
    recs[i.rec] <- vec[i.rec]
  }
  if (sum(i.rec.2,na.rm=TRUE)>0) {
    recs.2[i.rec.2] <- vec.2[i.rec.2]
  }  
  y[i] <- sum(i.rec,na.rm=TRUE)/n.series
  y.2[i] <- sum(i.rec.2,na.rm=TRUE)/n.series
}

chi2.test <- function(x,p) {
  # Press et al, (1989) Num. recipes, p. 517
  chi2 <- sum( (x-p)^2/(x+p), na.rm=TRUE )
  print(paste("Chi2=",round(chi2)))
  p <- pchisq(chi2,df=length(x)-1)
  invisible(list(p=p,chi2=chi2))
}

# Carry out a Monte-Carlo simulation for the analysis of coincidental
# record-events:

print("Monte-Carlo simulation for coincidence statistics..")
mc.counts <- matrix(rep(0,ny*N.test),ny,N.test)

for (i.mc in 1:N.test) {
  mc.res <- matrix(rnorm(ny*n.series),ny,n.series)
  recs.mc <- rep(-999,n.series)
  for (i in 1:ny) {
    ii <- ny-i+1
    vec <- mc.res[i,]
    i.rec <- (vec > recs.mc)
    if (sum(i.rec,na.rm=TRUE)>0) {
      recs.mc[i.rec] <- vec[i.rec]
    }
    mc.counts[i,i.mc] <- sum(i.rec,na.rm=TRUE)
  }
}

x11()
rec.h1 <- hist((y[1:ny]*n.series),breaks=0:300,plot=FALSE)
rec.h2 <- hist((y.2[1:ny]*n.series),breaks=0:300,plot=FALSE)
expc.no <- hist(as.vector(mc.counts[1:ny,]),breaks=0:300,plot=FALSE)
chi2.1 <- chi2.test(x=rec.h1$counts,p=expc.no$counts/N.test)
chi2.2 <- chi2.test(x=rec.h2$counts,p=expc.no$counts/N.test)
plot(rec.h2$breaks,c(rec.h2$counts,rec.h2$counts[length(rec.h1$counts)]),
      type="s",col="black",lwd=5,xlim=c(0,25),
     main="Coincidence of record-events",
     xlab="Number of coinciding events",ylab="Counts",
     sub=paste("Empirical data vs random incidences for n=",ny,
       "and",n.series,"series"))
grid()
lines(rec.h1$breaks,c(rec.h1$counts,rec.h1$counts[length(rec.h1$counts)]),
      type="s",col="grey",lwd=3)
lines(expc.no$mids,expc.no$counts/N.test,col="grey95")
lines(expc.no$mids,expc.no$counts/N.test,lty=2)
points(expc.no$mids,expc.no$counts/N.test,pch=20,cex=0.8)

legend(12,10.5,c(paste("original (chi2=",round(chi2.1$chi2,1),
                     ", p=",round(chi2.1$p,2),")",sep=""),
               paste("detrended (chi2=",round(chi2.2$chi2,1),
                     ", p=",round(chi2.2$p,2),")",sep=""),
               paste("Monte-Carlo:",mc.descr)),cex=0.8,
       col=c("black","grey","black"),lty=c(1,1,2),lwd=c(5,3,1),
       pch=c(26,26,20),bg="grey95")
dev.copy2eps(file=paste("paper17b_0.eps",sep=""))


p.1 <-1/(1:ny)
p.N <- 1 - (1-1/(1:ny))^n.series
p.10 <- 1 - (1-1/(1:ny))^10
p.50 <- 1 - (1-1/(1:ny))^50
p.100 <- 1 - (1-1/(1:ny))^100
p.300 <- 1 - (1-1/(1:ny))^300
p.2 <- 1 - (1-1/(1:ny))^2
p <- p.1

print("Monte Carlo")
mc.sums <- matrix(NA,ny,N.test)
q.025 <- rep(NA,ny)
q.975 <- rep(NA,ny)
rec <- matrix(rep(-999,n.series*N.test),n.series,N.test)
vec <- matrix(rnorm(ny*N.test*n.series),ny,n.series*N.test)
dim(vec) <- c(ny,n.series,N.test)
if (ar > 0) {
  print("with auto-correlation")
  vec[,2:n.series,] <- (1-ar)*vec[,2:n.series,] + ar*vec[,1:(n.series-1),]
}

for (i in 1:ny) {
  for (ii in 1:N.test) {
    i.rec <- (as.vector(vec[i,,ii]) > rec[,ii])
    if (sum(i.rec,na.rm=TRUE)>0) {
      rec[i.rec,ii] <- vec[i,i.rec,ii]
    }
    mc.sums[i,ii] <- sum(i.rec,na.rm=TRUE)/n.series
  }
  q.025[i] <- quantile(mc.sums[i,],0.025)
  q.975[i] <- quantile(mc.sums[i,],0.975)
}

Pr1 <- p[is.finite(y)]/sum(p[is.finite(y)])
Pr2 <- p[is.finite(y.2)]/sum(p[is.finite(y.2)])
chi.2 <- chisq.test(y[is.finite(y)],p=Pr1)
chi.2.2 <- chisq.test(y.2[is.finite(y.2)],p=Pr2)
chi.2.mc <- rep(NA,N.test)
e.sums <- rep(NA,N.test)
for (ii in 1:N.test) {
  test <- chisq.test(mc.sums[,ii],p=p/sum(p,na.rm=TRUE))$statistic
  chi.2.mc[ii] <- round(as.numeric(test),2)
  e.sums[ii] <- sum(mc.sums[,ii])
}
chi.squared.results <- paste("Chi-squared=",
                             round(as.numeric(chi.2$statistic),2),"/",
                             round(as.numeric(chi.2.2$statistic),2),
                             " (Monte-Carlo 95% conf.lev.=",
                             quantile(chi.2.mc,0.95),")",sep="")
y.log <- log(y[is.finite(y)])
y.log.2 <- log(y.2[is.finite(y.2)])
p.log <- log(p[is.finite(y)])
ii <- is.finite(y.log) & is.finite(p.log)
ii.2 <- is.finite(y.log.2) & is.finite(p.log)
pairs <- data.frame(y=y.log[ii],x=p.log[ii])
pairs.2 <- data.frame(y=y.log.2[ii.2],x=p.log[ii.2])
fit <- lm(y ~ x,data=pairs)
fit.2 <- lm(y ~ x,data=pairs.2)

x11()
plot(p,type="l",col="grey30",lty=3,
     main="New records",
     sub="25  x 12 monthly maximum temperature series",
     ylab="Record density",xlab="Time")
polygon(c(spline(q.025)$x,reverse(spline(q.975)$x)),
        c(spline(q.025)$y,reverse(spline(q.975)$y)),
        col="wheat")
lines(y.2,type="s",col="grey40",lwd=6)
lines(y,type="s",col="black",lwd=2)
lines(p,col="grey30",lty=2)
grid()
legend(length(p)/2,1,
       c("E/N detrended","E/N original","p","Null hyp. 95 conf."),
       col=c("black","grey40","grey30","grey85"),lwd=c(2,6,1,10),
       lty=c(1,1,2,1),bg="grey95")


x11()
plot(log(p),log(y),pch=20,col="black",
     main="Expected and observed records",
     sub="25  x 12 monthly maximum temperature",
     ylab="log(E/N)",
     xlab="log(theoretical record density)")
lines(log(p),log(q.975),lty=2,col="grey30")
lines(log(p),log(q.025),lty=2,col="grey30")
lines(c(-10,10),c(-10,10),lty=2,col="grey30")
abline(fit,col="black",lwd=1)
abline(fit.2,col="grey",lwd=2)
points(log(p),log(y.2),pch=5,col="black")
points(log(p),log(y),pch=20,col="grey")
points(log(p),log(y),pch=21)
grid()
text(-2.5,0,chi.squared.results)
#text(-2.5,-0.2,paste("P-value=",round(as.numeric(chi.2$p.value),2)))


lons <- c(oslo$lon,sodankylae$lon,bergen$lon,bodoe$lon,haparanda$lon,
          jokkmokk$lon,koebenhavn$lon,turku$lon,helsinki$lon,
          tampere$lon,jyvaeskylae$lon,lappeenranta$lon,
          torshavn$lon,mariehamn$lon,stockholm$lon,goeteborg$lon,
          stensele$lon,karesuando$lon,falun$lon,karlstad$lon,
           sveg$lon,oestersund$lon,vestervig$lon,varnes$lon,haernoesand$lon)
lats <- c(oslo$lat,sodankylae$lat,bergen$lat,bodoe$lat,haparanda$lat,
          jokkmokk$lat,koebenhavn$lat,turku$lat,helsinki$lat,
          tampere$lat,jyvaeskylae$lat,lappeenranta$lat,
          torshavn$lat,mariehamn$lat,stockholm$lat,goeteborg$lat,
          stensele$lat,karesuando$lat,falun$lat,karlstad$lat,
          sveg$lat,oestersund$lat,vestervig$lat,varnes$lat,haernoesand$lat)
          

X <- cbind(oslo$val,sodankylae$val,bergen$val,bodoe$val,haparanda$val,
               jokkmokk$val,koebenhavn$val,turku$val,
               helsinki$val,tampere$val,jyvaeskylae$val,lappeenranta$val,
               torshavn$val,mariehamn$val,stockholm$val,
               goeteborg$val,stensele$val,karesuando$val,falun$val,
               karlstad$val,sveg$val,oestersund$val,vestervig$val,
               varnes$val,haernoesand$val)
X[is.na(X)] <- 0
X.cor <- cor(X)
dims <- dim(X.cor)
X.cor[X.cor==1] <- NA
print(summary(as.vector(X.cor)))
print(paste("sum(y)=",sum(y),"log(length(y)) + 1=",round(log(length(y)) + 1,1),
            "M-C 95%:",quantile(e.sums,0.025),"-",quantile(e.sums,0.975)))



E <- function(n) {
  E <- log(n) + 1
  E
}

ns <- seq(1,nt,by=5)
plot(ns,ns,type="n",ylim=c(0,15),
     main="Expected number of record-events",
     xlab="Record length (n)",
     ylab="Number of record-events / Expectated number",
     sub=paste("Theory vs Monte-Carlo (N=",N.test,") simulation",sep=""),
     cex=1.2)
#lines(ns,E(ns))
#points(cumsum(1/(1:200)),pch=20,cex=0.4)
grid()
r.map <- matrix(rep(0,length(ns)*length(ns)),length(ns),length(ns))

iii <- 0
E.mc <- rep(NA,length(ns))
for (i in ns) {
    iii <- iii+1
    for (it in 1:N.test) {
      vec <- rnorm(i)
      rec.val <- vec[1]
      no.rec <- 1
      for (ii in 1:length(vec)) {
        if (rec.val > vec[ii]) {
          rec.val <- vec[ii]
          no.rec <- no.rec + 1
        }
      }
      r.map[iii,no.rec] <- r.map[iii,no.rec]+1
      points(i,no.rec,pch=20,cex=0.3,col="grey90")
    }
    N.iii <- sum(r.map[iii,],na.rm=TRUE)
    i.x <- 1:iii
    E.mc[iii] <- sum(i.x*r.map[iii,i.x]/N.iii,na.rm=TRUE)
    q025 <- min(i.x[cumsum(r.map[iii,i.x]) >= 0.025*sum(r.map[iii,i.x])])
    q975 <- min(i.x[cumsum(r.map[iii,i.x]) >= 0.975*sum(r.map[iii,i.x])])
#    points(i,E.mc[iii],pch=20,cex=0.3,col="black")
    lines(rep(i,2),c(q025,q975),lty=1,col="grey70")
}
contour(ns,1:length(ns),r.map,col="grey40",lwd=1,add=TRUE,
        levels=seq(0,N.test,by=10))
#lines(ns,E.mc,lwd=2,lty=2)
lines(cumsum(1/(1:200)),lwd=3,lty=1)
#lines(ns,E(ns),lwd=3)

legend(10,14.5,c("cumsum(1/n)",
                 paste("Monte-Carlo:",mc.descr)),
       cex=0.8,lty=c(1,1),lwd=c(3,1),
       col=c("black","grey40"),bg="grey95")







ar <- 0.0
N.test <- 100
if (ar == 0) mc.descr <- "white noise" else
             mc.descr <- paste("red noise ( ar=",ar,")")

yy0 <- 0
yy <- 1:92
ny <- length(yy)
n.series <- 25*12
recs <- rep(-999,n.series)
recs.2 <- rep(-999,n.series)
y <- rep(0,ny)
y.2 <- rep(0,ny)
for (i in 1:ny) {
  ii <- ny-i+1
  vec <- rnorm(n.series)
#  vec.2 <- rep(vec[1:12],25)
  vec.2 <- rep(vec[1:60],5)
  i.rec <- (vec > recs)
  i.rec.2 <- (vec.2 > recs.2)
  if (sum(i.rec,na.rm=TRUE)>0) {
    recs[i.rec] <- vec[i.rec]
  }
  if (sum(i.rec.2,na.rm=TRUE)>0) {
    recs.2[i.rec.2] <- vec.2[i.rec.2]
  }  
  y[i] <- sum(i.rec,na.rm=TRUE)/n.series
  y.2[i] <- sum(i.rec.2,na.rm=TRUE)/n.series
}


p <-1/(1:ny)

print("Monte Carlo")
mc.sums <- matrix(NA,ny,N.test)
q.025 <- rep(NA,ny)
q.975 <- rep(NA,ny)
rec <- matrix(rep(-999,n.series*N.test),n.series,N.test)
vec <- matrix(rnorm(ny*N.test*n.series),ny,n.series*N.test)
dim(vec) <- c(ny,n.series,N.test)
if (ar > 0) {
  print("with auto-correlation")
  vec[2:ny,,] <- (1-ar)*vec[2:ny,,] + ar*vec[1:(ny-1),,]
}

for (i in 1:ny) {
  for (ii in 1:N.test) {
    i.rec <- (as.vector(vec[i,,ii]) > rec[,ii])
    if (sum(i.rec,na.rm=TRUE)>0) {
      rec[i.rec,ii] <- vec[i,i.rec,ii]
    }
    mc.sums[i,ii] <- sum(i.rec,na.rm=TRUE)/n.series
  }
  q.025[i] <- quantile(mc.sums[i,],0.025)
  q.975[i] <- quantile(mc.sums[i,],0.975)
}

Pr1 <- p[is.finite(y)]/sum(p[is.finite(y)])
Pr2 <- p[is.finite(y.2)]/sum(p[is.finite(y.2)])
chi.2 <- chisq.test(y[is.finite(y)],p=Pr1)
chi.2.2 <- chisq.test(y.2[is.finite(y.2)],p=Pr2)
chi.2.mc <- rep(NA,N.test)
e.sums <- rep(NA,N.test)
for (ii in 1:N.test) {
  test <- chisq.test(mc.sums[,ii],p=p/sum(p))$statistic
  chi.2.mc[ii] <- round(as.numeric(test),2)
  e.sums[ii] <- sum(mc.sums[,ii])
}
chi.squared.results <- paste("Chi-squared=",
                             round(as.numeric(chi.2$statistic),2),"/",
                             round(as.numeric(chi.2.2$statistic),2),
                             " (Monte-Carlo 95% conf.lev.=",
                             quantile(chi.2.mc,0.95),")",sep="")
y.log <- log(y[is.finite(y)])
y.log.2 <- log(y.2[is.finite(y.2)])
p.log <- log(p[is.finite(y)])
ii <- is.finite(y.log) & is.finite(p.log)
ii.2 <- is.finite(y.log.2) & is.finite(p.log)
pairs <- data.frame(y=y.log[ii],x=p.log[ii])
pairs.2 <- data.frame(y=y.log.2[ii.2],x=p.log[ii.2])
fit <- lm(y ~ x,data=pairs)
fit.2 <- lm(y ~ x,data=pairs.2)

x11()
plot(p,type="l",col="black",lty=3,
     main="New records",
     sub=paste("25 x 12 monthly precipitation series.",
       "Monte-Carlo w. inter-station correlation & AR=",ar),
     ylab="Record density",xlab="Time")
polygon(c(spline(q.025)$x,reverse(spline(q.975)$x)),
        c(spline(q.025)$y,reverse(spline(q.975)$y)),
        col="grey85")
lines(y.2,type="s",col="grey40",lwd=4)
lines(y,type="s",col="black",lwd=2)
lines(p,col="grey30",lty=2)
grid()
legend(length(p)/2,1,
       c("E/N indendent","E/N 25 identical stations","p","Null hyp. 95 conf."),
       col=c("black","grey40","grey30","grey85"),lwd=c(2,4,1,10),
       lty=c(1,1,2,1),bg="grey95")
dev.copy2eps(file=paste("paper17d_1.eps",sep=""))

x11()
plot(log(p),log(y),pch=20,col="black",
     main="Expected and observed records",
     sub=paste("25 x 12 monthly precipitation.",
       "Monte-Carlo w. inter-station correlation & AR=",ar),
     ylab="log(E/N)",
     xlab="log(theoretical record density)")
lines(log(p),log(q.975),lty=2,col="grey30")
lines(log(p),log(q.025),lty=2,col="grey30")
lines(c(-10,10),c(-10,10),lty=2,col="grey30")
abline(fit.2,col="black",lwd=1)
abline(fit,col="grey",lwd=2)
points(log(p),log(y.2),pch=5,col="black")
points(log(p),log(y),pch=20,col="grey")
points(log(p),log(y),pch=21)
grid()
text(-2.5,0,chi.squared.results)
#text(-2.5,-0.2,paste("P-value=",round(as.numeric(chi.2$p.value),2)))
dev.copy2eps(file=paste("paper17d_2.eps",sep=""))

print(paste("sum(y)=",sum(y),"log(length(y)) + 1=",round(log(length(y)) + 1,1),
            "M-C 95%:",quantile(e.sums,0.025),"-",quantile(e.sums,0.975)))
e.sums2 <- colSums(mc.sums)
print(paste("sum(y.2)=",sum(y.2),"log(length(y.2)) + 1=",round(log(length(y.2)) + 1,1),
            "M-C 95%:",quantile(e.sums2,0.025),"-",quantile(e.sums2,0.975)))


