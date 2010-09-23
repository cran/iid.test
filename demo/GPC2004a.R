rm(list=ls())
do.monte.carlo <- FALSE


# Load nesescary packages:

library(iid.test)
library(clim.pact)
library(evd)
library(xtable)

# Set font
str(Hershey)

# Number of tests
N.test <- 1000
wd0 <- getwd()
data(giss)
ns <- length(giss$Location)
x11()
plot(c(-180,180),c(-90,90),type="n",xlab="degrees east",ylab="degrees north")
addland()
grid()

points(giss$lon,giss$lat,pch=20,col="grey90")
points(giss$lon,giss$lat,pch=21,col="grey30")
#isel <- giss$population < 500000
#points(giss$lon[isel],giss$lat[isel],pch=20,col="blue",cex=0.6)
#isel <- giss$first < 1875
isel <- c(1,6,8,10,16,18,19,23,24,25,28,29,30,32,35,36,38)
points(giss$lon[isel],giss$lat[isel],pch=20,cex=1.75)
points(giss$lon[isel]+0.2,giss$lat[isel]+0.2,pch=20,cex=1.90,col="grey30")
points(giss$lon[isel]+0.3,giss$lat[isel]+0.5,pch=20,cex=1.00,col="grey50")
points(giss$lon[isel]+0.4,giss$lat[isel]+0.7,pch=20,cex=0.50,col="grey70")
points(giss$lon[isel]+0.5,giss$lat[isel]+0.9,pch=20,cex=0.30,col="white")
#points(giss$lon[isel],giss$lat[isel],pch=20,col="darkred",cex=1)
for (i in 1:length(giss$lon[isel])) text(giss$lon[isel[i]],
                                         giss$lat[isel[i]]-2.5,
                                         as.character(giss$Location[isel[i]]),
                                         cex=0.5)
ns <- length(isel)*12
nt <- max(giss$last[isel]) - max(giss$first[isel]) + 1
t2m.all <- matrix(rep(NA,nt*ns),nt,ns)
yy.all <- seq(max(giss$first[isel]),max(giss$last[isel]),by=1)
  
i <- 1
for (fname in as.character(giss$file.name[isel])) {
  eval(parse(text=paste("data(",fname,")",sep="")))
  t2m <- eval(parse(text=fname))
  i1 <- is.element(yy.all,t2m$YEAR)
  i2 <- is.element(t2m$YEAR,yy.all)
  t2m.all[i1,i] <- t2m$JAN[i2]
  t2m.all[i1,i+1] <- t2m$FEB[i2]
  t2m.all[i1,i+2] <- t2m$MAR[i2]
  t2m.all[i1,i+3] <- t2m$APR[i2]
  t2m.all[i1,i+4] <- t2m$MAY[i2]
  t2m.all[i1,i+5] <- t2m$JUN[i2]
  t2m.all[i1,i+6] <- t2m$JUL[i2]
  t2m.all[i1,i+7] <- t2m$AUG[i2]
  t2m.all[i1,i+8] <- t2m$SEP[i2]
  t2m.all[i1,i+9] <- t2m$OCT[i2]
  t2m.all[i1,i+10] <- t2m$NOV[i2]
  t2m.all[i1,i+11] <- t2m$DEC[i2]  
  print(paste(i,fname," No. missing=",sum(t2m.all[i1,i:(i+11)]>= 999)))
  i <- i + 12
}
t2m.all[t2m.all >= 990] <- NA

print("records:")
x11()
plot(c(-10,nt),c(1,ns*2+25),type="n",xlab="Time (years)",ylab="Series",
     main="Record incidence",cex=1.5)
for (i in 1:ns) {
  lines(c(1,nt),rep(i,2),lwd=1,col="grey40")
  lines(c(1,nt),rep(i,2)+0.3,lwd=1,col="grey90")
  lines(c(1,nt),rep(i,2)+ns+20,lwd=1,col="grey40")
  lines(c(1,nt),rep(i,2)+ns+20.3,lwd=1,col="grey90")
}
lines(c(1,nt),rep(ns,2),lwd=3,col="black")
lines(c(1,nt),rep(ns+20,2),lwd=3,col="black")
par.0 <- par()
par(srt=90)
text(-8,100,"Forward",cex=2,vfont=c("sans serif","italic"))
text(-8,300,"Backward",cex=2,vfont=c("sans serif","italic"))
par(par.0)
recs.1 <- rep(-999,ns)
recs.2 <- rep(-999,ns)
y.1 <- rep(0,nt)
y.2 <- rep(0,nt)

for (it in 1:nt) {
  ii <- nt-it+1
  vec.1 <- t2m.all[it,]
  vec.2 <- t2m.all[ii,]
  i.rec.1 <- (vec.1 > recs.1) & is.finite(vec.1)
  i.rec.2 <- (vec.2 > recs.2) & is.finite(vec.2)
  if (sum(i.rec.1,na.rm=T)>0) {
    recs.1[i.rec.1] <- vec.1[i.rec.1]
  }
  if (sum(i.rec.2,na.rm=T)>0) {
    recs.2[i.rec.2] <- vec.2[i.rec.2]
  }  
  y.1[it] <- sum(i.rec.1,na.rm=T)/sum(is.finite(t2m.all[it,]))
  y.2[it] <- sum(i.rec.2,na.rm=T)/sum(is.finite(t2m.all[ii,]))

  points(rep(it,sum(i.rec.1,na.rm=T))+0.1,seq(1,ns,by=1)[i.rec.1]+0.1,
         pch=20,cex=1.50,col="grey30")
  points(rep(it,sum(i.rec.1,na.rm=T))+0.2,seq(1,ns,by=1)[i.rec.1]+0.2,
         pch=20,cex=0.70,col="grey50")
  points(rep(it,sum(i.rec.1,na.rm=T))+0.3,seq(1,ns,by=1)[i.rec.1]+0.3,
         pch=20,cex=0.50,col="grey70")
  points(rep(it,sum(i.rec.1,na.rm=T))+0.4,seq(1,ns,by=1)[i.rec.1]+0.4,
         pch=20,cex=0.30,col="white")

  points(rep(it,sum(i.rec.2,na.rm=T))+0.1,seq(1,ns,by=1)[i.rec.2]+ns+20.1,
         pch=20,cex=1.50,col="grey30")
  points(rep(it,sum(i.rec.2,na.rm=T))+0.2,seq(1,ns,by=1)[i.rec.2]+ns+20.2,
         pch=20,cex=0.70,col="grey50")
  points(rep(it,sum(i.rec.2,na.rm=T))+0.3,seq(1,ns,by=1)[i.rec.2]+ns+20.3,
         pch=20,cex=0.50,col="grey70")
  points(rep(it,sum(i.rec.2,na.rm=T))+0.4,seq(1,ns,by=1)[i.rec.2]+ns+20.4,
         pch=20,cex=0.30,col="white")
}

p.1 <-1/(1:nt)
p.N <- 1 - (1-1/(1:nt))^ns
p.10 <- 1 - (1-1/(1:nt))^10
p.50 <- 1 - (1-1/(1:nt))^50
p.100 <- 1 - (1-1/(1:nt))^100
p.300 <- 1 - (1-1/(1:nt))^300
p.2 <- 1 - (1-1/(1:nt))^2
p <- p.1

print(paste("Monte-Carlo:",nt,ns,N.test))
mc.sums <- matrix(rnorm(nt*N.test),nt,N.test)
q.025 <- rep(NA,nt)
q.975 <- rep(NA,nt)
rec <- matrix(rep(-999,ns*N.test),ns,N.test)

for (i in 1:nt) {
  for (ii in 1:N.test) {
    vec <- rnorm(ns)
    i.rec <- (vec > rec[,ii])
    if (sum(i.rec,na.rm=T)>0) {
      rec[i.rec,ii] <- vec[i.rec]
    }
    mc.sums[i,ii] <- sum(i.rec,na.rm=TRUE)/ns
  }
  q.025[i] <- quantile(mc.sums[i,],0.025)
  q.975[i] <- quantile(mc.sums[i,],0.975)
}

sumy <- sum(y.1[is.finite(y.1)]); sump <- sum(p[is.finite(y.1)])
chi.2.1 <- chisq.test(y.1[is.finite(y.1)]/sumy,p=p[is.finite(y.1)]/sump)
sumy <- sum(y.2[is.finite(y.2)]); sump <- sum(p[is.finite(y.2)])
chi.2.2 <- chisq.test(y.2[is.finite(y.2)]/sumy,p=p[is.finite(y.2)]/sump)
chi.2.mc <- rep(NA,N.test)
e.sums <- rep(NA,N.test)
for (ii in 1:N.test) {
  summc <- sum(mc.sums[,ii]); 
  test <- chisq.test(mc.sums[,ii]/summc,p=p/sum(p))$statistic
  chi.2.mc[ii] <- round(as.numeric(test),2)
  e.sums[ii] <- sum(mc.sums[,ii])
}
chi.squared.results <- paste("Chi-squared=",
                             round(as.numeric(chi.2.1$statistic),2),"/",
                             round(as.numeric(chi.2.2$statistic),2),
                             " (Monte-Carlo 95% conf.lev.=",
                             quantile(chi.2.mc,0.95),")",sep="")
y.log.1 <- log(y.1[is.finite(y.1)])
y.log.2 <- log(y.2[is.finite(y.2)])
p.log <- log(p[is.finite(y.1)])
ii.1 <- is.finite(y.log.1) & is.finite(p.log)
ii.2 <- is.finite(y.log.2) & is.finite(p.log)
pairs.1 <- data.frame(y=y.log.1[ii.1],x=p.log[ii.1])
pairs.2 <- data.frame(y=y.log.2[ii.2],x=p.log[ii.2])
fit.1 <- lm(y ~ x,data=pairs.1)
fit.2 <- lm(y ~ x,data=pairs.2)

x11()
plot(p,type="l",col="grey30",lty=3,
     main="New records",
     sub="25  x 12 monthly mean temperature series",
     ylab="Record density",xlab="Time",cex=1.5)
polygon(c(spline(q.025)$x,reverse(spline(q.975)$x)),
        c(spline(q.025)$y,reverse(spline(q.975)$y)),
        col="grey85")
lines(y.2,type="s",col="grey40",lwd=4)
lines(y.1,type="s",col="black",lwd=2)
lines(p,col="grey30",lty=2)
grid()
legend(length(p)/2,1,
       c("E/N forward","E/N backward","p","Null hyp. 95 conf."),
       col=c("black","grey40","grey30","grey85"),lwd=c(2,4,1,10),
       lty=c(1,1,2,1),bg="grey95")

x11()
plot(log(p),log(y.1),pch=20,col="black",
     main="Expected and observed records",
     sub="25  x 12 monthly mean temperature series",
     ylab="log(E/N)",
     xlab="log(theoretical record density)",cex=1.5)
lines(log(p),log(q.975),lty=2,col="grey30")
lines(log(p),log(q.025),lty=2,col="grey30")
lines(c(-10,10),c(-10,10),lty=2,col="grey30")
abline(fit.2,col="black",lwd=1)
abline(fit.1,col="grey",lwd=2)
points(log(p),log(y.2),pch=5,col="black")
points(log(p),log(y.1),pch=20,col="grey")
points(log(p),log(y.1),pch=21)
grid()
text(-2.5,0,chi.squared.results)

x11()
plot(p.1,type="l",col="black",lwd=3,
     main="Probability of new records",
     sub="Different number of independent series",
     ylab="Record density",xlab="Time",cex=1.5)
lines(p.N,col="grey40",lwd=3,lty=2)
grid()
legend(85,1,c("E_n/N   ",paste("N=",ns,"   ")),
       col=c("black","grey40"),
       lwd=c(3,3),lty=c(1,2),bg="grey95",cex=0.75)


X <- t2m.all
X[is.na(X)] <- 0
X.cor <- cor(X)
dims <- dim(X.cor)
X.cor[X.cor==1] <- NA
print(summary(as.vector(X.cor)))
print(paste("sum(y)=",sum(y.1),"log(length(y.1)) + 1=",
            round(log(length(y.1)) + 1,1),
            "M-C 95%:",quantile(e.sums,0.025),"-",quantile(e.sums,0.975)))




mc.descr <- "white noise"

  nl <- seq(1,nt,by=5)
  x11()
  par(col.axis="white",cex.axis=1.5)
  plot(nl,nl,type="n",ylim=c(-7,300),clim=c(0,130),axes=FALSE,
       main="Expected number of record-events",
       xlab="Record length (n)",
       ylab="Number of record-events / Expectated number",
       sub=paste("mean observed vs Monte-Carlo (N=",N.test,
         ") simulation",sep=""),
       cex=1.2)
  par(col.axis="black")
  axis(1,at=seq(0,nt,by=10), lty = 1, lwd = 0.5)
  yticks <- exp(seq(1,15,by=0.5))
  axis(2,at=yticks, labels=round(log(yticks),1), lty = 1, lwd = 0.5)
  lines(exp(cumsum(1/(1:200))),lwd=2,lty=2)

  iii <- 0
  i.X <- 1:length(nl)
  E.mc <- rep(NA,length(nl))
  q025 <- rep(NA,length(nl))
  q975 <- rep(NA,length(nl))
  for (i in nl) {
    iii <- iii+1
    points(i+0.0,exp(sum(y.2[1:i])),
         pch=20,cex=1.50,col="grey70")
    points(i+0.1,exp(sum(y.2[1:i]))+0.5,
         pch=20,cex=0.70,col="grey50")
    points(i+0.2,exp(sum(y.2[1:i]))+1,
         pch=20,cex=0.50,col="grey30")
    points(i+0.3,exp(sum(y.2[1:i]))+1.5,
         pch=20,cex=0.30,col="black")
      
    points(i+0.0,exp(sum(y.1[1:i])),
         pch=20,cex=1.50,col="grey30")
    points(i+0.1,exp(sum(y.1[1:i]))+0.5,
         pch=20,cex=0.70,col="grey50")
    points(i+0.2,exp(sum(y.1[1:i]))+1,
         pch=20,cex=0.50,col="grey70")
    points(i+0.3,exp(sum(y.1[1:i]))+1.5,
         pch=20,cex=0.30,col="white")

    if (i > 1) e.sums <- colSums(mc.sums[1:i,]) else
               e.sums <- 1
    E.mc[iii] <- mean(e.sums)
    q025[iii] <- quantile(e.sums,0.025)
    q975[iii] <- quantile(e.sums,0.975)
    lines(rep(i,2),exp(c(q025[iii],q975[iii])),lty=1,col="grey70")
  }
  lines(nl,exp(q025),lty=2,lwd=0.5,col="grey40")
  lines(nl,exp(q975),lty=2,lwd=0.5,col="grey40")

  upper.lim <- data.frame(y=exp(q975),x=nl)
  lower.lim <- data.frame(y=exp(q025),x=nl)
  best <- data.frame(y=exp(E.mc),x=nl)

  upper.fit <- lm(y ~ x,data=upper.lim)
  lower.fit <- lm(y ~ x,data=lower.lim)
  best.fit <- lm(y ~ x,data=best)

  abline(upper.fit,lty=3,lwd=0.5,col="grey60")
  abline(lower.fit,lty=3,lwd=0.5,col="grey60")
  abline(best.fit,lty=3,lwd=0.5,col="grey40")
  lines(nl,exp(E.mc),lwd=0.5,lty=1,col="grey40")
  points(nl,exp(E.mc),pch=19,col="grey40",cex=0.5)
  lines(exp(cumsum(1/(1:200))),lwd=2,lty=2)

  legend(5,297,c("E=cumsum(1/(1:n))",
                 paste("mean E Monte-Carlo:",mc.descr),
                 paste("95% C.I. Monte-Carlo:",mc.descr)),
       cex=0.75,lty=c(2,1,2),lwd=c(2,0.5,0.5),pch=c(26,19,26),
       col=c("black","grey40","grey40"),bg="grey95")

    points(20+0.0,-5+0.00,
         pch=20,cex=1.50,col="grey70")
    points(20+0.1,-5+0.2,
         pch=20,cex=0.70,col="grey50")
    points(20+0.2,-5+0.4,
         pch=20,cex=0.50,col="grey30")
    points(20+0.3,-5+0.6,
         pch=20,cex=0.30,col="black")
    text(30,-5,"Backward")
    
    points(60+0.0,-5+0.00,
         pch=20,cex=1.50,col="grey30")
    points(60+0.1,-5+0.2,
         pch=20,cex=0.70,col="grey50")
    points(60+0.2,-5+0.4,
         pch=20,cex=0.50,col="grey70")
    points(60+0.3,-5+0.6,
         pch=20,cex=0.30,col="white")
    text(50,-5,"Forward")
  

# Test stationarity of CET

data(cet)

#t2m.1 <- cet$JAN[cet$Year < median(cet$Year)]
#t2m.2 <- cet$JAN[cet$Year >= median(cet$Year)]
t2m.1 <- cet$AUG[cet$Year < median(cet$Year)]
t2m.2 <- cet$AUG[cet$Year >= median(cet$Year)]
hi.tail.1 <- t2m.1[t2m.1 > quantile(t2m.1,0.75)]
hi.tail.2 <- t2m.2[t2m.2 > quantile(t2m.2,0.75)]
t2m.gev.1 <- fgev(hi.tail.1)
t2m.gev.2 <- fgev(hi.tail.2)

x11()
rl.1 <- rl(t2m.gev.1)
rl.2 <- rl(t2m.gev.2,ylim=c(16,23))
lines(rl.2,lwd=3)
lines(rl.1,col="grey60",lwd=3)
points(-1/log(ppoints(t2m.gev.1$data)),sort(t2m.gev.1$data),col="grey60")
title(sub=paste("Central England August Temperature"))
grid(ny=16)
legend(0.5,23,c(paste("GEV:",min(cet$Year),"-",median(cet$Year)-1),
               paste("GEV:",median(cet$Year),"-",max(cet$Year))),
       col=c("grey60","black"),pch=c(" ","x"),lty=c(1,1),lwd=c(3,3),bg="grey95")


# Tables:

upper.tab <- xtable(upper.fit,
           caption=paste("Upper level of 95-percent conf. int.",
                         "for number of record-events in one series"))
#digits(upper.tab) <- c(2,2,2,2,2)
print.table(upper.tab,file="paper18_upper.tex")

lower.tab <- xtable(lower.fit,
           caption=paste("Upper level of 95-percent conf. int.",
                         "for number of record-events in one series"))
#digits(lower.tab) <- c(2,2,2,2,2)
print.table(lower.tab,file="paper18_lower.tex")


# coefficients for linear relationship between c.i. and n and N.
# Test1 for sensitivity to distribution

if (do.monte.carlo) {
  print("Test1: relationship between coefficients and N.test")
  print("also test for sensitivity to distribution type")
  print("This may take a while, please be patient...")

  NS <- c(204,10,20,30,50,68,75,100,500,600,1000)
  b.u.fit <- matrix(rep(NA,2*length(NS)),length(NS),2); b.l.fit <- b.u.fit;
  b.fit <- b.u.fit
  b.fit.1 <- b.fit; b.fit.2 <- b.fit; b.fit.3 <- b.fit; b.fit.4 <- b.fit
  b.fit.5 <- b.fit; b.fit.6 <- b.fit; b.fit.6 <- b.fit; b.fit.7 <- b.fit

  b.u.fit[1,] <- as.numeric(upper.fit$coefficients)
  b.l.fit[1,] <- as.numeric(lower.fit$coefficients)
  b.fit[1,]   <- as.numeric(best.fit$coefficients)
  nl <- seq(1,nt,by=5)

  mc.sums <- matrix(rnorm(nt*N.test),nt,N.test)
  mc.sums.1 <- matrix(rep(NA,nt*N.test),nt,N.test)
  mc.sums.2 <- matrix(rep(NA,nt*N.test),nt,N.test)
  mc.sums.3 <- matrix(rep(NA,nt*N.test),nt,N.test)
  mc.sums.4 <- matrix(rep(NA,nt*N.test),nt,N.test)
  mc.sums.5 <- matrix(rep(NA,nt*N.test),nt,N.test)
  mc.sums.6 <- matrix(rep(NA,nt*N.test),nt,N.test)
  mc.sums.7 <- matrix(rep(NA,nt*N.test),nt,N.test)
  for (iv in 2:length(NS)) {
    rec <- matrix(rep(-999,NS[iv]*N.test),NS[iv],N.test)
    rec.1 <- matrix(rep(-999,NS[iv]*N.test),NS[iv],N.test)
    rec.2 <- matrix(rep(-999,NS[iv]*N.test),NS[iv],N.test)
    rec.3 <- matrix(rep(-999,NS[iv]*N.test),NS[iv],N.test)
    rec.4 <- matrix(rep(-999,NS[iv]*N.test),NS[iv],N.test)
    rec.5 <- matrix(rep(-999,NS[iv]*N.test),NS[iv],N.test)
    rec.6 <- matrix(rep(-999,NS[iv]*N.test),NS[iv],N.test)
    rec.7 <- matrix(rep(-999,NS[iv]*N.test),NS[iv],N.test)

    for (i in 1:nt) {
      for (ii in 1:N.test) {
        vec <- rnorm(NS[iv])
        vec.1 <- rgamma(NS[iv], shape=1)
        vec.2 <- rgamma(NS[iv], shape=10)
        vec.3 <- rgamma(NS[iv], shape=0.1)
        vec.4 <- rgev(NS[iv])
        vec.5 <- rgev(NS[iv], shape=10)
        vec.6 <- rgev(NS[iv],shape=-1)
        vec.7 <- rbinom(NS[iv],size=1000,p=0.1)

        i.rec <- (vec > rec[,ii])
        i.rec.1 <- (vec.1 > rec.1[,ii])
        i.rec.2 <- (vec.2 > rec.2[,ii])
        i.rec.3 <- (vec.3 > rec.3[,ii])
        i.rec.4 <- (vec.4 > rec.4[,ii])
        i.rec.5 <- (vec.5 > rec.5[,ii])
        i.rec.6 <- (vec.6 > rec.6[,ii])
        i.rec.7 <- (vec.7 > rec.7[,ii])
        if (sum(i.rec,na.rm=T)>0) {
          rec[i.rec,ii] <- vec[i.rec]
        }
        if (sum(i.rec.1,na.rm=T)>0) {
          rec.1[i.rec.1,ii] <- vec.1[i.rec.1]
        }
        if (sum(i.rec.2,na.rm=T)>0) {
          rec.2[i.rec.2,ii] <- vec.2[i.rec.2]
        }
        if (sum(i.rec.3,na.rm=T)>0) {
          rec.3[i.rec.3,ii] <- vec.3[i.rec.3]
        }
        if (sum(i.rec.4,na.rm=T)>0) {
          rec.4[i.rec.4,ii] <- vec.4[i.rec.4]
        }
        if (sum(i.rec.5,na.rm=T)>0) {
          rec.5[i.rec.5,ii] <- vec.5[i.rec.5]
        }
        if (sum(i.rec.6,na.rm=T)>0) {
          rec.6[i.rec.6,ii] <- vec.6[i.rec.6]
        }
        if (sum(i.rec.7,na.rm=T)>0) {
          rec.7[i.rec.7,ii] <- vec.7[i.rec.7]
        }
        mc.sums[i,ii] <- sum(i.rec,na.rm=TRUE)/NS[iv]
        mc.sums.1[i,ii] <- sum(i.rec.1,na.rm=TRUE)/NS[iv]
        mc.sums.2[i,ii] <- sum(i.rec.2,na.rm=TRUE)/NS[iv]
        mc.sums.3[i,ii] <- sum(i.rec.3,na.rm=TRUE)/NS[iv]
        mc.sums.4[i,ii] <- sum(i.rec.4,na.rm=TRUE)/NS[iv]
        mc.sums.5[i,ii] <- sum(i.rec.5,na.rm=TRUE)/NS[iv]
        mc.sums.6[i,ii] <- sum(i.rec.6,na.rm=TRUE)/NS[iv]
        mc.sums.7[i,ii] <- sum(i.rec.7,na.rm=TRUE)/NS[iv]
      }
    }
          
    iii <- 0
    i.X <- 1:length(nl)
    E.mc <- rep(NA,length(nl)); E.mc.1 <- E.mc; E.mc.2 <- E.mc; E.mc.3 <- E.mc
    E.mc.4 <- E.mc; E.mc.5 <- E.mc; E.mc.6 <- E.mc; E.mc.7 <- E.mc
    q025 <- rep(NA,length(nl))
    q975 <- rep(NA,length(nl))
    for (i in nl) {
      iii <- iii+1
      if (i > 1) {
        e.sums <- colSums(mc.sums[1:i,])
        e.sums.1 <- colSums(mc.sums.1[1:i,])
        e.sums.2 <- colSums(mc.sums.2[1:i,])
        e.sums.3 <- colSums(mc.sums.3[1:i,])
        e.sums.4 <- colSums(mc.sums.4[1:i,])
        e.sums.5 <- colSums(mc.sums.5[1:i,])
        e.sums.6 <- colSums(mc.sums.6[1:i,])
        e.sums.7 <- colSums(mc.sums.7[1:i,])
      } else {
        e.sums   <- 1; e.sums.1 <- 1; e.sums.2 <- 1; e.sums.3 <- 1
        e.sums.4 <- 1; e.sums.5 <- 1; e.sums.6 <- 1; e.sums.7 <- 1       
      }
      E.mc[iii] <- mean(e.sums)
      q025[iii] <- quantile(e.sums,0.025)
      q975[iii] <- quantile(e.sums,0.975)
      
      E.mc.1[iii] <- mean(e.sums.1)
      E.mc.2[iii] <- mean(e.sums.2)
      E.mc.3[iii] <- mean(e.sums.3)
      E.mc.4[iii] <- mean(e.sums.4)
      E.mc.5[iii] <- mean(e.sums.5)
      E.mc.6[iii] <- mean(e.sums.6)
      E.mc.7[iii] <- mean(e.sums.7)
    }         
    upper.lim <- data.frame(y=exp(q975),x=nl)
    lower.lim <- data.frame(y=exp(q025),x=nl)
    best <- data.frame(y=exp(E.mc),x=nl)
    best.1 <- data.frame(y=exp(E.mc.1),x=nl)
    best.2 <- data.frame(y=exp(E.mc.2),x=nl)
    best.3 <- data.frame(y=exp(E.mc.3),x=nl)
    best.4 <- data.frame(y=exp(E.mc.4),x=nl)
    best.5 <- data.frame(y=exp(E.mc.5),x=nl)
    best.6 <- data.frame(y=exp(E.mc.6),x=nl)
    best.7 <- data.frame(y=exp(E.mc.7),x=nl)
    upper.fit <- lm(y ~ x,data=upper.lim)
    lower.fit <- lm(y ~ x,data=lower.lim)
    best.fit <- lm(y ~ x,data=best)
    best.fit.1 <- lm(y ~ x,data=best.1)
    best.fit.2 <- lm(y ~ x,data=best.2)
    best.fit.3 <- lm(y ~ x,data=best.3)
    best.fit.4 <- lm(y ~ x,data=best.4)
    best.fit.5 <- lm(y ~ x,data=best.5)
    best.fit.6 <- lm(y ~ x,data=best.6)
    best.fit.7 <- lm(y ~ x,data=best.7)
    b.u.fit[iv,] <- as.numeric(upper.fit$coefficients)
    b.l.fit[iv,] <- as.numeric(lower.fit$coefficients)
    b.fit[iv,]   <- as.numeric(best.fit$coefficients)
    b.fit.1[iv,]   <- as.numeric(best.fit$coefficients)
    b.fit.2[iv,]   <- as.numeric(best.fit$coefficients)
    b.fit.3[iv,]   <- as.numeric(best.fit$coefficients)
    b.fit.4[iv,]   <- as.numeric(best.fit$coefficients)
    b.fit.5[iv,]   <- as.numeric(best.fit$coefficients)
    b.fit.6[iv,]   <- as.numeric(best.fit$coefficients)
    b.fit.7[iv,]   <- as.numeric(best.fit$coefficients)
  }
  rm(mc.sums,rec, mc.sums.1,rec.1, mc.sums.2,rec.2, mc.sums.3,rec.3,
    mc.sums.4 ,rec.4, mc.sums.5,rec.5, mc.sums.6,rec.6, mc.sums.7,rec.7)

  i.sort <- order(NS); NS <- NS[i.sort]
  b.fit <- b.fit[i.sort,]
  b.u.fit <- b.u.fit[i.sort,]; b.l.fit <- b.l.fit[i.sort,]
  b.fit.1 <- b.fit.1[i.sort,]
  b.fit.2 <- b.fit.2[i.sort,]
  b.fit.3 <- b.fit.3[i.sort,]
  b.fit.4 <- b.fit.4[i.sort,]
  b.fit.5 <- b.fit.5[i.sort,]
  b.fit.6 <- b.fit.6[i.sort,]
  b.fit.7 <- b.fit.7[i.sort,]
  GPC2004=list(b.fit=b.fit,b.fit.1=b.fit.1,b.fit.2=b.fit.2,b.fit.3=b.fit.3,b.fit.4=b.fit.4,b.fit.5=b.fit.5,b.fit.6=b.fit.6,b.fit.7=b.fit.7,b.l.fit=b.l.fit,b.u.fit=b.u.fit,E.mc.1=E.mc.1,E.mc.2=E.mc.2,E.mc.3=E.mc.3,E.mc.4=E.mc.4,E.mc.5=E.mc.5,E.mc.6=E.mc.6,E.mc.7=E.mc.7,NS=NS)
  save(file="GPC2004.rda",GPC2004)
} else { data(GPC2004); attach(GPC2004)}

 x11()
 plot(range(NS),range(c(b.u.fit[,1],b.l.fit[,1]),na.rm=TRUE),type="n",
      main="Intercept for different N")
 grid()
 points(NS,b.u.fit[,1],pch=20,col="red")
 points(NS,b.l.fit[,1],pch=20,col="blue")
 lines(NS,b.u.fit[,1],lty=2,col="red")
 lines(NS,b.l.fit[,1],lty=2,col="blue")

 ord.fit <- lm(log(b.u.fit[,2]) ~ log(NS),weight=NS)
 a <- summary(ord.fit)
 ord.fit <- lm(log(b.l.fit[,2]) ~ log(NS),weight=NS)
 b <- summary(ord.fit)
 b.u.t <- NS^as.numeric(a$coefficients[2])
 b.u <- as.numeric(predict(lm(b.u.fit[,2] ~ b.u.t)))
 b.l.t <- NS^as.numeric(b$coefficients[2])
 b.l <- as.numeric(predict(lm(b.l.fit[,2] ~ b.l.t)))

 x11()
 plot(range(NS),range(c(b.u.fit[,2],b.l.fit[,2])),type="n",
      main="Slope for different N")
 grid()
 points(NS,b.u.fit[,2],pch=20,col="red")
 lines(NS,b.u,col="darkred",lty=3)
 lines(NS,b.l,col="darkblue",lty=3)
 points(NS,b.l.fit[,2],pch=20,col="blue")
 lines(NS,b.u.fit[,2],lty=2,col="red")
 lines(NS,b.l.fit[,2],lty=2,col="blue")

 x11()
 plot(NS,b.fit[,1],type="l",lwd=1,
      main="Intercept for different distributions and N")
 lines(NS,b.fit.1[,1],lwd=5,col="red",lty=1)
 lines(NS,b.fit.2[,1],lwd=4,col="blue",lty=1)
 lines(NS,b.fit.3[,1],lwd=4,col="grey20",lty=2)
 lines(NS,b.fit.4[,1],lwd=3,col="darkred",lty=2)
 lines(NS,b.fit.5[,1],lwd=3,col="darkblue",lty=2)
 lines(NS,b.fit.6[,1],lwd=2,col="darkgreen",lty=2)
 lines(NS,b.fit.7[,1],lwd=2,col="grey60",lty=2)

 x11()
 plot(NS,b.fit[,2],type="l",lwd=1,
      main="Slope for different distributions and N")
 lines(NS,b.fit.1[,2],lwd=5,col="red",lty=1)
 lines(NS,b.fit.2[,2],lwd=4,col="blue",lty=1)
 lines(NS,b.fit.3[,2],lwd=4,col="grey20",lty=2)
 lines(NS,b.fit.4[,2],lwd=3,col="darkred",lty=2)
 lines(NS,b.fit.5[,2],lwd=3,col="darkblue",lty=2)
 lines(NS,b.fit.6[,2],lwd=2,col="darkgreen",lty=2)
 lines(NS,b.fit.7[,2],lwd=2,col="grey60",lty=2)

# Test2 for serial correlation: subsample the data by taking every 3rd month:
#       for spatial correlation: subsample the data by taking fewer stations further apart:

print("Test 2: subsample data to reduce spatial and serial correlation")

JAJO <- c(1,4,7,10); FMAN <- c(2,5,8,11); MJSD <- c(3,6,9,12)
subsample <- c( FMAN +   0,    # Aberdeen: Jan, Apr, Jul, Oct
                FMAN +  12,    # Bismarck
                FMAN +  24,    # Buenos-Aires
                MJSD +  36,    # Concord
                MJSD +  48,    # Funchal
                FMAN +  60,    # Honolulu
                FMAN +  72,    # ishigakijima
                JAJO +  84,    # Lisboa
                JAJO +  96,    # Moskva
                FMAN + 108,    # Nassau
                MJSD + 120,    # Portland
                MJSD + 132,    # Saentis
                MJSD + 144,    # Sao Paulo
                FMAN + 156,    # Seychelles
                FMAN + 168,    # Thessaloniki
                MJSD + 180,    # Thiruvanantha
                FMAN + 192)    # Turuhansk

               
t2m.x <- t2m.all[,subsample]
ns.x <- length(subsample)
recs.1 <- rep(-999,ns.x)
recs.2 <- rep(-999,ns.x)
y.1x <- rep(0,nt)
y.2x <- rep(0,nt)

for (it in 1:nt) {
  ii <- nt-it+1
  vec.1 <- t2m.x[it,]
  vec.2 <- t2m.x[ii,]
  i.rec.1 <- (vec.1 > recs.1) & is.finite(vec.1)
  i.rec.2 <- (vec.2 > recs.2) & is.finite(vec.2)
  if (sum(i.rec.1,na.rm=T)>0) {
    recs.1[i.rec.1] <- vec.1[i.rec.1]
  }
  if (sum(i.rec.2,na.rm=T)>0) {
    recs.2[i.rec.2] <- vec.2[i.rec.2]
  }  
  y.1x[it] <- sum(i.rec.1,na.rm=T)/sum(is.finite(t2m.x[it,]))
  y.2x[it] <- sum(i.rec.2,na.rm=T)/sum(is.finite(t2m.x[ii,]))
}                 

  x11()
  par(col.axis="white",cex.axis=1.5)
  plot(nl,nl,type="n",ylim=c(-7,300),clim=c(0,130),axes=FALSE,
       main="Expected number of record-events",
       xlab="Record length (n)",
       ylab="Number of record-events / Expectated number",
       sub=paste("using temporal-spatial subsampling (N=",length(subsample),
         ") simulation",sep=""),cex=1.2)
  par(col.axis="black")
  axis(1,at=seq(0,nt,by=10), lty = 1, lwd = 0.5)
  yticks <- exp(seq(1,15,by=0.5))
  axis(2,at=yticks, labels=round(log(yticks),1), lty = 1, lwd = 0.5)
  lines(exp(cumsum(1/(1:200))),lwd=2,lty=2)

  for (i in nl) {
    points(i+0.0,exp(sum(y.2x[1:i])),
         pch=20,cex=1.50,col="grey70")
    points(i+0.1,exp(sum(y.2x[1:i]))+0.5,
         pch=20,cex=0.70,col="grey50")
    points(i+0.2,exp(sum(y.2x[1:i]))+1,
         pch=20,cex=0.50,col="grey30")
    points(i+0.3,exp(sum(y.2x[1:i]))+1.5,
         pch=20,cex=0.30,col="black")
      
    points(i+0.0,exp(sum(y.1x[1:i])),
         pch=20,cex=1.50,col="grey30")
    points(i+0.1,exp(sum(y.1x[1:i]))+0.5,
         pch=20,cex=0.70,col="grey50")
    points(i+0.2,exp(sum(y.1x[1:i]))+1,
         pch=20,cex=0.50,col="grey70")
    points(i+0.3,exp(sum(y.1x[1:i]))+1.5,
         pch=20,cex=0.30,col="white")
  }
  i68 <- NS==68
  eb <- b.fit[i68,1] + b.fit[i68,2]*nl
  eq025 <- b.l.fit[i68,1] + b.l.fit[i68,2]*nl
  eq975 <- b.u.fit[i68,1] + b.u.fit[i68,2]*nl
  lines(nl,eq025,lty=3,lwd=0.5,col="grey60")
  lines(nl,eq975,lty=3,lwd=0.5,col="grey60")
  lines(nl,eb,lty=3,lwd=0.5,col="grey40")

  legend(5,297,c("E=cumsum(1/(1:n))",
                 paste("95% C.I. Monte-Carlo:",mc.descr)),
       cex=0.75,lty=c(2,2),lwd=c(2,0.5),pch=c(26,26),
       col=c("black","grey40","grey40"),bg="grey95")

    points(20+0.0,-5+0.00,
         pch=20,cex=1.50,col="grey70")
    points(20+0.1,-5+0.2,
         pch=20,cex=0.70,col="grey50")
    points(20+0.2,-5+0.4,
         pch=20,cex=0.50,col="grey30")
    points(20+0.3,-5+0.6,
         pch=20,cex=0.30,col="black")
    text(30,-5,"Backward")
    
    points(60+0.0,-5+0.00,
         pch=20,cex=1.50,col="grey30")
    points(60+0.1,-5+0.2,
         pch=20,cex=0.70,col="grey50")
    points(60+0.2,-5+0.4,
         pch=20,cex=0.50,col="grey70")
    points(60+0.3,-5+0.6,
         pch=20,cex=0.30,col="white")
    text(50,-5,"Forward")
  

  c1 <- cor(t2m.x,use="pairwise.complete.obs")
  MC <- rep(NA,N.test)
  for (i in 1:N.test) MC[i] <- cor(rnorm(107),rnorm(107))
  CI <- c(quantile(MC,0.025),quantile(MC,0.975))
  c1m <- c1; c1m[(c1m > CI[1]) & (c1m < CI[2])] <- NA
  my.col <- rgb(c(seq(0,1,length=10)^0.5,1,seq(1,0,length=10)^0.5),
                c(seq(0,1,length=10)^0.5,1,seq(1,0,length=10)^0.5),
                c(seq(0,1,length=10)^0.5,1,seq(1,0,length=10)^0.5))

  x11()
  filled.contour(1:68,1:68,c1m,levels=seq(-1,1,by=0.1),col = my.col)
      mar.orig <- (par.orig <- par(c("mar","las","mfrow")))$mar
    on.exit(par(par.orig))

    w <- (3 + mar.orig[2]) * par('csi') * 2.54
    layout(matrix(c(2, 1), nc=2), widths=c(1, lcm(w)))
    
    par(las = 1)
    mar <- mar.orig
    mar[4] <- 1
    par(mar=mar)
  grid()
  contour(1:68,1:68,c1,levels=seq(-1,1,by=0.2),col,add=TRUE)
  title(main="Correlation matrix",sub="Testing for dependencies",
        ylab="series",xlab="series")

 x11()
 plot(nl,E.mc,type="l",lwd=6,xlab="n",ylab="Number of records",
      main="Slope for different distributions and N")
 grid()
 lines(nl,E.mc.1,lwd=5,col="red",lty=1)
 lines(nl,E.mc.2,lwd=4,col="blue",lty=1)
 lines(nl,E.mc.3,lwd=4,col="grey20",lty=2)
 lines(nl,E.mc.4,lwd=3,col="darkred",lty=2)
 lines(nl,E.mc.5,lwd=3,col="darkblue",lty=2)
 lines(nl,E.mc.6,lwd=1,col="green",lty=2)
 lines(nl,E.mc.7,lwd=2,col="grey60",lty=2)

 legend(60,3,c("Normal","gamma s=1","gamma s=10","gamma s=0.1",
               "GEV s=1","GEV s=10","GEV s=-1","Binomial"),
        col=c("black","red","blue","grey20","darkred","darkblue","green","grey60"),
        lwd=c(6,5,4,4,3,4,1,2),lty=c(1,1,1,2,2,2,2,2),bg="grey95")


# Simple plot

  x11()
  par(col.axis="white",cex.axis=1.5)
  plot(nl,nl,type="n",ylim=c(-7,400),xlim=c(0,100),axes=FALSE,
       main="Expected number of record-events",
       xlab="Record length (n)",
       ylab="Number of record-events / Expectated number",
       sub=paste("mean observed vs Monte-Carlo (N=",N.test,
         ") simulation",sep=""),
       cex=1.2)
  par(col.axis="black")
  axis(1,at=seq(0,nt,by=10), lty = 1, lwd = 0.5)
  yticks <- exp(seq(1,15,by=0.5))
  axis(2,at=yticks, labels=round(log(yticks),1), lty = 1, lwd = 0.5)
  lines(exp(cumsum(1/(1:200))),lwd=2,lty=2)

  iii <- 0
  i.X <- 1:length(nl)
  E.mc <- rep(NA,length(nl))
  q025 <- rep(NA,length(nl))
  q975 <- rep(NA,length(nl))
  for (i in nl) {
    iii <- iii+1
      
    points(i+0.0,exp(sum(y.1[1:i])),
         pch=20,cex=1.50,col="grey30")
    points(i+0.1,exp(sum(y.1[1:i]))+0.5,
         pch=20,cex=0.70,col="grey50")
    points(i+0.2,exp(sum(y.1[1:i]))+1,
         pch=20,cex=0.50,col="grey70")
    points(i+0.3,exp(sum(y.1[1:i]))+1.5,
         pch=20,cex=0.30,col="white")

    if (i > 1) e.sums <- colSums(mc.sums[1:i,]) else
               e.sums <- 1
    E.mc[iii] <- mean(e.sums)
    q025[iii] <- quantile(e.sums,0.025)
    q975[iii] <- quantile(e.sums,0.975)
    lines(rep(i,2),exp(c(q025[iii],q975[iii])),lty=1,col="grey70")
  }
  lines(nl,exp(q025),lty=2,lwd=0.5,col="grey40")
  lines(nl,exp(q975),lty=2,lwd=0.5,col="grey40")

  upper.lim <- data.frame(y=exp(q975),x=nl)
  lower.lim <- data.frame(y=exp(q025),x=nl)
  best <- data.frame(y=exp(E.mc),x=nl)

  upper.fit <- lm(y ~ x,data=upper.lim)
  lower.fit <- lm(y ~ x,data=lower.lim)
  best.fit <- lm(y ~ x,data=best)

  abline(upper.fit,lty=3,lwd=0.5,col="grey60")
  abline(lower.fit,lty=3,lwd=0.5,col="grey60")
  abline(best.fit,lty=3,lwd=0.5,col="grey40")
  lines(nl,exp(E.mc),lwd=0.5,lty=1,col="blue")
  points(nl,exp(E.mc),pch=19,col="blue",cex=0.5)
  lines(exp(cumsum(1/(1:200))),lwd=2,lty=2)

  legend(5,397,c("E=cumsum(1/(1:n))",
                 paste("mean E Monte-Carlo:",mc.descr),
                 paste("95% C.I. Monte-Carlo:",mc.descr),
                 "T(2m) from NASA/GISS"),
       cex=0.75,lty=c(2,1,2,0),lwd=c(2,0.5,0.5),pch=c(26,19,26,19),
       col=c("black","blue","grey40","grey30"),bg="grey95")
    
    points(8+0.0,345+0.00,
         pch=20,cex=1.50,col="grey30")
    points(8+0.1,345+0.2,
         pch=20,cex=0.70,col="grey50")
    points(8+0.2,345+0.4,
         pch=20,cex=0.50,col="grey70")
    points(8+0.3,345+0.6,
         pch=20,cex=0.30,col="white")
  
