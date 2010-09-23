rm(list=ls())

# Load nesescary packages:

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
isel <- c(1,6,8,10,16,18,19,23,24,25,28,29,30,32,35,36,38)
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
p.1 <-1/(1:nt)
p.N <- 1 - (1-1/(1:nt))^ns

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

               
t2m.all <- t2m.all[,subsample]
ns <- length(subsample)

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
text(-8,30,"Forward",cex=2,vfont=c("sans serif","italic"))
text(-8,120,"Backward",cex=2,vfont=c("sans serif","italic"))
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

x11()
plot(p.1,type="l",col="black",lwd=3,
     main="Probability of new records",
     sub="Different number of independent series",
     ylab="Record density",xlab="Time",cex=1.5)
lines(p.N,col="grey40",lwd=3,lty=2)
grid()
legend(85,1,c("E_n/N",paste("N=",ns)),
       col=c("black","grey40"),
       lwd=c(3,3),lty=c(1,2),bg="grey95",cex=0.75)

t2m.x <- t2m.all
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
NS <- c(204,10,20,30,50,68,75,100,500,600,1000)
nl <- seq(1,nt,by=5)

data(GPC2004); attach(GPC2004)


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
  lines(nl,eq025,lty=2,lwd=1,col="grey40")
  lines(nl,eq975,lty=2,lwd=1,col="grey40")
  lines(nl,eb,lty=2,lwd=0.5,col="grey40")
  mc.descr<- 'white noise'

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
  
