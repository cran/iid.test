library(clim.pact)

n.records <- function(y) {
  m <- length(y)
  y[!is.finite(y)] <- min(y,na.rm=TRUE)
  N <- 1; N.rev <- N
  t <- rep(1,m); t.rev <- rep(m,m)
  events <- rep(FALSE,m); events.rev <- events
  events[1] <- TRUE; events.rev[m] <- TRUE
  for (i in 2:m) {
    if (y[i] > max(y[1:(i-1)],na.rm=TRUE)) {
      N <- N + 1
      t[N] <- i
      events[i] <- TRUE
    }
    if (y[m-i+1] > max(y[(m-i+2):m],na.rm=TRUE)) {
      N.rev <- N.rev + 1
      t.rev[N.rev] <- length(y)-i
      events.rev[length(y)-i] <- TRUE
    }
  }
  t <- t[1:N]; t.rev <- t.rev[1:N.rev]
  records <- list(N=N,t=t,events=events,N.rev=N.rev, 
                  t.rev=t.rev, events.rev=events.rev)
  invisible(records)
} 

iid.test <- function(Y,plot=TRUE,Monte.Carlo=TRUE,N.test=200) {
  Y <- as.matrix(Y)
  Y[!is.finite(Y)] <- NA
  t.r <- dim(Y)
  events <- matrix(rep(FALSE,t.r[1]*t.r[2]),t.r[1],t.r[2])
  events.rev <- events
  N.records <- rep(NA,t.r[2])

  if (plot) {
    par(col.axis="white")
    plot(c(1,t.r[1]),c(1,2*t.r[2]),type="n",main="iid-test",
         xlab="time",ylab="location")
    par(col.axis="black")
    axis(1)
    lines(c(1,t.r[1]),rep(t.r[2],2),lwd=3)
    par.0 <- par(); par(srt=90)
    text(0,round(t.r[2]/2),"Forward",cex=1,vfont=c("sans serif","italic"))
    text(0,round(3*t.r[2]/2),"Backward",cex=1,vfont=c("sans serif","italic"))
    par(par.0)
  }

  for (ir in 1:t.r[2]) {
    record.stats <- n.records(Y[,ir])
    N.records[ir] <- record.stats$N
    events[,ir] <- as.numeric(record.stats$events)
    events.rev[,ir] <- as.numeric(record.stats$events.rev)
    if (plot) {
      lines(c(1,t.r[1]),rep(ir,2),col="grey70")
      points(record.stats$t+0.1,rep(ir,record.stats$N)+0.1,pch=20,cex=1.50,col="grey30")
      points(record.stats$t+0.1,rep(ir,record.stats$N)+0.2,pch=20,cex=0.70,col="grey50")
      points(record.stats$t+0.1,rep(ir,record.stats$N)+0.3,pch=20,cex=0.50,col="grey70")
      points(record.stats$t+0.1,rep(ir,record.stats$N)+0.4,pch=20,cex=0.30,col="white")

      lines(c(1,t.r[1]),rep(ir+t.r[2],2),col="grey70")
      points(record.stats$t.rev+0.1,rep(ir,record.stats$N.rev)+0.1+t.r[2],
      pch=20,cex=1.50,col="grey30")
      points(record.stats$t.rev+0.1,rep(ir,record.stats$N.rev)+0.2+t.r[2],
      pch=20,cex=0.70,col="grey50")
      points(record.stats$t.rev+0.1,rep(ir,record.stats$N.rev)+0.3+t.r[2],
      pch=20,cex=0.50,col="grey70")
      points(record.stats$t.rev+0.1,rep(ir,record.stats$N.rev)+0.4+t.r[2],
      pch=20,cex=0.30,col="white")
    }
  }

  events[!is.finite(Y)] <- NA
  events.rev[!is.finite(Y)] <- NA
  record.density <- rowMeans(events,na.rm=TRUE)
  record.density.rev <- reverse(rowMeans(events.rev,na.rm=TRUE))
  N <- length(record.density)

  if (Monte.Carlo) {
    print(paste("Please be patient -",N.test,"Monte Carlo runs in progress..."))
    record.mc <- rep(NA,2*N.test*N); dim(record.mc) <- c(N,N.test,2)
    for (ii in 1:N.test) {
      mc.stats <- test.iid.test(d=dim(events),plot=FALSE,Monte.Carlo=FALSE)  
      record.mc[,ii,1] <- cumsum(mc.stats$record.density)
      record.mc[,ii,2] <- cumsum(mc.stats$record.density.rev)
    } 

    q025=rep(NA,N); q975=q025    
    for (i in 1:N) {
      q025[i] <- quantile(record.mc[i,,],0.025)
      q975[i] <- quantile(record.mc[i,,],0.955)
    }
    sub <- paste("Shaded region= 95% conf.int. from Monte-Carlo with N=",N.test)
  } else sub <- ""

  if (plot) {
    newFig()
    par(col.axis="white")
    Time <- 1:N
    plot(Time,exp( cumsum( 1/(1:N)) ),type="l",lwd=3,col="grey60",
           xlab="Time",main="Observed & Expected number of record-events",
           sub=sub)
    par(col.axis="black")
    axis(1)
    axis(2,at=exp(1:(2*sum(record.density,na.rm=TRUE))),
         label=1:(2*sum(record.density,na.rm=TRUE)))
    legend(1,exp(sum(1/(1:N))),c("Theoretical","Forward","Backward"),
           pch=c(26,20,21),lwd=c(3,0,0),lty=c(1,0,0),col=c("grey60",rep("black",2)))

    if (Monte.Carlo) {
      polygon(c(Time,reverse(Time)),c(exp(q025),reverse(exp(q975))),
              col="grey90",border="grey85")
      lines(Time,exp( cumsum( 1/(1:N)) ),lwd=3,col="grey60")
    }
    grid()
    points(Time,exp(cumsum(record.density)),pch=20,cex=0.9)
    points(Time,exp(cumsum(record.density.rev)),pch=21,cex=0.9)
  }

  results <- list(record.density=record.density,
                  record.density.rev=record.density.rev)
  invisible(results)
}

test.iid.test <- function(distr="rnorm",d=c(100,30),plot=TRUE,Monte.Carlo=TRUE) {
  rnd <- eval(parse(text=paste(distr,"(",d[1]*d[2],")",sep="")))
  dim(rnd) <- c(d[1],d[2])
  test.results <- iid.test(rnd,plot=plot,Monte.Carlo=Monte.Carlo)
  invisible(test.results)
}

daily.station.records <- function(obs,element="precip",subsample=5,tolerance=2) {

  if (class(obs)[2] != "daily.station.record") 
     stop("Need a 'daily.station.record' object!")
  years <- as.numeric(rownames(table(obs$yy)))
  ny <- length(years)
  dat <- rep(NA,ny*366); dim(dat) <- c(ny,366)

  for (i in 1:ny) {
    iyear <- is.element(obs$yy,years[i])
    ii <- julday(obs$mm[iyear],obs$dd[iyear],obs$yy[iyear]) -
          julday(1,1,years[i])+1
    data.thisyear <- eval(parse(text=paste("obs$",element,"[iyear]",sep="")))
    print(c(years[i],sum(iyear),NA,range(ii),NA,length(data.thisyear)))
  
    if (sum(is.finite(data.thisyear)) > 2) plot(data.thisyear)
    dat[i,ii] <- data.thisyear
  }
  
  #plot(dat[1,],main=obs$location,ylab=element,xlab="Day in the year",
  #     pch=20,cex=0.8,col="grey50")
  #for (i in 1:ny) points(dat[i,],pch=20,cex=0.8,col="grey50")

  print(paste("sub-sample every",subsample,"points"))

  ykeep <- rep(TRUE,366)
  for (i in 1:length(ykeep)) {
    if (mod(i,subsample) != 1) ykeep[i] <- FALSE
    if (sum(!is.finite(dat[,i])) > tolerance) ykeep[i] <- FALSE
  }
  dat <- dat[,ykeep] # extra days during leap year will bias the results 

  #image(log(dat)); x11()


  newFig()
  iid.test(dat)
}