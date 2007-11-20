
par(col.axis="white",lab=c(1,1,1))

t<- seq(1,1000,by=1)
X <- rnorm(length(t))
e1 <- n.records(X)
plot(X,xlim=c(-50,max(t)+50),ylim=c(-2,30),type="l",col="grey",xlab="",ylab="",
     main=paste("n=",length(t)))

polygon(c(-50,rep(max(t)+50,2),-50,-50),c(-3,-3,3,3,-3),col="grey95",border="grey95")
polygon(c(-50,rep(max(t)+50,2),-50,-50),c(14,14,22,22,14),col="grey95",border="grey95")

lines(X,col="grey40")
points(t[e1$events],X[e1$events],pch=20,cex=1.75,col="black")
points(t[e1$events],X[e1$events]+0.025,pch=20,cex=1.5,col="grey30")
points(t[e1$events]+0.05,X[e1$events]+0.025,pch=20,cex=0.7,col="grey50")
points(t[e1$events]+0.07,X[e1$events]+0.025,pch=20,cex=0.5,col="grey70")
points(t[e1$events]+0.1,X[e1$events]+0.025,pch=20,cex=0.3,col="white")
text(-50,0,"A",cex=1.5,font=2)
text(max(t)+50,0,sum(e1$events),cex=1.2,font=2)

X <- rnorm(length(t),mean=0.0055*t)+6
e1 <- n.records(X)
lines(X,col="grey40")
points(t[e1$events],X[e1$events],pch=20,cex=1.75,col="black")
points(t[e1$events],X[e1$events]+0.025,pch=20,cex=1.5,col="grey30")
points(t[e1$events]+0.05,X[e1$events]+0.025,pch=20,cex=0.7,col="grey50")
points(t[e1$events]+0.07,X[e1$events]+0.025,pch=20,cex=0.5,col="grey70")
points(t[e1$events]+0.1,X[e1$events]+0.025,pch=20,cex=0.3,col="white")
text(-50,7,"B",cex=1.5,font=2)
text(max(t)+50,7,sum(e1$events),cex=1.2,font=2)

X <- rnorm(length(t),mean=-0.0045*t)+20
e1 <- n.records(X)
lines(X,col="grey40")
points(t[e1$events],X[e1$events],pch=20,cex=1.75,col="black")
points(t[e1$events],X[e1$events]+0.025,pch=20,cex=1.5,col="grey30")
points(t[e1$events]+0.05,X[e1$events]+0.025,pch=20,cex=0.7,col="grey50")
points(t[e1$events]+0.07,X[e1$events]+0.025,pch=20,cex=0.5,col="grey70")
points(t[e1$events]+0.1,X[e1$events]+0.025,pch=20,cex=0.3,col="white")
text(-50,15,"C",cex=1.5,font=2)
text(max(t)+50,15,sum(e1$events),cex=1.2,font=2)

X <- rnorm(length(t),sd=0.0017*t)+25
e1 <- n.records(X)
lines(X,col="grey40")
points(t[e1$events],X[e1$events],pch=20,cex=1.75,col="black")
points(t[e1$events],X[e1$events]+0.025,pch=20,cex=1.5,col="grey30")
points(t[e1$events]+0.05,X[e1$events]+0.025,pch=20,cex=0.7,col="grey50")
points(t[e1$events]+0.07,X[e1$events]+0.025,pch=20,cex=0.5,col="grey70")
points(t[e1$events]+0.1,X[e1$events]+0.025,pch=20,cex=0.3,col="white")
text(-50,25,"D",cex=1.5,font=2)
text(max(t)+50,25,sum(e1$events),cex=1.2,font=2)

newFig()
test.iid.test()

newFig()
x <- seq(-5,6,length=400)
plot(x,dnorm(x,sd=1,mean=0),type="l",lwd=3,xlab="",ylab="",main="pdf")
lines(x,dnorm(x,sd=1.5,mean=0),type="l",lwd=2,col="grey70")
lines(x, dnorm(x,sd=1,mean=1),type="l",lwd=2,col="grey50")

text(0,0.37,"A",cex=1.5)
text(1,0.37,"B",cex=1.5,col="grey50")
text(0,0.24,"D",cex=1.5,col="grey70")
