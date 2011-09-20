
source("Lyapunov-approx-fns.R")

# for exploration:
dvals <- 2^((-20:30)/5)

png(file="four-cases.png")
par(mfrow=c(2,2))

# example when no dispersal is best:
mu <- c(1,10)
sigma <- c(4,4)
Dmat <- rbind( c(-0.1,0.1), c(1,-1) )
exact <- sapply(dvals, function(d) { twopatch(mu,sigma,d*Dmat)$lyap } )
ptitle <- sprintf("mu=(%.0f,%.0f) sigma=(%.0f,%.0f), Dmat=(%.1f,%.1f)",mu[1],mu[2],sigma[1],sigma[2],Dmat[1,2],Dmat[2,1])
plot(dvals,exact,main=ptitle,cex.main=1)

# example with intermediate dispersal best:
mu <- c(1,10)
sigma <- c(3,5)
Dmat <- rbind( c(-0.1,0.1), c(1,-1) )
exact <- sapply(dvals, function(d) { twopatch(mu,sigma,d*Dmat)$lyap } )
ptitle <- sprintf("mu=(%.0f,%.0f) sigma=(%.0f,%.0f), Dmat=(%.1f,%.1f)",mu[1],mu[2],sigma[1],sigma[2],Dmat[1,2],Dmat[2,1])
plot(dvals,exact,main=ptitle,cex.main=1)

# example with maximal dispersal best:
mu <- c(8,8)
sigma <- c(1,10)
Dmat <- rbind( c(-0.1,0.1), c(1,-1) )
exact <- sapply(dvals, function(d) { twopatch(mu,sigma,d*Dmat)$lyap } )
ptitle <- sprintf("mu=(%.0f,%.0f) sigma=(%.0f,%.0f), Dmat=(%.1f,%.1f)",mu[1],mu[2],sigma[1],sigma[2],Dmat[1,2],Dmat[2,1])
plot(dvals,exact,main=ptitle,cex.main=1)

# example with intermediate dispersal worst:
mu <- c(15,1)
sigma <- c(1.5,10)
Dmat <- rbind( c(-0.1,0.1), c(.2,-.2) )
exact <- sapply(dvals, function(d) { twopatch(mu,sigma,d*Dmat)$lyap } )
ptitle <- sprintf("mu=(%.0f,%.0f) sigma=(%.0f,%.0f), Dmat=(%.1f,%.1f)",mu[1],mu[2],sigma[1],sigma[2],Dmat[1,2],Dmat[2,1])
plot(dvals,exact,main=ptitle,cex.main=1)

dev.off()

# interpolate through both intermediate dispersal best and worst
png(file="interpolation.png")
interps <- sapply( (0:8)/8, function(a) {
        mu <- c(15,1)
        sigma <- a*c(6,1) + (1-a)*c(1,6)
        Dmat <- rbind( c(-0.1,0.1), c(1,-1) )
        sapply(dvals, function(d) { twopatch(mu,sigma,d*Dmat)$lyap } )
    } )
par(mfrow=c(3,3),mai=c(0,0.25,0,0))
apply(interps,2,function (x) plot(dvals,x,xaxt="n"))
dev.off()

