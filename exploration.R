
source("Lyapunov-approx-fns.R")
h <- .01
tmax <- 100  # should be at least 1/sqrt(h)
nsamps <- 200

# compute exact and approximate values of chi
# as a function of 1/delta


# for two-patch case
idvals <- (1:20)/(2*80)
mu <- c(1,10)
sigma <- c(4,4)
Dmat <- rbind(c(-0.1,0.1),c(1,-1))
# exact two-patch values
exact <- sapply(idvals, function (id) { twopatch(mu,sigma,(1/id)*Dmat)$lyap } )
# high-dispersal approximation
approxchi <- hdlyap(mu=mu, Gmat=diag(sigma), Q=Dmat)
# by simulation, now
simchi <- lapply( idvals, function (id) { lyapprox( mu, Dmat/id, G=diag(sigma), h, tmax ) } )
simchis <- sapply( simchi, function (x) { c( x$lyap, x$lyapsd ) } )
plot( idvals, exact, xlab=expression(1/delta), ylab=expression(chi) )
lines( idvals, approxchi$a - idvals*approxchi$b )
lines( idvals, simchis[1,], col="red")
lines( idvals, simchis[1,]+(-2)*simchis[2,], col="red", lty=2)
lines( idvals, simchis[1,]+2*simchis[2,], col="red", lty=2)
# save( idvals, mu, sigma, Dmat, exact, simchi, simchis, approxchi, h, tmax, file="sims1.RData" )
# load("sims1.RData")

# hm.  what effect has h?
#  -> no standard effect.
exchi <- twopatch(mu, sigma, Dmat)$lyap
schi <- lapply( c(1, .1, .01, .001), function (h) { lyapprox(mu, Dmat, G=diag(sigma), h, 10000*h) } )
sschi <- sapply(schi, function (x) { c(x$lyap, x$lyapsd) } )
sschi <- rbind(sschi, sschi[1,]+2*sschi[2,], sschi[1,]-2*sschi[2,])
plot( 1:4, sschi[1,], ylim=range(sschi[1,],exchi) )
abline(h=exchi)
lines(1:4, sschi[3,], lty=2)
lines(1:4, sschi[4,], lty=2)

# for more patches
# idvals <- (1:80)/(80)
idvals <- c( (1:20)/160, (11:40)/80, (11:40)/20 )
npatches <- 8 
badprop <- .75  # proportion of bad patches
nbad <- floor(badprop*npatches)
mu <- c( rep(1,nbad), rep(10,npatches-nbad) )
Gmat <- diag(rep(4,npatches))
Dmat <- matrix( 1/10, nrow=npatches, ncol=npatches )
Dmat[(nbad+1):npatches,] <- 1
diag(Dmat) <- - apply(Dmat-diag(diag(Dmat)),1,sum)
# exact two-patch values  (not so useful?)
twoDmat <- matrix( c( -sum(Dmat[1,(nbad+1):npatches]), sum(Dmat[npatches,1:nbad]), sum(Dmat[1,(nbad+1):npatches]), -sum(Dmat[npatches,1:nbad]) ), nrow=2 )
twoDmat <- nbad*twoDmat/twoDmat[2,1]
# twosig <- diag(Gmat)[c(1,length(mu))]^2 / c(nbad, npatches-nbad); 
twosig <- diag(Gmat)[c(1,length(mu))]
twomu <- mu[c(1,length(mu))]
exact <- sapply(idvals, function (id) { twopatch(twomu,twosig,(1/id)*twoDmat)$lyap } )
# by simulation, now
# simchi <- lyapprox( mu, Dmat, Gmat, h, tmax )
# plotsim(simchi)  # look at it
simchi  <- lapply(idvals, function (id) { 
            # set.seed(1)   # cheating =)
            lyapprox(mu, (1/id)*Dmat, Gmat, h, tmax, dsamp=1) 
        } )
simchis <- sapply(simchi, function (x) { c(x$lyap, x$lyapsd) } )
# high-dispersal approximation
approxchi <- hdlyap(mu=mu, Gmat=Gmat, Q=Dmat)
# save( idvals, npatches, mu, Gmat, Dmat, exact, simchi, simchis, approxchi, h, tmax, file="sims2.RData" ) # with idvals = (1:80)/80
# save( idvals, npatches, mu, Gmat, Dmat, exact, simchi, simchis, approxchi, h, tmax, file="sims3.RData" ) # with wider (log) range of idvals
# save( idvals, npatches, mu, Gmat, Dmat, exact, simchi, simchis, approxchi, h, tmax, file="sims-cheating.RData" ) # set seed same every time

# produce plots
# load("sims-cheating.RData")
# pdf(file="eight-patch-same-seed.pdf", width=6.5, height=5)
load("sims3.RData")
pdf(file="eight-patch.pdf", width=6.5, height=5)
plot( idvals, simchis[1,], xlab=expression(paste("Inverse dispersal rate", 1/delta)), ylab=expression(paste("Lyapunov exponent ", chi)), ylim=range(zeroinf(mu,Dmat,Gmat),simchis[1,]+2*simchis[2,],simchis[1,]-2*simchis[2,]))
lines(lowess( idvals, simchis[1,], f=1/4 ) )
lines( lowess(idvals, simchis[1,]+(-2)*simchis[2,], f=1/4), lty=2)
lines( lowess(idvals, simchis[1,]+2*simchis[2,], f=1/4), lty=2)
# lines( idvals, exact, col="blue" )
abline( approxchi$a, -approxchi$b, col="red" )
lines( c(1.5*mean(range(idvals)), 1.5*max(idvals)), rep(zeroinf(mu,Dmat,Gmat)["zero"],2), lty=3, lwd=1.5 )
text( 1.5*mean(range(idvals))+0.1, zeroinf(mu,Dmat,Gmat)["zero"], expression(delta == 0), adj=c(0,-.25) )
lines( c(0,mean(range(idvals))/3), rep(zeroinf(mu,Dmat,Gmat)["infty"],2), lty=3, lwd=1.5 )
text( 0.2, zeroinf(mu,Dmat,Gmat)["infty"], expression(delta == infinity), adj=c(0,-.25) )
legend("topleft", lty=c(1,1,1,3), col=c("black","blue","red","black"), legend=c("Simulation, 8 patches", "Exact value, 2 patches", "High-dispersal approximation","Asymptotic values"), cex=.75)
dev.off()

# produce plots
# load("sims4.RData")
load("sims5.RData")
pdf(file="fourty-patch-misc.pdf", width=6.5, height=5)
plot( idvals, simchis[1,], xlab=expression(paste("Inverse dispersal rate", 1/delta)), ylab=expression(paste("Lyapunov exponent ", chi)), ylim=range(zeroinf(mu,Dmat,Gmat),simchis[1,]+2*simchis[2,],simchis[1,]-2*simchis[2,]), log="x")
lowline <- lowess( log(idvals), simchis[1,], f=1/3 ); lines( exp(lowline$x), lowline$y )
lowline <- lowess( log(idvals), simchis[1,]+(-2)*simchis[2,], f=1/3); lines( exp(lowline$x), lowline$y, lty=2)
lowline <- lowess( log(idvals), simchis[1,]+2*simchis[2,], f=1/3); lines( exp(lowline$x), lowline$y, lty=2)
lines(idvals, approxchi$a-idvals*approxchi$b, col="red")
# lines( idvals, exact, col="blue" )
lines( c(10^5, 1.5*max(idvals)), rep(zeroinf(mu,Dmat,Gmat)["zero"],2), lty=3, lwd=1.5 )
text( 10^5, zeroinf(mu,Dmat,Gmat)["zero"], expression(delta == 0), adj=c(0,-.25) )
lines( c(min(idvals)/10,1), rep(zeroinf(mu,Dmat,Gmat)["infty"],2), lty=3, lwd=1.5 )
text( .55, zeroinf(mu,Dmat,Gmat)["infty"], expression(delta == infinity), adj=c(0,-.25) )
legend("topright", lty=c(1,1,1,3), col=c("black","blue","red","black"), legend=c("Simulation, 40 patches", "Exact value, 2 patches", "High-dispersal approximation","Asymptotic values"), cex=.75)
dev.off()


# produce plots: 8 patches, log scale
load("sims6.RData")
pdf(file="eight-patch-log.pdf", width=6.5, height=5)
plot( idvals, simchis[1,], xlab=expression(paste("Inverse dispersal rate", 1/delta)), ylab=expression(paste("Lyapunov exponent ", chi)), ylim=range(zeroinf(mu,Dmat,Gmat),simchis[1,]+2*simchis[2,],simchis[1,]-2*simchis[2,]), log="x")
lowline <- lowess( log(idvals), simchis[1,], f=1/3 ); lines( exp(lowline$x), lowline$y )
lowline <- lowess( log(idvals), simchis[1,]+(-2)*simchis[2,], f=1/3); lines( exp(lowline$x), lowline$y, lty=2)
lowline <- lowess( log(idvals), simchis[1,]+2*simchis[2,], f=1/3); lines( exp(lowline$x), lowline$y, lty=2)
lines(idvals, approxchi$a-idvals*approxchi$b, col="red")
# lines( idvals, exact, col="blue" )
lines( c(10^5, 1.5*max(idvals)), rep(zeroinf(mu,Dmat,Gmat)["zero"],2), lty=3, lwd=1.5 )
text( 10^5, zeroinf(mu,Dmat,Gmat)["zero"], expression(delta == 0), adj=c(0,-.25) )
lines( c(min(idvals)/10,1), rep(zeroinf(mu,Dmat,Gmat)["infty"],2), lty=3, lwd=1.5 )
text( .55, zeroinf(mu,Dmat,Gmat)["infty"], expression(delta == infinity), adj=c(0,-.25) )
legend("topright", lty=c(1,1,1,3), col=c("black","blue","red","black"), legend=c("Simulation, 8 patches", "Exact value, 2 patches", "High-dispersal approximation","Asymptotic values"), cex=.75)
dev.off()

# produce plots
load("sims7.RData")
pdf(file="fourty-patch-log.pdf", width=6.5, height=5)
plot( idvals, simchis[1,], xlab=expression(paste("Inverse dispersal rate", 1/delta)), ylab=expression(paste("Lyapunov exponent ", chi)), ylim=range(zeroinf(mu,Dmat,Gmat),simchis[1,]+2*simchis[2,],simchis[1,]-2*simchis[2,]), log="x")
lowline <- lowess( log(idvals), simchis[1,], f=1/3 ); lines( exp(lowline$x), lowline$y )
lowline <- lowess( log(idvals), simchis[1,]+(-2)*simchis[2,], f=1/3); lines( exp(lowline$x), lowline$y, lty=2)
lowline <- lowess( log(idvals), simchis[1,]+2*simchis[2,], f=1/3); lines( exp(lowline$x), lowline$y, lty=2)
lines(idvals, approxchi$a-idvals*approxchi$b, col="red")
# lines( idvals, exact, col="blue" )
lines( c(10^5, 1.5*max(idvals)), rep(zeroinf(mu,Dmat,Gmat)["zero"],2), lty=3, lwd=1.5 )
text( 10^5, zeroinf(mu,Dmat,Gmat)["zero"], expression(delta == 0), adj=c(0,-.25) )
lines( c(min(idvals)/10,1), rep(zeroinf(mu,Dmat,Gmat)["infty"],2), lty=3, lwd=1.5 )
text( .55, zeroinf(mu,Dmat,Gmat)["infty"], expression(delta == infinity), adj=c(0,-.25) )
legend("topright", lty=c(1,1,1,3), col=c("black","blue","red","black"), legend=c("Simulation, 40 patches", "Exact value, 2 patches", "High-dispersal approximation","Asymptotic values"), cex=.75)
dev.off()

# produce plots: 8 patches, log scale, adaptive
load("sims8.RData")
pdf(file="eight-patch-adaptive.pdf", width=6.5, height=5)
plot( idvals, simchis[1,], xlab=expression(paste("Inverse dispersal rate", 1/delta)), ylab=expression(paste("Lyapunov exponent ", chi)), ylim=range(zeroinf(mu,Dmat,Gmat),simchis[1,]+2*simchis[2,],simchis[1,]-2*simchis[2,]), log="x")
lowline <- lowess( log(idvals), simchis[1,], f=1/3 ); lines( exp(lowline$x), lowline$y )
lowline <- lowess( log(idvals), simchis[1,]+(-2)*simchis[2,], f=1/3); lines( exp(lowline$x), lowline$y, lty=2)
lowline <- lowess( log(idvals), simchis[1,]+2*simchis[2,], f=1/3); lines( exp(lowline$x), lowline$y, lty=2)
lines(idvals, approxchi$a-idvals*approxchi$b, col="red")
# lines( idvals, exact, col="blue" )
lines( c(10^5, 1.5*max(idvals)), rep(zeroinf(mu,Dmat,Gmat)["zero"],2), lty=3, lwd=1.5 )
text( 10^5, zeroinf(mu,Dmat,Gmat)["zero"], expression(delta == 0), adj=c(0,-.25) )
lines( c(min(idvals)/10,1), rep(zeroinf(mu,Dmat,Gmat)["infty"],2), lty=3, lwd=1.5 )
text( .55, zeroinf(mu,Dmat,Gmat)["infty"], expression(delta == infinity), adj=c(0,-.25) )
legend("topright", lty=c(1,1,1,3), col=c("black","blue","red","black"), legend=c("Simulation, 8 patches", "Exact value, 2 patches", "High-dispersal approximation","Asymptotic values"), cex=.75)
dev.off()

# produce plots
load("sims9.RData")
pdf(file="fourty-patch-adaptive.pdf", width=6.5, height=5)
plot( idvals, simchis[1,], xlab=expression(paste("Inverse dispersal rate", 1/delta)), ylab=expression(paste("Lyapunov exponent ", chi)), ylim=range(zeroinf(mu,Dmat,Gmat),simchis[1,]+2*simchis[2,],simchis[1,]-2*simchis[2,]), log="x")
lowline <- lowess( log(idvals), simchis[1,], f=1/3 ); lines( exp(lowline$x), lowline$y )
lowline <- lowess( log(idvals), simchis[1,]+(-2)*simchis[2,], f=1/3); lines( exp(lowline$x), lowline$y, lty=2)
lowline <- lowess( log(idvals), simchis[1,]+2*simchis[2,], f=1/3); lines( exp(lowline$x), lowline$y, lty=2)
lines(idvals, approxchi$a-idvals*approxchi$b, col="red")
# lines( idvals, exact, col="blue" )
lines( c(10^5, 1.5*max(idvals)), rep(zeroinf(mu,Dmat,Gmat)["zero"],2), lty=3, lwd=1.5 )
text( 10^5, zeroinf(mu,Dmat,Gmat)["zero"], expression(delta == 0), adj=c(0,-.25) )
lines( c(min(idvals)/10,1), rep(zeroinf(mu,Dmat,Gmat)["infty"],2), lty=3, lwd=1.5 )
text( .55, zeroinf(mu,Dmat,Gmat)["infty"], expression(delta == infinity), adj=c(0,-.25) )
legend("topright", lty=c(1,1,1,3), col=c("black","blue","red","black"), legend=c("Simulation, 40 patches", "Exact value, 2 patches", "High-dispersal approximation","Asymptotic values"), cex=.75)
dev.off()


# produce plots
pdf(file="combined-patches-smooth.pdf", width=7, height=4, pointsize=14)
par(mai=par("mai")[c(1,2,3,2)]*c(1,1,1/2,1))
## fourty patches
load("sims-40patch-log-smooth.RData")  # eight patches smooth
plot( 1/idvals, simchis[1,], xlab=expression(paste("Dispersal rate ", delta)), ylab=expression(paste("Lyapunov exponent ", chi)), type="l", ylim=range(zeroinf(mu,Dmat,Gmat),simchis[1,]+2*simchis[2,],simchis[1,]-2*simchis[2,]), log="x", lwd=2)
lines( 1/idvals, simchis[1,]+(-2)*simchis[2,], lty=2, lwd=2)
lines( 1/idvals, simchis[1,]+(+2)*simchis[2,], lty=2, lwd=2)
lines(1/idvals, approxchi$a-idvals*approxchi$b, col="red", lwd=2)
# lines( c(.5*min(1/idvals),10^-2), rep(zeroinf(mu,Dmat,Gmat)["zero"],2), lty=3, lwd=1.5 )
# text( 10^-2, zeroinf(mu,Dmat,Gmat)["zero"], expression(delta == 0), adj=c(0,-.25) )
axis(side=4, at=zeroinf(mu,Dmat,Gmat)["zero"], labels=expression(delta == 0), las=1)
# lines( c(10^-2,1.5*max(1/idvals)), rep(zeroinf(mu,Dmat,Gmat)["infty"],2), lty=3, lwd=1.5 )
# text( 10^-2, zeroinf(mu,Dmat,Gmat)["infty"], expression(delta == infinity), adj=c(0,-.25) )
axis(side=4, at=zeroinf(mu,Dmat,Gmat)["infty"], labels=expression(delta[40] == infinity), las=1)
## eight patches
load("sims-8patch-log-smooth.RData")  # eight patches smooth
lines( 1/idvals, simchis[1,], col="blue", lwd=2 )
lines( 1/idvals, simchis[1,]+(-2)*simchis[2,], lty=2, col="blue", lwd=2)
lines( 1/idvals, simchis[1,]+(+2)*simchis[2,], lty=2, col="blue", lwd=2)
lines(1/idvals, approxchi$a-idvals*approxchi$b, col="red", lwd=2)
# lines( c(10^-2,1.5*max(1/idvals)), rep(zeroinf(mu,Dmat,Gmat)["infty"],2), lty=3, lwd=1.5 )
# text( 10^-2, zeroinf(mu,Dmat,Gmat)["infty"], expression(delta == infinity), adj=c(0,-.25) )
axis(side=4, at=zeroinf(mu,Dmat,Gmat)["infty"], labels=expression(delta[8] == infinity), las=1)
## legend
legend("topleft", lty=c(1,1,1,2), col=c("black","blue","red","black"), legend=c("40 patches", "8 patches", "High dispersal", "+/- 2 SE"), cex=.75, lwd=2)
dev.off()

# produce plots
pdf(file="combined-patches-smooth-mal.pdf", width=7, height=4, pointsize=14)
par(mai=par("mai")[c(1,2,3,2)]*c(1,1,1/2,1))
## fourty patches
load("sims-40patch-log-smooth-mal.RData")  # eight patches smooth
plot( 1/idvals, simchis[1,], xlab=expression(paste("Dispersal rate ", delta)), ylab=expression(paste("Lyapunov exponent ", chi)), ylim=c(-1.5,5.5), log="x", type="l", lwd=2)
lowline <- lines( 1/idvals, simchis[1,]+(-2)*simchis[2,], lty=2, lwd=2)
lowline <- lines( 1/idvals, simchis[1,]+(+2)*simchis[2,], lty=2, lwd=2)
lines(1/idvals, approxchi$a-idvals*approxchi$b, col="red", lwd=2)
axis(side=4, at=zeroinf(mu,Dmat,Gmat)["zero"], labels=expression(delta == 0), las=1)
axis(side=4, at=zeroinf(mu,Dmat,Gmat)["infty"], labels=expression(delta[40] == infinity), las=1)
## eight patches
load("sims-8patch-log-smooth-mal.RData")  # eight patches smooth
lines( 1/idvals, simchis[1,], col="blue", lwd=2 )
lowline <- lines( 1/idvals, simchis[1,]+(-2)*simchis[2,], lty=2, col="blue", lwd=2)
lowline <- lines( 1/idvals, simchis[1,]+(+2)*simchis[2,], lty=2, col="blue", lwd=2)
lines(1/idvals, approxchi$a-idvals*approxchi$b, col="red", lwd=2)
axis(side=4, at=zeroinf(mu,Dmat,Gmat)["infty"], labels=expression(delta[8] == infinity), las=1)
## legend
legend("topright", lty=c(1,1,1,2), col=c("black","blue","red","black"), legend=c("40 patches", "8 patches", "High dispersal", "+/- 2 SE"), cex=.75, lwd=2)
dev.off()
