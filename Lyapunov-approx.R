source("Lyapunov-approx-fns.R")

mu <- c(3,3)
sigma <- c(1,1)
Dmat <- cbind( c(-1,1), c(1,-1) )
h <- .001
tmax <- 1000  # should be at least 1/sqrt(h)
nsamps <- 200

# At these parameters:
exact <- twopatch(mu,sigma,Dmat)
yy <- (1:99)/100; plot(yy, exact$rho(yy), type="l")

ox <- ourex(mu, Dmat, diag(sigma))
eapprox <- lyapprox(ox$B,ox$G,h,tmax,dsamp=0,method="eulermethod")
happrox <- lyapprox(ox$B,ox$G,h,tmax,dsamp=0,method="milshteinmethod",P=talayhelp(ox$B,ox$G))
tapprox <- lyapprox(ox$B,ox$G,h,tmax,dsamp=0,method="talaymethod",P=talayhelp(ox$B,ox$G))
ezz <- log( apply( eapprox$Xhist, 1, sum ) )/eapprox$h
hzz <- log( apply( happrox$Xhist, 1, sum ) )/happrox$h
tzz <- log( apply( tapprox$Xhist, 1, sum ) )/tapprox$h

eacf <- acf(ezz)
hacf <- acf(hzz)
tacf <- acf(tzz)

# not so useful?
elow <- lowess(ezz, f=1/32)
hlow <- lowess(ezz, f=1/32)
tlow <- lowess(tzz, f=1/32)
plot(elow, type="l"); lines(hlow, col="green"); lines(tlow, col="red")
abline(h=exact$lyap)

# the approximation as time goes on
cezz <- cumsum(ezz); cezz <- cezz/1:length(ezz)
chzz <- cumsum(hzz); chzz <- chzz/1:length(hzz)
ctzz <- cumsum(tzz); ctzz <- ctzz/1:length(tzz)
plot(cezz,type="l", ylim=range(cezz[-(1:100)]))
abline(h=eapprox$lyap+c(-1,1)*eapprox$lyapsd, lty=2)
abline(h=happrox$lyap+c(-1,1)*happrox$lyapsd, lty=2, col="green")
abline(h=tapprox$lyap+c(-1,1)*tapprox$lyapsd, lty=2, col="red")
lines(chzz,col="green")
lines(ctzz,col="red")
abline(h=exact$lyap)

save(eapprox,happrox,tapprox,ezz,hzz,tzz,file="Lyapprox-run-1.Rdata")

# compare independent time intervals?
nint <- 20
dk <- floor(length(ezz)/nint)
ezzsect <- sapply(1:nint, function (k) { cumsum(ezz[(k-1)*dk + 1:dk])/(1:dk) })
hzzsect <- sapply(1:nint, function (k) { cumsum(hzz[(k-1)*dk + 1:dk])/(1:dk) })
tzzsect <- sapply(1:nint, function (k) { cumsum(tzz[(k-1)*dk + 1:dk])/(1:dk) })
plot(0,xlim=c(1,dk),ylim=range(cezz[-(1:100)]))
apply(ezzsect,2,lines)
apply(hzzsect,2,lines,col="green")
apply(tzzsect,2,lines,col="red")
abline(h=exact$lyap,col="blue",lwd=2)
abline(h=exact$lyap+c(-2,2)*eapprox$lyapsd,lty=2,col="red",lwd=2)

lapply( list(eapprox, happrox, tapprox), function (approx) {
        zz <- log( apply( approx$Xhist, 1, sum ) )/approx$h
        plot(0, ylim=c(0,2), xlim=c(0,length(zz)))
        abline(h=mean(zz))
        for (nn in 2^(1:7)) {
            nwind <- 20
            windpoints <- (0:(nwind-1))*length(zz)/nwind
            window <- floor(length(zz)/nn)
            running <- sapply(windpoints, function (k) {
                    c(mean( zz[k + 1:window], na.rm=TRUE ), median( zz[k + 1:window], na.rm=TRUE ) )
                } )
            lines(windpoints+length(zz)/(2*nwind), running[1,], col="red")
            lines(windpoints+length(zz)/(2*nwind), running[2,], col="green")
        }
    } )


