
  
source("Lyapunov-approx-fns.R")

mu <- c(1,10)
sigma <- c(4,4)
Dmat <- rbind( c(-0.1,0.1), c(1,-1) )
h <- .0001
tmax <- 30  # should be at least 1/sqrt(h)
nsamps <- 200


# Across different dispersal levels:
#   evaulate at dvals*Dmat
dvals <- 1+4*(0:5)
exact <- sapply(dvals, function(d) { twopatch(mu,sigma,d*Dmat)$lyap } )
approx <- lapply(dvals, function (d) { 
        ox <- ourex(mu, d*Dmat, diag(sigma))
        P <- talayhelp(ox$B,ox$G)
        euler <- lyapprox(ox$B,ox$G,h,tmax,dsamp=0.1,method="eulermethod")
        milshtein <- lyapprox(ox$B,ox$G,h,tmax,dsamp=0.1,method="milshteinmethod",P=P)
        talay <- lyapprox(ox$B,ox$G,h,tmax,dsamp=0.1,method="talaymethod",P=P)
        return( list( euler=euler, milshtein=milshtein, talay=talay ) )
    } )
estimates <- sapply(approx, function (x) c(x$euler$lyap, x$euler$lyapsd, x$milshtein$lyap, x$milshtein$lyapsd,x$talay$lyap, x$talay$lyapsd))

save(approx, estimates, mu, sigma, Dmat, h, tmax, nsamps, dvals, exact, file="Lcheck.Rdata")

pdf(file="Lcheck.pdf")
plot(dvals, estimates[1,],ylim=range(c(estimates[1,]-estimates[2,],estimates[1,]+estimates[2,],exact)), col="red")
arrows(dvals,estimates[1,], dvals, estimates[1,]+estimates[2,], length=0, col="red")
arrows(dvals,estimates[1,], dvals, estimates[1,]-estimates[2,], length=0, col="red")
points(dvals, estimates[3,],ylim=range(c(estimates[3,]-estimates[4,],estimates[3,]+estimates[4,],exact)), col="green")
arrows(dvals,estimates[3,], dvals, estimates[3,]+estimates[4,], length=0, col="green")
arrows(dvals,estimates[3,], dvals, estimates[3,]-estimates[4,], length=0, col="green")
points(dvals, estimates[5,],ylim=range(c(estimates[5,]-estimates[6,],estimates[5,]+estimates[6,],exact)), col="purple")
arrows(dvals,estimates[5,], dvals, estimates[5,]+estimates[6,], length=0, col="purple")
arrows(dvals,estimates[5,], dvals, estimates[5,]-estimates[6,], length=0, col="purple")
lines(dvals, exact)
dev.off()

