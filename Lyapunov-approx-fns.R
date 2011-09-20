
require(MASS)    # for ginv
require(Matrix)  # for expm

# X in an n-vector
# mu is a n-vector
# Dmat is an (n x n) array
# G is an (r x n) array
# S = t(G)%*%G is a (n x n) array
# U is an r-vector of random variables
# h is a time increment
# R is sum(X) .

lyapprox <- function( mu, Dmat, G, h, tmax, method="explmethod", noise="rnorm", dsamp=0, X=NULL, nsegs=10, pathest=FALSE, ... ) {
    # Starting at X (randomly chosen if not specified),
    # repeatedly evolve forward in time steps of h
    # using method(X,mu,Dmat,G,h,...) to approximate the time steps
    # renormalizing after each
    #
    # recording state every dsamp time units (zero => every step, NA => not at all)
    #   should make dsamp small but large enough there is some independence between time points.
    #
    # To get the variance down to the same order as the approximation,
    # need total number of samples of order h^(-3/2), e.g. tmax = h^(-1/2).
    #
    # if pathest is TRUE, estimate Y*mu and Y Sigma Y along the path.

    n <- dim(Dmat)[1]; r <- dim(G)[1]
    S <- t(G) %*% G
    diagS <- diag(S)
    # where to begin?
    if (is.null(X)) { X <- runif(n); X <- X/sum(X) }

    if (!is.function(method)) { method <- get(method) }
    if (!is.function(noise)) { noise <- get(noise) }

    # how many times will we step forward?
    nreps <- ceiling( tmax/h );
    # how often will we record the state?
    if (!is.finite(dsamp)) { dnsamp <- NA }  # never
    else if (dsamp==0) { dnsamp <- 1; }   # all the time
    else { dnsamp <- floor( dsamp/h ); }
    # here's where we store the states:
    if (!is.na(dnsamp)) {
        Xhist <- matrix(0, nrow=floor(nreps/dnsamp), ncol=n)
        Xhist[1,] <- X
    } else { Xhist <- X }

    # pre-compute!
    exphDmat <- expm(h*Dmat)


    # also record Y * mu and Y * Sigma * Y as we go along:
    Ymu <- 0; Ymuvals <- rep(0,nsegs)
    ySy <- 0; ySyvals <- rep(0,nsegs)

    # the approximation accumulates here
    lyap <- 0

    # estimate variability by splitting out at these points
    segpts <- floor(nreps/nsegs)
    segvals <- rep(0,nsegs)

    for (k in 1:nreps) {
        U <- noise(r)
        X <- method(X,mu,exphDmat,G,diagS,h,n,r,U,...)
        if (!is.na(dnsamp) & (k %% dnsamp == 0)) {
            Xhist[k/dnsamp,] <- X
        }
        # renormalize X
        R <- sum(X); X <- X/R
        if (pathest) {
            Ymu <- Ymu + X %*% mu
            ySy <- ySy + t(X) %*% S %*% X
        }
        if (R <= 0 | any(X<0) ) { message("Went to (or below) zero:", X); break; }
        lR <- log(R)
        lyap <- lyap + lR
        if (k%%segpts == 0) { 
            segvals[k/segpts] <- lyap 
            Ymuvals[k/segpts] <- Ymu
            ySyvals[k/segpts] <- ySy
        }
    }
    lyap <- lyap/(k*h)

    # variability estimate:
    #   note that our estimate lyap is the mean of segvars/(segpts*h),
    #   so to estimate the SD we divide its SD by sqrt(nsegs)
    segvars <- diff(c(0,segvals))
    lyapsd <- sqrt( var(segvars/(segpts*h)) / nsegs )

    # and means and variances of Y*mu and Y^T S T along the path
    Ymu <- Ymu / nreps
    ySy <- ySy / nreps
    Ymuvals <- diff(c(0,Ymuvals))
    Ymusd <- sqrt( var(Ymuvals/segpts) / nsegs )
    ySyvals <- diff(c(0,ySyvals))
    ySysd <- sqrt( var(ySyvals/segpts) / nsegs )

    return(list(lyap=lyap, lyapsd=lyapsd, Xhist=Xhist, h=h, tmax=h*k, Ymu=Ymu, ySy=ySy, Ymusd=Ymusd, ySysd=ySysd, mu=mu, Dmat=Dmat, G=G))
}


explmethod <- function (X, mu, exphDmat, G, diagS, h, n, r, U) {
    # do a step of the two-step "exponential" method
    # first the growth:
    X <- X * exp( sqrt(h) * apply(G*U,2,sum) + h * (mu - diagS/2) )   # yes, entry-wise exponential
    X <- as.vector( exphDmat %*% X )   # matrix exponential
    # return( X + h*Dmat %*% X )
    return(X)
}

# asymptotic values
zeroinf <- function (mu, Dmat, G) {
    # compute chi(0) and chi(infty)
    # here G is a (r x n) matrix.
    S <- t(G) %*% G
    J <- Dmat
    J[,] <- 1
    pi <- solve( t(Dmat+J), rep(1,dim(J)[1]) )
    chi <- c(NA,NA)
    names(chi) <- c("zero", "infty")
    chi[1] <- max( mu - diag(S)/2 )
    chi[2] <- sum(pi*mu) - t(pi)%*%S%*%pi/2
    return(chi)
}

#####
# Fancier version
# now G is a matrix

# X in an n-vector
# B is an (n x n) array
# G is an (r x n x n) array
# U is an r-vector of random variables
# h is a time increment
# R is sum(X) .

lyapprox.general <- function( B, G, h, tmax, method="eulermethod", noise="rnorm", dsamp=0, X=NULL, nsegs=10, gdriftfn=NULL, ... ) {
    # Starting at X (randomly chosen if not specified),
    # repeatedly evolve forward in time steps of h
    # using method(X,B,G,h,...) to approximate the time steps
    # renormalizing after each
    #
    # recording state every dsamp time units (zero => every step, NA => not at all)
    #
    # To get the variance down to the same order as the approximation,
    # need total number of samples of order h^(-3/2), e.g. tmax = h^(-1/2).
    #
    # if gdriftfn is not NULL then will store and report change-of-measure information
    #   gdriftfn should be a function as returned by e.g. gdrift below

    n <- dim(B)[1]; r <- dim(G)[1]
    # where to begin?
    if (is.null(X)) { X <- runif(n); X <- X/sum(X) }

    if (!is.function(method)) { method <- get(method) }
    if (!is.function(noise)) { noise <- get(noise) }

    # how many times will we step forward?
    nreps <- ceiling( tmax/h );
    # how often will we record the state?
    if (!is.finite(dsamp)) { dnsamp <- NA }  # never
    else if (dsamp==0) { dnsamp <- 1; }   # all the time
    else { dnsamp <- floor( dsamp/h ); }
    # here's where we store the states:
    if (!is.na(dnsamp)) {
        Xhist <- matrix(0, nrow=floor(nreps/dnsamp), ncol=n)
        Xhist[1,] <- X
    } else { Xhist <- X }

    # the approximation
    lyap <- 0
    # girsanov change-of-measure functional
    girfun <- c(0,0)

    # estimate variability by splitting out at these points
    segpts <- floor(nreps/nsegs)
    segvals <- rep(0,nsegs)
    # and record girsanov information here
    girstep <- 0
    girvals <- rep(0,2*nsegs)
    dim(girvals) <- c(nsegs,2)

    for (k in 1:nreps) {
        U <- noise(r)
        X <- method(X,B,G,h,n,r,U,...)
        if (!is.na(dnsamp) & (k %% dnsamp == 0)) {
            Xhist[k/dnsamp,] <- X
        }
        # renormalize X
        R <- sum(X); X <- X/R
        if (R <= 0 | any(X<0) ) { message("Went to (or below) zero:", X); break; }
        lR <- log(R)
        lyap <- lyap + lR
        if (k%%segpts == 0) { segvals[k/segpts] <- lyap }
        # compute Girsanov step
        if (!is.null(gdriftfn)) {
            girstep <- gdriftfn(X)
            girfun <- girfun + c(U%*%girstep, sum(girstep^2))
            if (k%%segpts == 0) { girvals[k/segpts,] <- girfun }
        }
    }
    lyap <- lyap/(k*h)

    # variability estimate:
    #   note that our estimate lyap is the mean of segvars/(segpts*h),
    #   so to estimate the SD we divide its SD by sqrt(nsegs)
    segvars <- diff(c(0,segvals))
    lyapsd <- sqrt( var(segvars/(segpts*h)) / nsegs )

    # Girsanov information
    if (!is.null(gdriftfn)) {
        girsanov <- cbind(segvals, girvals)
        dimnames(girsanov)[[2]] <- c("lyap","dB","dt")
    } else { girsanov <- NULL }

    return(list(lyap=lyap, lyapsd=lyapsd, Xhist=Xhist, h=h, tmax=h*k, girsanov=girsanov))
}

eulermethod <- function (X, B, G, h, n, r, U) {
    # do a step of the "Euler" method
    H <- sqrt(h) * apply(G*U,c(2,3),sum) + h * B
    return( X + H%*%X )
}

milshteinmethod <- function (X, B, G, h, n, r, U, P=talayhelp(B,G)) {
    # do a step of the Milshtein method
    #  for comparison to the others
    # P$GG is such that GG[i,j,k,l] = G[i,k,] %*% G[j,,l]
    Z <- diag(0,r)
    Z[upper.tri(Z)] <- 2*rbinom(r*(r-1)/2, 1, 1/2) - 1
    Z <- Z - t(Z) + diag(1,r)
    S <- rep(outer(U,U) - Z, n^2)
    H <- sqrt(h) * apply(G*U,c(2,3),sum) + h * ( B + (1/2)*apply(P$GG*S,c(3,4),sum) )
    return( X + H%*%X )
}


talaymethod <- function (X, B, G, h, n, r, U, P=talayhelp(B,G)) {
    # do the h^2 method from Talay
    # P$GG is such that GG[i,j,k,l] = G[i,k,] %*% G[j,,l]
    # P$GB is such that GB[i,j,k] = G[i,j,] %*% B[,k]
    # P$BG is such that BG[i,j,k] = B[j,] %*% G[i,,k]
    # P$BB = B%*%B
    Z <- diag(0,r)
    Z[upper.tri(Z)] <- 2*rbinom(r*(r-1)/2, 1, 1/2) - 1
    Z <- Z - t(Z) + diag(1,r)  # Z should be antisymmetric but with 1 on the diagonal
    S <- rep(outer(U,U) - Z, n^2)
    H <- sqrt(h) * apply(G*U,c(2,3),sum) + h * ( B + (1/2)*apply(P$GG*S,c(3,4),sum) ) + h^(3/2) * (1/2) * apply(P$BG + P$GB, c(2,3), sum) + h^2 * (1/2) * P$BB
    return( X + H%*%X )
}

talayhelp <- function (B,G) {
    # help precompute some things for Talay's method:
    # return:
    #   GG is such that GG[i,j,k,l] = G[i,k,] %*% G[j,,l]
    #   GB is such that GB[i,j,k] = G[i,j,] %*% B[,k]
    #   BG is such that BG[i,j,k] = B[j,] %*% G[i,,k]
    #   BB = B%*%B
    r <- dim(G)[1]
    n <- dim(G)[2]
    GG <- array(1,c(r,r,n,n))
    BG <- GB <- array(1, c(r,n,n))
    BB <- B%*%B
    for (i in 1:r) {
        for (j in 1:r) {
            GG[i,j,,] <- G[i,,] %*% G[j,,]
        }
        GB[i,,] <- G[i,,] %*% B
        BG[i,,] <- B %*% G[i,,]
    }
    return( list( GG=GG, GB=GB, BG=BG, BB=BB ) )
}


# plot the results
plotsim <- function( x ) {
    # x is what is returned from lyapprox
    layout(1:2)
    n <- dim(x$Xhist)[2]
    # plot marginal values
    plot(0, xlim=c(1,dim(x$Xhist)[1]), ylim=c(0,max(x$Xhist)))
    sapply(1:n, function(k) { lines(x$Xhist[,k],col=rainbow(n)[k]) } )
    # And plot the log growth rate
    zz <- log(x$Xhist %*% rep(1,n))/h
    wpoints <- (0:99)*length(zz)/100
    window <- floor(length(zz)/50)
    # with a running mean (should converge to estimated value)
    running <- sapply(wpoints, function (k) { mean( zz[k+1:window], na.rm=TRUE) } )
    plot(zz, xlim=c(0,length(zz)), ylim=x$lyap + c(-1,1)*(max(running)-min(running)), pch=".")
    lines(wpoints + length(zz)/200, running, col="red")
    abline(h=x$lyap, col="green", lwd=2)
    abline(h=x$lyap+x$lyapsd*c(-2,2), col="green", lty=2, lwd=2)
    # Also add in Ymu-ySy/2 if this exists
    if (!is.null(x$mu)) {
        S <- t(x$G) %*% x$G
        Ymuvals <- (x$Xhist %*% x$mu )/( x$Xhist %*% rep(1,n) )
        ySyvals <- ( x$Xhist %*% S %*% t(x$Xhist) )/( x$Xhist %*% rep(1,n) )^2
        zzz <- Ymuvals - ySyvals/2
        running <- sapply(wpoints, function (k) { mean( zzz[k+1:window], na.rm=TRUE) } )
        points(zzz, col="blue", pch=".")
    lines(wpoints + length(zzz)/200, running, col="red")
    }
}

# Translate parameters for our SDE
ourex <- function (mu, Dmat, Gmat) {
    # translate from our notation into the more general notation
    # Gmat is our \Gamma, Dmat is our D, and mu is our \mu.
    n <- dim(Gmat)[2]; r <- dim(Gmat)[1]
    B <- diag(mu) + t(Dmat)

    # G is an array whose nonzero entries are Gmat,
    # and those are such that G[i,j,j] = Gmat[i,j].
    G <- rep(Gmat,n)
    x <- rep(row(diag(n)), each=r)
    y <- rep(col(diag(n)), each=r)
    G <- G*( x==y )
    dim(G) <- c(r,n,n)
    return( list(B=B, G=G) )
}

hdlyap <- function (mu,Gmat,Q) {
    # compute high-dispersal approximation:
    # chi(delta) is approximately
    #  (\mu^T \pi - (1/2)\pi^T \Sigma \pi) 
    #   - 1/delta * (  (\mu - \Sigma \pi)^T (Q^T)^- \diag(\pi) (\mu - \Sigma \pi)
    #       + (1/2) \int_0^\infty Tr( \exp(s Q^T) (\diag(\pi)-\pi\pi^T) 
    #           \Sigma (\diag(\pi)-\pi\pi^T) \exp(s Q) \Sigma ) ds
    #
    # Recall \Sigma = \Gamma^T \Gamma
    #  \pi Q = 0 and \sum \pi = 1
    #  and D = Dmat = \delta Q .
    n <- dim(Gmat)[2]; r <- dim(Gmat)[1]
    Sigma <- t(Gmat) %*% Gmat

    # find pi and (Q^T)^-
    pi <- qr.solve( t(cbind(Q,1)), c(rep(0,n),1) )
    QTinv <- ginv(t(Q))
    dpi <- diag(as.vector(pi))

    # value at delta=infty
    lyapinf <- pi %*% (mu - Sigma %*% pi / 2)

    # slope
    lyapslope <- t(mu - Sigma %*% pi) %*% QTinv %*% dpi %*% (mu - Sigma %*% pi)
    # integrand must be vectorized
    MM <- (dpi - pi%*%t(pi)) %*% Sigma %*% (dpi - pi%*%t(pi))
    integrand <- function (sv) { sapply(sv, function (s) {
            expsQ <- expm(s*Q)
            sum( diag( t(expsQ) %*% MM %*% expsQ %*% Sigma ) )
        } ) }
    integral <- integrate(integrand, 0, Inf)
    lyapslope <- lyapslope + integral$value

    return(list(a=lyapinf,b=lyapslope,integral=integral))
}


# Girsanov code (using our notation again)
gdrift <- function (Qmat, Gmat) {
    # function to compute the "drift" term to be used in Girsanov calculations
    Gmatinv <- ginv(Gmat)
    return ( function (X) { Gmatinv %*% diag(as.vector(1/X)) %*% t(Qmat) %*% X } )
}

girlyap <- function ( gir, delta ) {
    # use girsanov information as recorded by lyapprox
    #  to estimate Lyapunov exponent
    #
    # RN derivative is exp( (delta-1) gir$dB[k] + (delta-1)^2 gir$dt[k] )
    #
    # for now using smallest window size
    lyapvals <- diff(gir$lyap)
    RN <- exp( (delta-1) * diff(gir$dB) + (delta-1)^2 * diff(gir$dt) )
    return( mean( lyapvals * RN ) )
}

# Exact results
twopatch <- function (mu, sigma, Dmat) {
    # mu, and sigma are vectors of length 2; Dmat is a 2x2 matrix.
    a <- 2 * sigma^2 / sum(sigma^2)
    b <- (mu[1]-mu[2]+Dmat[2,1]-Dmat[1,2])*2/sum(sigma^2)

    # the density
    frho <- function (y) {
        exp( -(Dmat[2,1]/y + Dmat[1,2]/(1-y))*2/sum(sigma^2) ) * y^(b-a[1]) * (1-y)^(-b-a[2])
    }
    Z <- integrate( frho, 0, 1 )

    f <- function (y) { 
        ((mu[1]-mu[2]+sigma[2]^2)*y - sum(sigma^2)*y^2/2) * #
        exp( -(Dmat[2,1]/y + Dmat[1,2]/(1-y))*2/sum(sigma^2) ) * y^(b-a[1]) * (1-y)^(-b-a[2]) / Z$value
    }
    chi <- mu[2] - sigma[2]^2/2 + integrate( f, 0, 1 )$value

    return(list(lyap=chi,rho=function(y) { frho(y)/Z$value } ))
}


###
# Use package sde

use.sde.sim <- function( mu, Dmat, sigma, h, tmax, X=runif(1), ... ) {
    # Simulate using sde.sim
    # and compute mu*Y and Y * Sigma * Y.

    require(sde)

    G <- diag(sigma)
    S <- t(G) %*% G

    drift.expr <- as.expression( substitute( x * (1-x) * ( mu1 - mu2 - x * s1^2 - (1-x) * s2^2 ) - x * D12 + (1-x) * D21, list( s1=sigma[1], s2=sigma[2], D12=Dmat[1,2], D21=Dmat[2,1], mu1=mu[1], mu2=mu[2] ) ) )
    A.expr <- as.expression( substitute( 
            x * D21 + (1/2)*x^2 * (mu1-mu2-D12-D21-s2^2)
            + (1/3)*x^3 * (-mu1+mu2-s1^2)
            + (1/4)*x^4 * (s1^2-s2^2)
        , list( s1=sigma[1], s2=sigma[2], D12=Dmat[1,2], D21=Dmat[2,1], mu1=mu[1], mu2=mu[2] ) ) )
    A.fn <- function(x) { eval(A.expr) }
    drift.expr.x <- as.expression( substitute( 
            (1-x) * ( mu1 - mu2 - x * s1^2 - (1-x) * s2^2 ) 
            - x * ( mu1 - mu2 - x * s1^2 - (1-x) * s2^2 ) 
            + x * (1-x) * ( s1^2 - s2^2 )
            - D12 - D21
        , list( s1=sigma[1], s2=sigma[2], D12=Dmat[1,2], D21=Dmat[2,1], mu1=mu[1], mu2=mu[2] ) ) )
    drift.expr.xx <- as.expression( substitute( 
            - ( mu1 - mu2 - x * s1^2 - (1-x) * s2^2 ) 
            + (1-x) * ( - s1^2 + s2^2 )
            - ( mu1 - mu2 - x * s1^2 - (1-x) * s2^2 ) 
            - x * ( - s1^2 + s2^2 ) 
            + (1-x) * ( s1^2 - s2^2 )
            - x * ( s1^2 - s2^2 )
        , list( s1=sigma[1], s2=sigma[2], D12=Dmat[1,2], D21=Dmat[2,1], mu1=mu[1], mu2=mu[2] ) ) )
    sigma.expr <- as.expression( substitute( x^2 * (1-x)^2 * (s1^2 + s2^2), list( s1=sigma[1], s2=sigma[2], D12=Dmat[1,2], D21=Dmat[2,1] ) ) )
    sigma.expr.x <- as.expression( substitute( 
            2*x * (1-x)^2 * (s1^2 + s2^2)
            - x^2 * 2*(1-x) * (s1^2 + s2^2)
        , list( s1=sigma[1], s2=sigma[2], D12=Dmat[1,2], D21=Dmat[2,1] ) ) )
    sigma.expr.xx <- as.expression( substitute( 
            2 * (1-x)^2 * (s1^2 + s2^2)
            - 2 * x * 2*(1-x) * (s1^2 + s2^2)
            - 2 * x * 2*(1-x) * (s1^2 + s2^2)
            + x^2 * 2 * (s1^2 + s2^2)
        , list( s1=sigma[1], s2=sigma[2], D12=Dmat[1,2], D21=Dmat[2,1] ) ) )
    # psi is (1/2)*( drift(x)^2 + drift.x(x) )
    # k1 and k2 are lower and upper bounds on it.
    k2 <- 10 * ( sum(mu) + sum(sigma^2) + Dmat[1,2] + Dmat[2,1] )^2

    N <- tmax / h

    Xhist <- sde.sim(t0=0, N=N, T=tmax, X0=runif(1), drift=drift.expr, sigma=sigma.expr, drift.x=drift.expr.x, drift.xx=drift.expr.xx, sigma.x=sigma.expr.x, sigma.xx=sigma.expr.xx, drift.t=expression(0), A=A.fn, k2=k2 ) 

    # and means and variances of Y*mu and Y^T S T along the path
    Y <- cbind(Xhist, 1-Xhist)
    Ymuvals <- apply(Y, 1, function (y) { y %*% mu } )
    ySyvals <- apply(Y, 1, function (y) { y%*%S%*%y } )
    Ymu <- mean(Ymuvals)
    ySy <- mean(ySy)
    Ymusd <- sqrt( var(Ymuvals) )
    ySysd <- sqrt( var(ySyvals) )

    lyap <- Ymu - ySy/2
    lyapsd <- Ymusd + ySysd

    return(list(lyap=lyap, lyapsd=lyapsd, Xhist=Xhist, h=h, tmax=h*N, Ymu=Ymu, ySy=ySy, Ymusd=Ymusd, ySysd=ySysd, mu=mu, Dmat=Dmat, G=G))

}


