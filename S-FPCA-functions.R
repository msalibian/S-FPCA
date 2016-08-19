

library(MASS)
library(mgcv)
library(pcaPP)
library(robustbase)


# Tukey's Rho function
# Rho <- function(r, tuning.rho=1.54764) tukeyChi(x=r, cc=tuning.rho, deriv=0)
Rho <- function(r, tuning.rho=1.54764) Mchi(x=r, psi='bisquare', cc=tuning.rho, deriv= 0)

# Tukey's Rho' function
# Rhop <- function(r, tuning.rho=1.54764) tukeyChi(x=r, cc=tuning.rho, deriv=1)
Rhop <- function(r, tuning.rho=1.54764) Mchi(x=r, psi='bisquare', cc=tuning.rho, deriv= 1)


# Given a vector of residuals "r", solve (for "s") the
# equation: mean( Rho(r/s) ) = b
s.scale <- function(r, tuning.rho=1.54764, b=.5, max.it=1000, ep=1e-4) 
{
  s1 <- mad(r)
  if(abs(s1)<1e-10) return(s1)
  s0 <- s1 + 1
  it <- 0
  while( ( abs(s0-s1) > ep ) && (it < max.it) ) {
    it <- it + 1
    s0 <- s1
    s1 <- s0*mean(Rho(r/s0,tuning.rho=tuning.rho))/b
  }
  return(s1)
}


# Compute the S-estimator for PCA
# x = n x p matrix (each row is an observation)
# mu = n-vector with the initial mean 
# q = dimension of the subspace to be estimated
# Ncand = number of random starting points for the iterative algorithm
# seed = seed
# init.it = number of iterations is initial candidate is tried
# max.it = number of iterations applied to the best of the initial candidates 
# tol = convergence criterion (stop when consecutive changes are less than "tol")
# trace = detailed info on interations?
# tuning.rho = tuning constant for Tukey's "rho" function
# bb = corresponding right-hand side constant. To obtain a consistent
#      scale estimator one needs to set bb = E( Rho(Z) ) where Z ~ N(0,1)
# 
# function returns:
# 
# mu       = final vector of means
# sig      = vector of \hat{\sigma}_j, 1 <= j <= p
# best.obj = best objective function = sum( \hat{\sigma}_j^2 )
# x.ls     = LS PCA solution = colMeans(x) + a.ls %*% b.ls
# a.ls     = n x q = matrix of LS loadings = i-th row contains the "loadings" 
#                of x_i on the q-dim basis in b.ls
# b.ls     = p x q = basis of the PCA subspace (each column is a vector)
# x.ls, a, b = same as above, but for the S estimator 
#                  (using mu instead of colMeans(x))
# 
# error.flag = did the last set of iterations end prematurely?

sfpca <- function(x, mu, q, Ncand=50, seed=191, init.it=20, max.it=500, tol=1e-6, 
                  trace=TRUE, tuning.rho=3, bb=0.24) {
  this.call <- match.call()
  # LS solution
  sih <- var(x)
  b.ls <- svd(sih)$u[,1:q]
  mu.ls <- colMeans(x)
  a.ls <- scale(x, center=mu.ls, scale=FALSE) %*% b.ls
  x.ls <- a.ls %*% t(b.ls)
  x.ls <- scale(x.ls, center=-mu.ls, scale=FALSE)
  best.obj <- +Inf
  #mu.med <- apply(x, 2, median)
  #ache <- 0.1
  #mu<- ksmooth(te, mu.med, kernel = "normal", bandwidth = ache,x.points=te)$y
  xc <- scale(x, center=mu, scale=FALSE)
  for(j in 1:Ncand) {
    tmp <- coo.fit(mu=mu, x=x, xc=xc, q=q, seed=seed+17*j, max.it=init.it, tol=tol, 
                   trace=trace, tuning.rho=tuning.rho, bb=bb)
    if(!tmp$error.flag) {
      ob <- sum(tmp$sig^2)
      if(trace) print( ob )
      if( ob < best.obj ) {
        best.obj <- ob
        best.j <- j
      }
    } else {print(paste('error, cand ', j, sep='')) }
  }
  if(trace) print(paste('Best: ', best.obj, sep=''))
  best.tmp <- coo.fit(mu=mu, x=x, xc=xc, q=q, seed=seed+17*best.j, max.it=max.it, tol=tol, 
                      trace=trace, tuning.rho=tuning.rho, bb=bb)
  tmp <- c(best.tmp, best.obj=sum(best.tmp$sig^2), list(a.ls=a.ls, b.ls=b.ls, x.ls=x.ls, call=this.call))
  return(tmp)
}


sfpca.par <- function(x, mu, q, Ncand=50, seed=191, init.it=20, max.it=500, tol=1e-6, 
                      trace=TRUE, tuning.rho=3, bb=0.24) {
  # doParallel version of sfpca()
  this.call <- match.call()
  # LS solution
  sih <- var(x)
  b.ls <- svd(sih)$u[,1:q]
  mu.ls <- colMeans(x)
  a.ls <- scale(x, center=mu.ls, scale=FALSE) %*% b.ls
  x.ls <- a.ls %*% t(b.ls)
  x.ls <- scale(x.ls, center=-mu.ls, scale=FALSE)
  best.obj <- +Inf
  #mu.med <- apply(x, 2, median)
  #ache <- 0.1
  #mu<- ksmooth(te, mu.med, kernel = "normal", bandwidth = ache,x.points=te)$y
  xc <- scale(x, center=mu, scale=FALSE)
  tmp.par <- foreach(j=1:Ncand, .combine=rbind, 
              .inorder=FALSE, .packages='robustbase', 
              .export=c('coo.fit', 's.scale', 'Rho', 'Rhop')) %dopar% {
        tmp <- coo.fit(mu=mu, x=x, xc=xc, q=q, seed=seed+17*j, max.it=init.it, tol=tol, 
                  trace=trace, tuning.rho=tuning.rho, bb=bb)
        if(!tmp$error.flag) {
           ( ob <- sum(tmp$sig^2) )
         } else { NA } 
  }
  best.j <- which.min(tmp.par[,1])
  # if(trace) print(paste('Best: ', best.obj, sep=''))
  best.tmp <- coo.fit(mu=mu, x=x, xc=xc, q=q, seed=seed+17*best.j, max.it=max.it, tol=tol, 
                      trace=trace, tuning.rho=tuning.rho, bb=bb)
  tmp <- c(best.tmp, best.obj=sum(best.tmp$sig^2), list(a.ls=a.ls, b.ls=b.ls, x.ls=x.ls, call=this.call))
  return(tmp)
}



coo.fit <- function(mu, x, xc, q, seed=191, max.it=1000, tol=1e-6, 
                    trace=TRUE, tuning.rho=3, bb=0.24) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  # initial estimators
  set.seed(seed)
  di <- 10*tol
  it <- 0
  # initial
  #mu <- apply(x, 2, median)
  #xc <- scale(x, center=mu, scale=FALSE)
  #b <- matrix(rnorm(p*q), p, q)
  #b <- qr.Q(qr(b))
  #a <- xc %*% b
  ## mu <- rep(0, p)
  #for(i in 1:n) {
  #a[i,] <- solve( t(b) %*% b, t(b) %*% xc[i,])
  #}
  # subsampling for initial
  xc <- scale(x, center=mu, scale=FALSE)
  ii <- sample(n, q)
  b <- t(x[ii,,drop=FALSE])
  b <- qr.Q(qr(b)) # gramm-schimdt de b
  a <- xc %*% b
  r0 <- r <- x - (x.s <- scale(a %*% t(b), center=-mu, scale=FALSE))
  rr <- x - a %*% t(b)
  ne <- sum( (sig <- apply(r, 2, s.scale, tuning.rho=tuning.rho, b=bb))^2 )
  r.sd <- scale(r, center=FALSE, scale=sig)
  ol <- ne + 100*tol
  error.flag <- FALSE
  while( ( (it <- it+1) < max.it ) & ( max(abs(ol-ne)/ol, abs(ol-ne)) > tol) ) {
    if(trace) print(c(ne, ol, ol-ne))
    ol <- ne
    ome <- Rhop(r.sd, tuning.rho=tuning.rho)
    hs <- colSums(ome*r.sd)/sig
    vs <- scale(ome/r, center=FALSE, scale=hs)
    ws <- Rhop(r.sd, tuning.rho=tuning.rho)/r
    if(any(is.na(vs))) {
      vs[is.na(vs)] <- 6/tuning.rho^2
      ws[is.na(ws)] <- 6/tuning.rho^2
    }
    xc <- scale(x, center=mu, scale=FALSE)
    tt <- try( {
      for(i in 1:n) {
        v <- as.vector( vs[i,] )
        a[i,] <- solve( t(b) %*% diag(v) %*% b, t(b) %*% diag(v) %*% xc[i,])
      }
      for(j in 1:p) {
        w <- ws[,j]
        b[j,] <- solve( t(a) %*% diag(w) %*% a, t(a) %*% diag(w) %*% xc[,j])
        mu[j] <- weighted.mean(rr[,j], w=w)
      }
    }, silent=TRUE)
    if(class(tt)=='try-error') it <- max.it
    r0 <- r <- x - (x.s <- scale(a %*% t(b), center=-mu, scale=FALSE))
    rr <- x - a %*% t(b)
    ne <- sum( (sig <- apply(r, 2, s.scale, tuning.rho=tuning.rho, b=bb))^2 )
    rho.p <- Rhop(r.sd <- scale(r, center=FALSE, scale=sig), tuning.rho=tuning.rho) 
  }
  if(class(tt)=='try-error') error.flag <- TRUE 
  r.sd <- scale(r, center=FALSE, scale=sig)
  rho.p <- Rhop(r.sd, tuning.rho=tuning.rho) 
  ome <- Rhop(r.sd, tuning.rho=tuning.rho)
  hs <- colSums(ome*r.sd)/sig
  eq1 <- t(a) %*% rho.p / n
  eq2 <- scale(rho.p, center=FALSE, scale=hs) %*% b / n
  tmp <- list(a=a, b=b, mu=mu, x.s=x.s, eq1=eq1, eq2=eq2, sig=sig,
              tuning.rho=tuning.rho, error.flag=error.flag)
  return(tmp)
}




generar.mesh.splines <- function(k)
{
  if (k == 5)
    part <- 2
  else
    part <- k/2
  res <- ((-part):part)/part
  return(res)
}




devolver.base01.sieves <- function(mesh,k)
{    
  aa <- generar.mesh.splines(k)
  knots <- (aa+1)*1.0/2
  base.estim.disc <- cSplineDes(mesh,knots)
  return(base.estim.disc)
}



devolver.base01.sieves.ortonormal <- function(mesh,k)
{
  base.estim.disc <- devolver.base01.sieves(mesh,k)
  base.ortonorm.disc <- qr.Q(qr(base.estim.disc))
  return(base.ortonorm.disc)
}



