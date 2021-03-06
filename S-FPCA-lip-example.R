

# http://www.stats.ox.ac.uk/~silverma/fdacasebook/lipemg.html

#Applied Functional Data Analysis: Methods and Case Studies
#by James O. Ramsay and Bernard W. Silverman
#Chapter 10: Predicting Lip Acceleration from Electromyography
#    Data on the position of the lower lip, through time for each of 32 records (plain text)
#    Data on the acceleration of the lower lip, through time for each of 32 records (plain text)
#    Data on the EMG activity in the depressor labii inferior, through time for each of 32 records (plain text)
#    Description of the analysis of the lip acceleration and EMG data in MATLAB (pdf)
#    Additional Matlab functions required for the analysis of the data (zip)
#    Background paper by Malfait, Ramsay and Froda (2001) (pdf)
#    Background paper by Malfait, Ramsay and Froda (2001) (postscript)


# read functions
source('S-FPCA-functions.R')

# read data
lippo <- read.table('lip-data.txt')
norm <- function(a) sqrt(sum(a^2))

# prepare data
x <- t(lippo)
te <- seq(0,0.69,length=501)

# plot data
matplot(x=te, y=lippo, ylab='', xlab='Time', type='l', col='gray', lty=1, lwd=2)


# functions to prepare a spline basis in (0, 0.69) (instead of (0,1))
devolver.base069.sieves <- function(mesh,k)
{    
  aa <- generar.mesh.splines(k)
  knots<-(aa+1)*0.69/2
  base.estim.disc <- cSplineDes(mesh,knots)
  return(base.estim.disc)
}

devolver.base069.sieves.ortonormal <- function(mesh,k)
{
  base.estim.disc <-devolver.base069.sieves(mesh,k)
  base.ortonorm.disc<-qr.Q(qr(base.estim.disc))
  return(base.ortonorm.disc)
}

# prepare the spline basis (with 20 functions)
mesh <- te
n <- nrow(x)
dimension.Bspline <- 20 

aa <- generar.mesh.splines(dimension.Bspline )
knots <- (aa+1)*0.69/2
base1 <- devolver.base069.sieves(mesh, dimension.Bspline)
base.estim.Bspline <- devolver.base069.sieves.ortonormal(mesh, dimension.Bspline)

# find the representation of the data on this basis
y <- x %*% base.estim.Bspline

# initial "mean" function (L1-estimate)
mui <- l1median(X=y,trace=-1) # pcaPP::
# dimension of the "best" subspace to be estimated 
q <- 5
# number of random starts for the iterative algorithm
Ncand <- 1000 # 1000

y.sfpca <- sfpca(x=y, mu=mui, q=q, Ncand=Ncand, seed=123, init.it=50, max.it=500, 
                 tol=1e-10, trace=TRUE, tuning.rho=3, bb = 0.2426) 


mu.hat <- as.vector( base.estim.Bspline %*% y.sfpca$mu )
x.hat.ls <- base.estim.Bspline %*% t( y.sfpca$x.ls ) 

yc <- scale(y, center=y.sfpca$mu, scale=FALSE)
bon <- qr.Q(qr(y.sfpca$b))
y.rulo.s <- scale( (yc %*% bon) %*% t( bon ), center=-y.sfpca$mu, scale=FALSE)
x.hat.rulo.s <- base.estim.Bspline %*% t( y.rulo.s ) 

# robust predictions
matplot(x=te, y=x.hat.rulo.s, lty=1, type='l', col='gray', lwd=2, xlab='Time', ylab='')
lines(mu.hat ~ te, lwd=4, col='black')


# S- and LS-prediction residuals
re.s <- colMeans((t(x)-x.hat.rulo.s)^2)
re.ls <- colMeans((t(x)-x.hat.ls)^2)

re.s <- re.s * 1e6
re.ls <- re.ls * 1e6

plot(re.s, type='b', pch=19, lwd=3, cex=2.5, col='gray70',ylab="Residual squared norm")
lines(re.ls, type='b', pch=19, lwd=3, cex=2.5, col='black')
abline(h=0.65, lwd=4, lty=2, col='gray30')

ous <- (1:n)[re.s > .65]

matplot(x=te, y=t(x), lty=1, type='l', col='gray', xlab='Time', ylab='', lwd=2)
for(i in ous) lines(x=te, y=x[i,], lwd=4, col='gray30')






