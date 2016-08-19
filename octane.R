
## Octane dataset
## Esbensen, Sch\"onkopf, and Midtgaard (1994)
##      Multivariate Analysis in Practice, Trondheim: Camo
## n = 39 near-infrared (NIR) absorbance spectra over p = 226
## wavelengths
## samples 25, 26, 36:39 contain added alcohol

# Combinations that work
# octane-q2.R

source('S-FPCA-functions.R')
library(mvoutlier)
library(doParallel)

# Start cluster
cl <- makeCluster(4)
registerDoParallel(cl)

x <- read.csv('octane.csv', header=TRUE, dec=',')
x <- x[ ,-(1:2)]

mui <- l1median(X=x, trace=-1) # pcaPP::
q <- 6
# number of random starts for the iterative algorithm
Ncand <- 1000
x.s <- sfpca.par(x=x, mu=mui, q=q, Ncand=Ncand, seed=123, init.it=50, max.it=500, 
                 tol=1e-6, trace=FALSE, tuning.rho=3, bb = 0.2426) 
# cc=1.54764, b=.5
mu.hat <- x.s$mu
xc <- scale(x, center=x.s$mu, scale=FALSE)
bon <- qr.Q(qr(x.s$b))
x.hat.s <- scale( (xc %*% bon) %*% t( bon ), center=-x.s$mu, scale=FALSE)
re.s <- rowMeans((x-x.hat.s)^2)

#known outliers
outs.true <- c(25, 26, 36:39)

names(re.s) <- 1:nrow(x)
ouS <- as.numeric(names(adjbox(re.s, plot=FALSE)$out))
inliers <- re.s[ouS] < median(re.s)
(outs.s <- ouS[!inliers])

