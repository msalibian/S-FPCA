# read functions
source('S-FPCA-functions.R')


# www.mortality.org
# Human Mortality Database. University of California, Berkeley (USA), 
# and Max Planck Institute for Demographic Research (Germany). Available 
# at www.mortality.org or www.humanmortality.de (data downloaded on 8 Feb 2013). 

# (FRANCE)
# read male mortality rates, by age group and year 
a <- read.table('FRA_Death_rates_1x1.txt', header=TRUE)

yrs <- unique(a$Year)

# only use data until 1948
# pattern changes in later years...
yrs <- yrs[ yrs < 1949 ] 

# find which death-rates are available in each year
# (there were missing entries, data was cleaned by hand)
ags <- rep(0, length(yrs))
for(i in 1:length(yrs))
    ags[i] <- length(unique(a$Age[ a$Year == yrs[i] ]))

# there are at least 103 age groups with data for each year
# 103 == min(ags)

ma <- min(ags)

# set up the data matrix
mo <- matrix(0, length(yrs), ma)
for(i in 1:length(yrs))
    mo[i, ] <- (a$Male[ a$Year == yrs[i] ])[1:ma]

# take log
mo <- log(mo[, 1:100])

# plot the data
matplot(x=0:99, y=t(mo), xlab='Age', ylab='log(mortality)', 
  lty=1, type='l', col='gray', lwd=2)

# Prepare mesh for spline basis
# move it away from the boundary of [0,1]
mesh <- (1:100)/110

# assign the data to the object "x"
x <- mo

# Analysis
n <- dim(x)[1]
dimension.Bspline <- 20
aa <- generar.mesh.splines(dimension.Bspline )
base.estim.Bspline <- devolver.base01.sieves.ortonormal(mesh, dimension.Bspline)

y <- x %*% base.estim.Bspline

mui <- l1median(X=y,trace=-1) # pcaPP::
q <- 2
Ncand <- 50 # 1000

mu.hat <- as.vector( base.estim.Bspline %*% mui )

# plot the data and the initial location function
matplot(x=mesh, y=t(x), lty=1, lwd=2, col='gray', type='l') #, ylim=c(1, 7))
lines(x=mesh, y=mu.hat, lwd=2, col='black')

# S-estimators
y.sfpca <- sfpca(x=y, mu=mui, q=q, Ncand=Ncand, seed=123, init.it=50, max.it=500, tol=1e-10, 
trace=TRUE, tuning.rho=3, bb=0.2426) 


# build the predictions

mu.hat <- as.vector( base.estim.Bspline %*% y.sfpca$mu )
#phi.hat <- base.estim.Bspline %*% y.sfpca$b
#x.hat.s <- base.estim.Bspline %*% t( y.sfpca$x.s ) 
x.hat.ls <- base.estim.Bspline %*% t( y.sfpca$x.ls ) 

yc <- scale(y, center=y.sfpca$mu, scale=FALSE)
bon <- qr.Q(qr(y.sfpca$b))
y.rulo.s <- scale( (yc %*% bon) %*% t( bon ), center=-y.sfpca$mu, scale=FALSE)
x.hat.rulo.s <- base.estim.Bspline %*% t( y.rulo.s ) 

# show the predicted curves
matplot(x=mesh*110, y=x.hat.rulo.s, 
lty=1, type='l', lwd=4, col='gray',xlab="Age", ylab="Robust Predictors", cex.axis=1.5, cex.lab=1.5, ylim=c(-10, 0)) 
lines(x=mesh*110, y=mu.hat, lwd=4, col='black')


## differences between S & LS predictions

# S- and LS-prediction residuals
re.s <- colMeans((t(x)-x.hat.rulo.s)^2)
re.ls <- colMeans((t(x)-x.hat.ls)^2)

plot(yrs, re.s, type='b', pch=19, lwd=3, cex=2.5, col='gray70',ylab="Residual squared norm", cex.axis=1.5, cex.lab=1.5)
lines(yrs, re.ls, type='b', pch=19, lwd=3, cex=2.5, col='black')

# focus on the LS prediction residuals 
plot(yrs, re.ls, type='b', pch=19, lwd=3, cex=2.5, col='black',ylab="Residual squared norm", cex.axis=1.5, cex.lab=1.5)
lines(yrs, re.s, type='b', pch=19, lwd=3, cex=2.5, col='gray70')


# cut-off for outliers
plot(yrs, re.ls, type='b', pch=19, lwd=3, cex=2.5, col='black', ylab="Residual squared norm", cex.axis=1.5, cex.lab=1.5)
lines(yrs, re.s, type='b', pch=19, lwd=3, cex=2.5, col='gray70')
abline(h=0.035, lwd=4, col='gray50', lty=2)


n <- length(re.s)
ous <- (1:n)[re.s > 0.035] 
ous.ls <- (1:n)[re.ls > 0.035]
mesh110 <- mesh*110

#> yrs[ous]
# [1] 1855 1871 1914 1915 1916 1917 1918 1919 1940 1942 1943 1944 1945 1946 1947
#[16] 1948 
#> yrs[ous.ls]
#[1] 1914 1915 1943 1944 1945 1946 1947 1948

 
# identify the curves labelled as outliers
# Event in these years
yrs[ous]
# 1855      Crimean War
# 1871      Franco-Prussian War
# 1914-18   WWI
# 1918-19   Spanish Flu
# 1940-1945 WWII  
# 1946-     post-war change


# wars, epidemics!
matplot(x=mesh*110, y=t(x), lty=1, lwd=.7, col='gray', type='l', ylab="X(t)",xlab="Age") # ylim=c(1, 7))
lines(x=mesh*110, y=mu.hat, lwd=2, col='black')
for(i in ous) lines(x[i,] ~ mesh110, lwd=4, col='gray40')

# crimean war - S
matplot(x=mesh*110, y=t(x), lty=1, lwd=.7, col='gray', type='l', ylab="X(t)",xlab="Age", 
cex.axis=2, cex.lab=2) # ylim=c(1, 7))
i <- ous[1]; lines(x[i,] ~ mesh110, lwd=4, col='gray40', lty=1)
lines(x.hat.rulo.s[,i] ~ mesh110, lwd=4, col='black', lty=2)

# crimean war - LS
matplot(x=mesh*110, y=t(x), lty=1, lwd=.7, col='gray', type='l', ylab="X(t)",xlab="Age", 
cex.axis=2, cex.lab=2) # ylim=c(1, 7))
i <- ous[1]; lines(x[i,] ~ mesh110, lwd=4, col='gray40', lty=1)
lines(x.hat.ls[,i] ~ mesh110, lwd=4, col='black', lty=2)

# prusian war - S
matplot(x=mesh*110, y=t(x), lty=1, lwd=.7, col='gray', type='l', ylab="X(t)",xlab="Age", 
cex.axis=2, cex.lab=2) # ylim=c(1, 7))
i <- ous[2]; lines(x[i,] ~ mesh110, lwd=4, col='gray40', lty=1)
lines(x.hat.rulo.s[,i] ~ mesh110, lwd=4, col='black', lty=2)

# prusian war - LS
matplot(x=mesh*110, y=t(x), lty=1, lwd=.7, col='gray', type='l', ylab="X(t)",xlab="Age", 
cex.axis=2, cex.lab=2) # ylim=c(1, 7))
i <- ous[2]; lines(x[i,] ~ mesh110, lwd=4, col='gray40', lty=1)
lines(x.hat.ls[,i] ~ mesh110, lwd=4, col='black', lty=2)


# WWI - Spanish flu - S
matplot(x=mesh*110, y=t(x), lty=1, lwd=.7, col='gray', type='l', ylab="X(t)",xlab="Age", 
cex.axis=2, cex.lab=2) # ylim=c(1, 7))
for(i in ous[3:8]) lines(x[i,] ~ mesh110, lwd=4, col='gray40')
for(i in ous[3:8]) lines(x.hat.rulo.s[,i] ~ mesh110, lwd=4, col='black', lty=2)

# WWI - Spanish flu - LS
matplot(x=mesh*110, y=t(x), lty=1, lwd=.7, col='gray', type='l', ylab="X(t)",xlab="Age", 
cex.axis=2, cex.lab=2) # ylim=c(1, 7))
for(i in ous[3:8]) lines(x[i,] ~ mesh110, lwd=4, col='gray40')
for(i in ous[3:8]) lines(x.hat.ls[,i] ~ mesh110, lwd=4, col='black', lty=2)

# WWII - S
matplot(x=mesh*110, y=t(x), lty=1, lwd=.7, col='gray', type='l', ylab="X(t)",xlab="Age", 
cex.axis=2, cex.lab=2) # ylim=c(1, 7))
for(i in ous[9:13]) lines(x[i,] ~ mesh110, lwd=4, col='gray40')
for(i in ous[9:13]) lines(x.hat.rulo.s[,i] ~ mesh110, lwd=4, col='black', lty=2)

# WWII - LS
matplot(x=mesh*110, y=t(x), lty=1, lwd=.7, col='gray', type='l', ylab="X(t)",xlab="Age", 
cex.axis=2, cex.lab=2) # ylim=c(1, 7))
for(i in ous[9:13]) lines(x[i,] ~ mesh110, lwd=4, col='gray40')
for(i in ous[9:13]) lines(x.hat.ls[,i] ~ mesh110, lwd=4, col='black', lty=2)



# outliers identified by LS
matplot(x=mesh, y=t(x), lty=1, lwd=.7, col='gray', type='l') # ylim=c(1, 7))
for(i in ous.ls) lines(x[i,] ~ mesh, lwd=3, col='darkblue')


# outliers identified by S
matplot(x=mesh, y=t(x), lty=1, lwd=.7, col='gray', type='l') # ylim=c(1, 7))
for(i in ous) lines(x[i,] ~ mesh, lwd=3, col='darkblue')


# Outliers found by S but missed by LS
ous.d <- setdiff(ous, ous.ls) 
matplot(x=mesh, y=t(x), lty=1, lwd=.7, col='gray', type='l') # ylim=c(1, 7))
for(i in ous.d) lines(x[i,] ~ mesh, lwd=3, col='springgreen4')









