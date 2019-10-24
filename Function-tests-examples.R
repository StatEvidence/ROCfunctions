
###############################################################
## Set of R commands for testing/illustration ROC functions
## AUC.R, ROC.R, ROCplot.R
## JD Blume
###############################################################

#######################################
## Data Generation
## y.0 ~ N[mu.0,s.0^2] and y.1 ~ N[mu.1, s.1^2]
## equivalent to (under some monotonic transformation)
## y.0 ~ N[0, 1] and y.1 ~ N[a/b,1/b^2]
## where
## a = (mu.1-mu.0)/s.1 and b = s.0/s.1
## because (y.1-mu.0)/s.0 ~ N[ (mu.1-mu.0)/s.0, s.1^2/s.0^2] = N[a/b, 1/b^2]
## ROC curve is sens = pnorm(a+b*qnorm(1-spec))
## AUC is pnorm(a/sqrt(1+b^2)) 
#######################################

##### Test AUC computation
a <- 2		# a = (mu.1-mu.0)/s.1
b <- 2		# b = s.0/s.1
auc.true <- pnorm(a/sqrt(1+b^2)) 	

n <- 400
p.0 <- 0.5

n.0 <- ceiling(p.0*n)
n.1 <- n-n.0

x.0 <- rnorm(n.0, mean=0, sd=1)
x.1 <- rnorm(n.1, mean=a/b, sd=1/b)

bins=quantile(c(x.0,x.1),seq(0,1,0.2))
d.0 <- as.numeric(cut(x.0, breaks=bins, include.lowest=TRUE)) # Avoids NA for min
d.1 <- as.numeric(cut(x.1, breaks=bins, include.lowest=TRUE)) # Avoids NA for min

c(AUC(x.0,x.1,as.groups=TRUE)$auc, AUC(d.0,d.1,as.groups=TRUE)$auc, auc.true)
AUC(x.0,x.1,as.groups=TRUE)
AUC(d.0,d.1,as.groups=TRUE)

# ROC(x.0,x.1,as.groups=TRUE)  # long output

ROC.plot(x.0,x.1,as.groups=TRUE)

ROC.plot(d.0,d.1,as.groups=TRUE,show.points=TRUE,flip.xaxis=TRUE,show.plot=FALSE,ci.region=FALSE,line.col="red")

z=seq(0,1,0.001)
sens.z <- pnorm(a+b*qnorm(z))
lines(z,sens.z,lty=1,col="firebrick")

##################################################

###
##
#