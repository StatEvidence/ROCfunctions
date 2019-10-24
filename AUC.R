AUC <- function( x, grp, as.groups=FALSE) {

################################################
## Function:	AUC
## Version:		1.0
## Purpose:		Compute area under ROC curve
##				Compute CIs for AUC (approx and Robust)
## 
## Author: 		JD Blume
## Date: 		June 2019 (v1)
################################################

require(asht) # https://rdrr.io/cran/asht/man/wmwTest.html

if (as.groups==TRUE) {
	n <- 	c(length(x),length(grp))
	x <-	c(x,grp)
	grp <-  c(rep("a",n[1]),rep("b",n[1]))
			}

#### Errors / Warnings
if (length(levels(as.factor(grp)))!=2) {
	stop("Grouping variable is not dichotomous.")
}

#### Set Data Frame
status <- 1*(grp==levels(as.factor(grp))[2])

d <- as.data.frame(cbind("score" = x, "status" = status))
d <- d[complete.cases(d),]

#### Compute AUC
s.0 <- d$score[d$status==0]
s.1 <- d$score[d$status==1]

w.min <- wilcox.test(s.0, s.1, exact=FALSE)$statistic/length(s.0)/length(s.1)
area  <- max(w.min, 1-w.min)

#### Compute CI 
####	Uses max possible varaince over all continuous distributions ; CLT
var.upbd <- area*(1-area)/min(length(s.0),length(s.1))
ci.rob <- round(area + c(-1,1) * 1.96 * sqrt(var.upbd),4)
ci.rob <- pmin(ci.rob, 1)
ci.rob <- pmax(ci.rob, 0)

####	Mann-Whitney-U based CIs for AUC (Asht package)
wmw <- if (w.min > 0.5) {wmwTest(s.1, s.0)} else {wmwTest(s.0, s.1)}
ci.mwu <- round(wmw$conf.int,4)

keep <- list("auc"=round(area,4), "ci.mwu"=sort(ci.mwu), "ci.rob"=sort(ci.rob))

return(keep) 
}



###
##
#