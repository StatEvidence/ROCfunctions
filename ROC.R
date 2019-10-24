ROC <- function(x, grp, as.groups=FALSE, digits=10, ci.lvl=0.95){

################################################
## Function:	ROC
## Version:		1.0
## Purpose:		Compute Empirical ROC points	
## 
## Author: 		JD Blume
## Date: 		June 2019
################################################

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

#### Seperate scores by group
	Neg <- d$score[status==0]
	Pos <- d$score[status==1]

	N.neg <- length(Neg)
	N.pos <- length(Pos)

#### Find unique scores (these are ROC point coordinates)
	Atoms <- unique(sort(d$score))

	N.atoms=length(Atoms)+1

	Score.range=range(d$score)
 
####   (i)   Compute TP and FP values for each distinct atom.
####   (ii)  Make sure that (0,0) and (1,1) lie on ROC.
####   (iii) Make sure that CDF has zero (and another 1) in its range.  
	FPC <- N.neg - cumsum(tabulate(
					match(Neg, c(Score.range[1]-10, Atoms)), 
						nbins = N.atoms))
    FP <- round(FPC/N.neg, digits=digits)
	
	TPC <- N.pos - cumsum(tabulate(
					match(Pos, c(Score.range[1]-10, Atoms)), 
						nbins=N.atoms))
    TP <- round(TPC/N.pos, digits=digits)

###	Compute Clopper-Pearson CIs for each ROC point (in x and y directions)

	FP.ci=sapply(FPC,
			function(x) binom.test(x, N.neg, conf.level=ci.lvl)$conf.int)

	TP.ci=sapply(TPC,
			function(x) binom.test(x, N.pos, conf.level=ci.lvl)$conf.int)

#### Return FPF, TPF, Cutpoint ; returns P(score > cutpoint | status)
keep <- list("fp"=FP, "tp"=TP, "cutpoint"=c(Atoms,tail(Atoms,n=1)+0.001), 
		"fpci"=round(FP.ci,4), "tpci"=round(TP.ci,4))

return(keep)
}

###
##
#