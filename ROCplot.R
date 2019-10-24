ROC.plot <- function(x, grp, as.groups=FALSE, digits=10,
			line.col='black', line.type=1, points.col='black', points.pch=20,
			show.plot=TRUE, show.lines=TRUE, show.points=FALSE,
			flip.xaxis = FALSE, ci.region=TRUE, ci.border=NA, ci.color="aliceblue"){

################################################
## Function:	ROC.plot
## Version:		1.0
## Purpose:		Plot ROC curve with CI Region
##
## Dependecies: ROC function
## 
## Author: 		JD Blume
## Date: 		June 2019
################################################

#### Get ROC Data
roc.pts <- ROC(x=x, grp=grp, as.groups=as.groups, digits=digits)

#### Get CI limits for ROC curve
## 	 Compute maximum CI region from CI limits above the cruve
top <- data.frame(rbind(cbind(x=roc.pts$fpci[1,],y=roc.pts$tp),
				 cbind(x=roc.pts$fp,y=roc.pts$tpci[2,])))

top <- 		top[order(top[,'x'],top[,'y']),]
top$ymax <- cummax(top[,'y'])

##   Compute minimum CI region from CI limits below the cruve
bot <- data.frame(rbind(cbind(x=roc.pts$fpci[2,],y=roc.pts$tp),
					 cbind(x=roc.pts$fp,y=roc.pts$tpci[1,])))
bot <- bot[order(-bot[,'x'],-bot[,'y']),]
bot$ymin <- cummin(bot[,'y'])

##   Compute ROC CI region using polygon function
region <- data.frame(x=c(top$x, bot$x), y=c(top$ymax, bot$ymin))

#### Plotting 
if (show.plot==TRUE) {
	plot(roc.pts$fp,roc.pts$tp,type="n",pty="s",xaxt="n",yaxt="n",
	ylab="Sensitivity", xlab="")
	axis(side=2,at=seq(0,1,0.2),labels=format(seq(0,1,0.2),digits=2),las=1)
		
if (flip.xaxis==TRUE) {
	axis(side=1,at=seq(0,1,0.2),labels=format(seq(1,0,-0.2),digits=2))
	mtext(side=1,line=3,"Specificity")
		}		
		else{ axis(side=1); mtext(side=1, line=3, "1-Specificity")}

		}
		
if (ci.region==TRUE) {
polygon(region$x, region$y, col=ci.color, border=ci.border)
		}

if (show.plot==TRUE) {abline(0,1,lty=2,lwd=0.5)}

if (show.points==TRUE) {
	points(roc.pts$fp,roc.pts$tp, pch=points.pch, col=points.col)
		}

if (show.lines==TRUE) {
	lines(roc.pts$fp,roc.pts$tp, lty=line.type, col=line.col)
		}
}

###
##
#