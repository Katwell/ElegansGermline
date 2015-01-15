#Plots gonadData.txt as in Figure 5, taking mean of 5 runs. NaN warning will be generated in the standard deviation
#calculation due to early time points where there is no first meiotic cell and so some data is undefined. These points 
#are ignored (not plotted).

#SETUP PLOT PARAMETERS
par(mfrow=c(3,3),lwd=2,cex.lab=1.2,cex.axis=1,mar=c(4.2,4.2,4.2,4.2))
total<-data.frame()
square<-data.frame()
std<-data.frame()

#LIST OF EXISTING FILE NUMBERS
fileNumbers = seq(1,5)
#BASE FILE NAME
fileName = "Baseline"



#LOOP THROUGH FILES, CALCULATE SUM AND SUM OF SQUARES FOR EACH DATA POINT
for(n in fileNumbers){

	fileCurrent = paste(fileName, n, sep="");
	dataCurrent<-read.table(paste(paste("/tmp/kathryn/testoutput/",fileCurrent,sep=""),"/GonadData.txt",sep=""), as.is=TRUE, header=FALSE, sep="\t");

	if(n==fileNumbers[1]){
		total = dataCurrent;
		square = dataCurrent^2;
	}else{
		total =  total + dataCurrent;
		square = square + dataCurrent^2;
	}
}

#CALCULATE MEAN AND STANDARD DEVIATION
data = total/length(fileNumbers);
std = sqrt((square/length(fileNumbers)) - data^2);
timepoints = (data$V1+18.5) # <- correct for fact that initial time is 18.5 hrs post-hatching


#PLOT TOTAL CELLS, +- 2 STDDEV
plot(0,0, xlab='Time (hph)', ylab="Total cells (count)", xlim=c(18.5,118.5), ylim=c(0,1100))                   # <- initialise plot
polygon(c(timepoints, rev(timepoints)), c(data$V7 + 2*std[,7], rev(data$V7 - 2*std[,7])), col=gray(0.4), lty=0) # <- confidence region
points(timepoints, data$V7, type='l')																		 # <- simulated mean
points(c(0,2.5,5.5,7.5,11.5,12.5,17,30.5,42.5,54.5)+18.5, c(16,40,63,110,154,177,308,700,750,900), type='l', col='red', lty=1) # <- expected result
mtext("a)", 3, at=7)


#PLOT PROLIFERATIVE CELLS, +- 2 STDDEV
plot(0,0, xlab='Time (hph)', ylab="Proliferative cells (count)", xlim=c(18.5,118.5), ylim=c(0,750)) 
polygon(c(timepoints,rev(timepoints)), c(data$V5 + 2*std[,5], rev(data$V5 - 2*std[,5])), col=gray(0.4),lty=0)
points(timepoints, data$V5, type='l')
points(c(0,2.5,5.5,7.5,11.5,12.5,17)+18.5, c(16,40,63,103,126,134,204), type='l', lty=1, col='red')
points(c(35.5,118), c(250,250), type='l', lty=2, col='dodgerblue')
mtext("b)", 3, at=7)


#PLOT SPERM COUNT, +- 2 STDDEV
plot(0,0, xlab='Time (hph)', ylab="Sperm cells (count)", xlim=c(18.5,118.5), ylim=c(0,170)) 
polygon(c(timepoints,rev(timepoints)), c(data$V4 + 2*std[,4], rev(data$V4 - 2*std[,4])), col=gray(0.4), lty=0)
points(timepoints, data$V4, type='l')
points(c(0,2.5,5.5,7.5,11.5,12.5,17,19.5)+18.5, c(0,0,0,0,4,41,85,140), type='l', lty=1, col='red')
points(seq(19.5,118,by=0.1)+18.5, 140-2.7*(seq(19.5,118,by=0.1)-19.5), type='l', lty=2, col='dodgerblue')
mtext("c)", 3, at=7)


#PLOT ORGAN LENGTH, +- 2 STDDEV
plot(0,0, xlab='Time (hph)', type='l', ylab=expression(paste("Length of gonad arm (", mu, "m)")), xlim=c(18.5,118.5), ylim =c(0,420))
polygon(c(timepoints,rev(timepoints)), c(data$V2 + 2*std[,2], rev(data$V2 - 2*std[,2])), col=gray(0.4), lty=0)
points(timepoints, data$V2, type='l')
points(c(18.5,22,26,31,35.5), c(32,61.7,91.4,198,405), lty=1, col="red", type='l')
points(c(35.5,118), c(405,405), lty=1, col="red", type='l')
mtext("d)",3, at=5)


#PLOT CELL ROW DISTANCES TO POINTS OF INTEREST, +- 2 STDDEV
#NOTE: the automated row counting algorithm in our code counts the endcap as one row. We therefore add 2 to the output
#to correct for this systematic error (the endcap usually contains about 3 cell rows by eye).
plot(0,0,xlab='Time (hph)',ylab=expression(paste("Last row containing a prolif. cell")),xlim=c(18.5,118),ylim=c(0,31))
polygon(c(timepoints,rev(timepoints)), c(data$V16+2 + 2*std[,16], rev(data$V16+2 - 2*std[,16])), col=gray(0.4), lty=0)
points(timepoints, data$V16 + 2,type='l')
points(c(35.5,118), c(28,28), type='l',lty=2,col='dodgerblue')
mtext("e)",3, at=4, adj=0)


#PLOT CELL ROW DISTANCES TO POINTS OF INTEREST, +- 2 STDDEV
#NOTE: the automated row counting algorithm in our code counts the endcap as one row. We therefore add 2 to the output
#to correct for this systematic error (the endcap usually contains about 3 cell rows by eye).
plot(0,0, xlab='Time (hph)', ylab="First row containing two meiotic cells", xlim=c(18.5,118), ylim=c(0,24))
polygon(c(timepoints,rev(timepoints)), c(data$V15+2 + 2*std[,15], rev(data$V15+2 - 2*std[,15])), col=gray(0.4), lty=0)
points(timepoints, data$V15+2, type='l')
points(c(35.5,118), c(19,19), type='l', lty=2, col='dodgerblue')
mtext("f)",3, at=5)


#PLOT MICRON DISTANCES TO POINTS OF INTEREST, +- 2 STDDEV. Commented out sections add an approximate 
#cell diameters axis.
plot(0,0, xlab='Time (hph)', ylab=expression(paste("Furthest mitotic cell from DTC (", mu, "m)")), xlim=c(18.5,118), ylim=c(0,180))
polygon(c(timepoints,rev(timepoints)), c(data$V8 +2*std[,8], rev(data$V8 - 2*std[,8])), col=gray(0.4),lty=0)
points(timepoints, data$V8, type='l')
points(c(35.5,118), c(66,66), type='l', lty=2, col='dodgerblue')
#par(new = T)
#CellDiam = 5.4;
#axis(side = 4, at=seq(0,210,by=50),  labels = round(seq(0,210,by=50)/CellDiam))
#mtext(side = 4, line = 3, "Cell diameters", cex =0.8)
mtext("g)", 3, at=5)


#PLOT MICRON DISTANCES TO POINTS OF INTEREST, +- 2 STDDEV
plot(0,0, xlab='Time (hph)', ylab=expression(paste("Closest meiotic cell to DTC (", mu, "m)")), xlim=c(18.5,118), ylim=c(0,110))
polygon(c(timepoints[data$V9<1000],rev(timepoints[data$V9<1000])), c(data$V9[data$V9<1000] +2*std[data$V9<1000,9], rev(data$V9[data$V9<1000] -2*std[data$V9<1000,9])),col=gray(0.4),lty=0)
points(timepoints[data$V9<1000], data$V9[data$V9<1000], type='l')
points(c(35.5,118), c(42,42), type='l', lty=2, col='dodgerblue')
#par(new = T)
#axis(side = 4, at=seq(0,110,by=20),  labels = round(seq(0,110,by=20)/CellDiam))
#mtext(side = 4, line = 3, "Cell diameters", cex =0.8)
mtext("h)",3, at=5)


#PLOT MICRON DISTANCES TO POINTS OF INTEREST, +- 2 STDDEV
plot(0,0, xlab='Time (hph)', ylab=expression(paste("Meiotic entry region length (", mu, "m)")), xlim=c(18.5,118), ylim=c(0,160))
variance = std^2;
combined=2*sqrt(variance[,8]-variance[,9])
polygon(c(timepoints[!is.nan(combined)],rev(timepoints[!is.nan(combined)])), c(data$V8[!is.nan(combined)]-data$V9[!is.nan(combined)] + combined[!is.nan(combined)], rev(data$V8[!is.nan(combined)]-data$V9[!is.nan(combined)] - combined[!is.nan(combined)])), col=gray(0.4),lty=0)
points(timepoints[data$V9<1000], data$V8[data$V9<1000] - data$V9[data$V9<1000],type='l')
points(c(35.5,118), c(24,24), type='l', lty=2, col='dodgerblue')
#par(new = T)
#axis(side = 4, at=seq(0,160,by=30),  labels = round(seq(0,160,by=30)/CellDiam))
#mtext(side = 4, line = 3, "Cell diameters",cex =0.8)
mtext("i)",3, at=4, adj=0)
