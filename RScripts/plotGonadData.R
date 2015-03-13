# Plots some of the output of a C. elegans germline simulation, as recorded in the file gonadData.txt
#
# A NaN warning will be generated when this script runs, because the distance to the first meiotic cell
# is undefined before meiotic cells appear. Usually the warning can be safely ignored.

# UNCOMMENT HERE AND AT END OF FILE FOR POSTSCRIPT OUTPUT-----------------------------------------------
# library("extrafont")
# postscript("GonadDataPrelim.eps", height = 4.5, width = 6.83, family = "Arial", paper = "special", onefile = FALSE, horizontal = FALSE)
#------------------------------------------------------------------------------------------------------


#USER INPUT - FILL IN THESE FIELDS:
resultsDirectory = "../exampleOutput/production_2.8/"
plotMultipleFiles = TRUE				# If true we plot the mean and standard deviation of gonadData.txt 
										# files from a collection of directories
fileNumbers = seq(0,29)					# Input directories should be located under resultsDirectory and
										# should be named <baseFileName> + <number>
baseFileName = "" 
useManuallyCountedCellRows = TRUE		# If true, manual cell row count data is expected to be 
										# provided in each results directory, in a file called 
										# ZoneLengths.txt (see our example data)
useManuallyCountedProlifCells = FALSE 	# If true, manual proliferative cell counts are expected 
										# to be provided in each results directory, in a file called 
										# CorrectedProlifCounts.txt (see our example data)
showMicronZoneLengths = FALSE 			# Whether to output the position of certain significant cells
										# in microns from the DTC, in addition to cell rows
										# from the DTC

nStandardDeviations = 1 				# Number of standard deviations to shade either side of the mean
MAXTIME = 82							# Only plot simulations with data up to time 82. Useful for
										# visualising incomplete sets of runs.		
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------



# SETUP GRAPH PROPERTIES AND GRAPH PLOTTING FUNCTIONS
if(showMicronZoneLengths){
	rowsOfGraphs = 1
}else{
	rowsOfGraphs = 2
}
par(mfrow=c(rowsOfGraphs,3), lwd=2, ps=12, bg="white", mar=c(5,5,5,5))

plotLarvalStageMarkers = function(ytop){
	polygon(c(31.5,35.5,35.5,31.5), c(0,0,ytop,ytop), col=gray(0.8), lty=0)
	segments(18.5,0,18.5, ytop, lty = 3,lwd = 1)
	segments(26,0,26, ytop, lty = 3,lwd = 1)
	segments(35.5,0,35.5, ytop, lty = 3,lwd = 1)
}



# LOOP THROUGH RESULTS FILES AND CALCULATE A MEAN AND STANDARD DEVIATION FOR EACH PIECE OF DATA
sum       <-data.frame()
sumSquares<-data.frame()
count =0 
if(!plotMultipleFiles){
	fileNumbers = 0;
	meandata = read.table(paste(resultsDirectory,baseFileName,"/results_from_time_0/GonadData.txt",sep=""), as.is=TRUE, header=FALSE, sep="\t");
	stddevdata = meandata - meandata #zero this out
}else{
	for(n in fileNumbers){
	
		dataCurrent<-read.table(paste(resultsDirectory,baseFileName,n,"/results_from_time_0/GonadData.txt",sep=""), as.is=TRUE, header=FALSE, sep="\t");

		print(length(dataCurrent$V1))
		dataCurrent=dataCurrent[dataCurrent$V1 <= MAXTIME,]
		
		if(sum(dataCurrent$V1<=MAXTIME)<MAXTIME){
		}else{
			count = count+1
			if(length(sum)==0){
				sum = dataCurrent;
				sumSquares = dataCurrent^2;
				timepoints = (dataCurrent$V1 + 18.5) # <- corrects for the fact that the simulation starts 
										 # at 18.5 hrs post-hatching
			}else{
				sum = sum + dataCurrent;
				sumSquares = sumSquares + dataCurrent^2;
			}
		}
	}
	meandata = sum / count;#length(fileNumbers);
	stddevdata = sqrt((sumSquares/count) - meandata^2)#length(fileNumbers)) - meandata^2 );
	print("No. runs = ")
	print(count)
}


# READ IN MANUALLY COUNTED CELL ROW MEASUREMENTS IF REQUESTED. THE ALTERNATIVE IS TO TRUST THE SIMULATION'S
# BUILT IN CELL ROW COUNTING ALGORITHM (USUALLY UNDERESTIMATES THE TRUE ROW NUMBER)
if(useManuallyCountedCellRows){
	if(!plotMultipleFiles){
		rowData <- read.table(paste(resultsDirectory,baseFileName,"/ZoneLengths.txt",sep=""), as.is=TRUE, header=FALSE, sep="\t");
		rowTimePoints = rowData[,1]
		rowMeioticMean = rowData[,2]
		rowMitoticMean = rowData[,3]
		rowMeioticStd = rowMeioticMean - rowMeioticMean
		rowMitoticStd = rowMitoticMean - rowMitoticMean
	}else{
		for(n in fileNumbers){
			rowDataCurrent <- read.table(paste(resultsDirectory,baseFileName,n,"/ZoneLengths.txt",sep=""), as.is=TRUE, header=FALSE, sep="\t");
			rowTimePoints = rowDataCurrent[,1] + 18.5
			if(n == fileNumbers[1]){
				rowMeioticMean = rowDataCurrent[,2]
				rowMitoticMean = rowDataCurrent[,3]
				rowMeioticSquare = rowDataCurrent[,2]*rowDataCurrent[,2]
				rowMitoticSquare = rowDataCurrent[,3]*rowDataCurrent[,3]	
			}else{
				rowMeioticMean = rowMeioticMean + rowDataCurrent[,2]
				rowMitoticMean = rowMitoticMean + rowDataCurrent[,3]
				rowMeioticSquare = rowMeioticSquare + rowDataCurrent[,2]*rowDataCurrent[,2]
				rowMitoticSquare = rowMitoticSquare + rowDataCurrent[,3]*rowDataCurrent[,3]
			}
		}
		rowMeioticMean = rowMeioticMean/length(fileNumbers)
		rowMitoticMean = rowMitoticMean/length(fileNumbers)
		rowMeioticStd = sqrt((rowMeioticSquare/length(fileNumbers)) - rowMeioticMean*rowMeioticMean)
		rowMitoticStd = sqrt((rowMitoticSquare/length(fileNumbers)) - rowMitoticMean*rowMitoticMean)	
	}
}



# READ IN MANUAL PROLIFERATIVE CELL COUNTS IF REQUESTED. THE SIMULATION COUNTS ALL NON-MEIOTIC CELLS AS 
# PROLIFERATIVE BY DEFAULT. MANUAL COUNTS MADE IN PARAVIEW CAN BE PROVIDED INSTEAD, FOR EXAMPLE, TO BETTER
# MIMIC THE WAY PROLIFERATIVE CELLS ARE COUNTED EXPERIMENTALLY (counting stops after first row containing 
# 2 meiotic cells)
if(useManuallyCountedProlifCells){
	if(!plotMultipleFiles){
		cellDataMean <- read.table(paste(resultsDirectory,baseFileName,"/CorrectedProlifCounts.txt",sep=""), as.is=TRUE, header=FALSE, sep="\t");
		cellTimePoints = cellDataMean$V1
		cellDataStd = cellDataMean - cellDataMean
	}else{
		for(n in fileNumbers){
			cellDataCurrent <- read.table(paste(resultsDirectory,baseFileName,n,"/CorrectedProlifCounts.txt",sep=""), as.is=TRUE, header=FALSE, sep="\t");
			cellTimePoints = cellDataCurrent$V1 + 18.5
			if(n == fileNumbers[1]){
				cellDataMean = cellDataCurrent$V2
				cellDataSquare = (cellDataCurrent$V2)*(cellDataCurrent$V2)		
			}else{
				cellDataMean = cellDataMean+cellDataCurrent$V2
				cellDataSquare = cellDataSquare+(cellDataCurrent$V2)*(cellDataCurrent$V2)
			}
		}
		cellDataMean = cellDataMean/length(fileNumbers)
		cellDataStd = sqrt((cellDataSquare/length(fileNumbers)) - cellDataMean*cellDataMean)	
	}
}



if(!showMicronZoneLengths){


#PLOT TOTAL CELLS, +- 2 STDDEV
ytop = max(meandata$V7 + nStandardDeviations*stddevdata[,7])*1.05 # <- max y
plot(0,0, xlab='Time (hph)', ylab="Total cells (count)", xlim=c(18.5,100), ylim=c(0,ytop)) # <- initialise plot
plotLarvalStageMarkers(ytop)

polygon(c(timepoints, rev(timepoints)), c(meandata$V7 + nStandardDeviations*stddevdata[,7], rev(meandata$V7 - nStandardDeviations*stddevdata[,7]) ), col=gray(0.4), lty=0) # <- confidence region
points(timepoints, meandata$V7, type='l')																   # <- simulated mean
points(c(0,2.5,5.5,7.5,11.5,12.5,17)+18.5, c(16,40,63,110,154,177,308), type='l', col='firebrick2', lty=1) # <- experimental result
mtext("(A)", 3, at=15, padj=-0.5)

print(timepoints)
print(meandata$V7)


#PLOT PROLIFERATIVE CELLS, +- 2 STDDEV
if(useManuallyCountedProlifCells){
	ytop = max(cellDataMean + nStandardDeviations*cellDataStd)*1.05
	plot(0,0, xlab='Time (hph)', ylab="Proliferative cells (count)", xlim=c(18.5,100), ylim=c(0,ytop)) 
	plotLarvalStageMarkers(ytop)
	polygon(c(cellTimePoints, rev(cellTimePoints)), c(cellDataMean + nStandardDeviations*cellDataStd, rev(cellDataMean - nStandardDeviations*cellDataStd)), col=gray(0.4),lty=0)
	points(cellTimePoints, cellDataMean, type='l')
}else{
	ytop = max(meandata$V5 + nStandardDeviations*stddevdata[,5])*1.05
	plot(0,0, xlab='Time (hph)', ylab="Proliferative cells (count)", xlim=c(18.5,100), ylim=c(0,ytop)) 
	plotLarvalStageMarkers(ytop)
	polygon(c(timepoints,rev(timepoints)), c(meandata$V5 + nStandardDeviations*stddevdata[,5], rev(meandata$V5 - nStandardDeviations*stddevdata[,5])), col=gray(0.4),lty=0)
	points(timepoints, meandata$V5, type='l')
}
points(c(0,2.5,5.5,7.5,11.5,12.5,17)+18.5, c(16,40,63,103,126,134,204), type='l', lty=1, col='firebrick2')
points(c(35.5,41.5,48.5,62.5,91.5), c(201,244,251,220,157), type='l', lty=1, col='firebrick2')
mtext("(B)", 3, at=15, padj=-0.5)


#PLOT SPERM COUNT, +- 2 STDDEV
ytop = max((meandata$V4 + nStandardDeviations*stddevdata[,4]),140)*1.05
plot(0,0, xlab='Time (hph)', ylab="Sperm cells (count)", xlim=c(18.5,100), ylim=c(0,ytop)) 
plotLarvalStageMarkers(ytop)
polygon(c(timepoints,rev(timepoints)), c(meandata$V4 + nStandardDeviations*stddevdata[,4], rev(meandata$V4 - nStandardDeviations*stddevdata[,4])), col=gray(0.4), lty=0)
points(timepoints, meandata$V4, type='l')
points(c(0,2.5,5.5,7.5,11.5,12.5,17,19.5)+18.5, c(0,0,0,0,4,41,85,140), type='l', lty=1, col='firebrick2')
points(seq(19.5,118,by=0.1)+18.5, 140-2.7*(seq(19.5,118,by=0.1)-19.5), type='l', lty=2, col='firebrick2')
mtext("(C)", 3, at=15, padj=-0.5)

print(meandata$V4)


#PLOT ORGAN LENGTH, +- 2 STDDEV
ytop = max((meandata$V2 + nStandardDeviations*stddevdata[,2]),405)*1.05
plot(0,0, xlab='Time (hph)', type='l', ylab=expression(paste("Length of gonad arm (", mu, "m)")), xlim=c(18.5,118.5), ylim =c(0,ytop))
plotLarvalStageMarkers(ytop)
polygon(c(timepoints,rev(timepoints)), c(meandata$V2 + nStandardDeviations*stddevdata[,2], rev(meandata$V2 - nStandardDeviations*stddevdata[,2])), col=gray(0.4), lty=0)
points(timepoints, meandata$V2, type='l')
points(c(18.5,22,26,31,35.5), c(32,61.7,91.4,198,405), lty=1, col="firebrick2", type='l')
points(c(35.5,118), c(405,405), lty=2, col="firebrick2", type='l')
mtext("(D)",3, at=15, padj=-0.5)


# PLOT CELL ROW DISTANCES TO POINTS OF INTEREST, +- 2 STDDEV
# NOTE: the automated row counting algorithm in our code counts the endcap as one row. We therefore add 2 to the output
# to correct for this systematic error (the endcap usually contains about 3 cell rows by eye).
if(useManuallyCountedCellRows){
	ytop = 35#max(rowMeioticMean + nStandardDeviations*rowMeioticStd)*1.05
	plot(0,0, xlab='Time (hph)', ylab="Distance to 1st row containing\n >1 meiotic cell (CD)", xlim=c(18.5,100), ylim=c(0,ytop))
	plotLarvalStageMarkers(ytop)
    polygon(c(rowTimePoints,rev(rowTimePoints)), c(rowMeioticMean-nStandardDeviations*rowMeioticStd, rev(rowMeioticMean+nStandardDeviations*rowMeioticStd)), col=gray(0.4), lty=0)
    points(rowTimePoints, rowMeioticMean, type='l')
}else{
	ytop = max(meandata$V15[meandata$V1>10] + stddevdata[,15][meandata$V1>10])*1.05
	plot(0,0, xlab='Time (hph)', ylab="1st row containing two meiotic cells", xlim=c(18.5,100), ylim=c(0,ytop))
	plotLarvalStageMarkers(ytop)
    polygon(c(timepoints,rev(timepoints)), c(meandata$V15+2 + nStandardDeviations*stddevdata[,15], rev(meandata$V15+2 - nStandardDeviations*stddevdata[,15])), col=gray(0.4), lty=0)
    points(timepoints, meandata$V15+2, type='l')
}
points(c(35.5,118), c(20,20), type='l', lty=2, col='firebrick2')
mtext("(E)",3, at=15, padj=-0.5)


#PLOT CELL ROW DISTANCES TO POINTS OF INTEREST, +- 2 STDDEV
#NOTE: the automated row counting algorithm in our code counts the endcap as one row. We therefore add 2 to the output
#to correct for this systematic error (the endcap usually contains about 3 cell rows by eye).
if(useManuallyCountedCellRows){
	ytop = 35#max(rowMitoticMean + nStandardDeviations*rowMitoticStd)*1.05
	plot(0,0,xlab='Time (hph)',ylab=expression(paste("Distance to last row containing\n a proliferative cell (CD)")),xlim=c(18.5,100),ylim=c(0,ytop))
	plotLarvalStageMarkers(ytop)
    polygon(c(rowTimePoints, rev(rowTimePoints)), c(rowMitoticMean-nStandardDeviations*rowMitoticStd,rev(rowMitoticMean+nStandardDeviations*rowMitoticStd)), col=gray(0.4), lty=0)
    points(rowTimePoints, rowMitoticMean, type='l')
}else{
	ytop = max(meandata$V16 + nStandardDeviations*stddevdata[,16])*1.05
	plot(0,0,xlab='Time (hph)',ylab=expression(paste("Last row containing a prolif. cell")),xlim=c(18.5,100),ylim=c(0,ytop))
	plotLarvalStageMarkers(ytop)
    polygon(c(timepoints,rev(timepoints)), c(meandata$V16+2 + nStandardDeviations*stddevdata[,16], rev(meandata$V16+2 - nStandardDeviations*stddevdata[,16])), col=gray(0.4), lty=0)
    points(timepoints, meandata$V16+2, type='l')
}
points(c(35.5,118), c(28,28), type='l',lty=2,col='firebrick2')
mtext("(F)",3, at=15, padj=-0.5)

}else{

#PLOT MICRON DISTANCES TO POINTS OF INTEREST, +- 2 STDDEV
ytop = max(meandata$V8 + nStandardDeviations*stddevdata[,8])*1.05
plot(0,0, ann=FALSE, xlab='Time (hph)', ylab=expression(paste("Furthest mitotic cell from DTC (", mu, "m)")), xlim=c(18.5,100), ylim=c(0,ytop))
plotLarvalStageMarkers(ytop)
polygon(c(timepoints,rev(timepoints)), c(meandata$V8 + nStandardDeviations*stddevdata[,8], rev(meandata$V8 - nStandardDeviations*stddevdata[,8])), col=gray(0.4),lty=0)
points(timepoints, meandata$V8, type='l')
points(c(35.5,118), c(70,70), type='l', lty=2, col='firebrick2')
mtext(side = 2, line = 3, expression(paste("Furthest mitotic cell from DTC (", mu, "m)")), padj=1, cex=0.7)
mtext(side = 1, line = 3, "Time (hph)", padj=-0.5, cex=0.7)
mtext("(A)", 3, at=15, padj=-0.5)


#PLOT MICRON DISTANCES TO POINTS OF INTEREST, +- 2 STDDEV
ytop = max(meandata$V9[meandata$V9<1000] + nStandardDeviations*stddevdata[,9][meandata$V9<1000])*1.05
plot(0,0, ann=FALSE, xlab='Time (hph)', ylab=expression(paste("Closest meiotic cell to DTC (", mu, "m)")), xlim=c(18.5,100), ylim=c(0,ytop))
plotLarvalStageMarkers(ytop)
polygon(c(timepoints[meandata$V9<1000],rev(timepoints[meandata$V9<1000])), c(meandata$V9[meandata$V9<1000] + nStandardDeviations*stddevdata[meandata$V9<1000,9],
        rev(meandata$V9[meandata$V9<1000] - nStandardDeviations*stddevdata[meandata$V9<1000,9])),col=gray(0.4),lty=0)
points(timepoints[meandata$V9<1000], meandata$V9[meandata$V9<1000], type='l')
points(c(35.5,118), c(40,40), type='l', lty=2, col='firebrick2')
mtext(side = 2, line = 3, expression(paste("Closest meiotic cell to DTC (", mu, "m)")), padj=1, cex=0.7)
mtext(side = 1, line = 3, "Time (hph)", padj=-0.5, cex=0.7)
mtext("(B)",3, at=15, padj=-0.5)


#PLOT MICRON DISTANCES TO POINTS OF INTEREST, +- 2 STDDEV
ytop = max(meandata$V8[meandata$V9<1000] - meandata$V9[meandata$V9<1000])*1.05
plot(0,0, ann=FALSE, xlab='Time (hph)', ylab=expression(paste("Meiotic entry region length (", mu, "m)")), xlim=c(18.5,100), ylim=c(0,ytop))
variance = stddevdata^2;
combined=nStandardDeviations*sqrt(variance[,8] + variance[,9])
plotLarvalStageMarkers(ytop)
polygon(c(timepoints[!is.nan(combined)],rev(timepoints[!is.nan(combined)])),
        c(meandata$V8[!is.nan(combined)]-meandata$V9[!is.nan(combined)] + combined[!is.nan(combined)], 
          rev(meandata$V8[!is.nan(combined)]-meandata$V9[!is.nan(combined)] - combined[!is.nan(combined)])), col=gray(0.4),lty=0)
points(timepoints[meandata$V9<1000], meandata$V8[meandata$V9<1000] - meandata$V9[meandata$V9<1000],type='l')
points(c(35.5,118), c(30,30), type='l', lty=2, col='firebrick2')
mtext(side = 2, line = 3, expression(paste("Meiotic entry region length (", mu, "m)")), padj=1, cex=0.7)
mtext(side = 1, line = 3, "Time (hph)", padj=-0.5, cex=0.7)
mtext("(C)",3, at=15, padj=-0.5)

}

# UNCOMMENT FOR POSTSCRIPT OUTPUT-----------------------------------------------------------------------
# dev.off()
# embed_fonts("./GonadDataPrelim.eps", outfile = "./GonadData.eps",options = "-dEPSCrop")
#------------------------------------------------------------------------------------------------------
