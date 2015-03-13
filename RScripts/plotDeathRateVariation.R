# Plots proliferative cell counts over time for two batches of results directories. Directories1 and directories2
# each get their own plot. The sample data we're plotting here represents three different cell death rates, with
# 2 different germ cell cycle lengths.

# UNCOMMENT HERE AND AT END OF FILE FOR POSTSCRIPT OUTPUT-----------------------------------------------
# library("extrafont")
# postscript("deathRateVariationPrelim.eps",width=6.83,height=4.2,family = "Arial", paper = "special", onefile = FALSE,horizontal = FALSE)
#-------------------------------------------------------------------------------------------------------

#USER INPUT: if not using our sample data, fill in the path to your results directories:
resultsDirectory = "../exampleOutput/"

#Results directories for plot 1. In this case, represents 3 different cell death rates with a cell cycle length of 8 hours
directories1=c("FastCycleD0.025",
			   "FastCycleD0.2",
			   "FastCycleD0.5")
endTime1 = 44
runNumbers1 = c("000","001","003","004","005") #exclude failed or slow runs

#Results directories for plot 2. In this case, represents 3 different cell death rates with a cell cycle length of 24 hours
directories2=c("SlowCycleD0.025",
	           "SlowCycleD0.1",
			   "SlowCycleD0.2")
endTime2 = 80 

runNumbers2 = c("000","001","002","003","004") #exclude failed or slow runs

#-------------------------------------------------------------------------------------------------------

makeplot = function(dir, endt, runs, colors, maint){
	
	for(i in seq(1,length(dir),by=1)){
	
		mean = c()
		square = c() 
		for( j in seq(1,length(runs),by=1) ){ 
			
			currentData<-read.table(paste(resultsDirectory, dir[i], "/output/",
				 runs[j],"/results_from_time_0/GonadData.txt",sep=""),as.is=TRUE,header=FALSE,sep="\t");

			if(length(mean) == 0){
				mean = currentData$V5[1:endt]
				square = currentData$V5[1:endt]*currentData$V5[1:endt]
			}else{
				mean = mean + currentData$V5[1:endt]
				square = square + currentData$V5[1:endt]*currentData$V5[1:endt]
			}
		}
		times = currentData$V1[1:endt] +18.5
		mean = mean/length(runs)
		std = sqrt(square/length(runs) - mean*mean) 
	
		if(i==1){
		plot(0, 0, xlab='Time (hours post-hatching)',type='l',main=maint,
				 ylab="Number of proliferating cells", xlim=c(18,max(times)), ylim=c(0,1.05*max(mean)))
		}
		points(times, mean, type='l',col=colors[i] )
	}
}



# Make plot 1"
par(mfrow=c(1,2),lwd=2, ps=10)
colors=c("darkgreen","orange","red")

makeplot(directories1, endTime1, runNumbers1, colors, "8 hour adult cell cycle")
legend("topleft",legend=c("death rate p = 0.025", "death rate p = 0.2", "death rate p = 0.5"),fill=colors, bty = "n")
abline(h=250, lty=2, col="red") #Experimental data for comparison
mtext("(A)", 3, at=15, padj=-0.5)

makeplot(directories2, endTime2, runNumbers2, colors, "24 hour adult cell cycle")
legend("topleft",legend=c("death rate p = 0.025", "death rate p = 0.1", "death rate p = 0.2"),fill=colors, bty = "n")
abline(h=250, lty=2, col="red") #Experimental data for comparison
mtext("(B)", 3, at=15, padj=-0.5)

# UNCOMMENT HERE FOR POSTSCRIPT OUTPUT------------------------------------------------------------------
#dev.off()
#embed_fonts("./deathRateVariationPrelim.eps", outfile = "./deathRateVariation.eps",options = "-dEPSCrop")
#-------------------------------------------------------------------------------------------------------
