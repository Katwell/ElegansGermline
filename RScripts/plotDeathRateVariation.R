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
directories1=c("FastCellCycleNoCI",
			   "FastCellCycleNoCI02",
			   "FastCellCycleNoCI05")

#Results directories for plot 2. In this case, represents 3 different cell death rates with a cell cycle length of 24 hours
directories2=c("SlowCellCycleNoCI",
	           "SlowCellCycleNoCI01",
			   "SlowCellCycleNoCI02")


# Make plot 1"
par(mfrow=c(1,2),lwd=2, ps=10)
colors=c("darkgreen","orange","red")
for(i in seq(1,length(directories1),by=1)){
	
	data<-read.table(paste(resultsDirectory, directories1[i],"/GonadData.txt",sep=""),as.is=TRUE,header=FALSE,sep="\t");

	if(i==1){
		plot(data$V1+18.5, data$V5, xlab='Time (hours post-hatching)',type='l',main="8 hour adult cell cycle",
			 ylab="Number of proliferating cells", xlim=c(18,60), ylim=c(0,1100), col=colors[i] )
	}else{
		points(data$V1+18.5, data$V5, type='l',col=colors[i] )
	}
}
legend("topleft",legend=c("death rate p = 0.025", "death rate p = 0.2", "death rate p = 0.5"),fill=colors, bty = "n")
abline(h=250, lty=2, col="red") #Experimental data for comparison
mtext("A)", 3, at=15, padj=-0.5)


# Make plot 2
for(i in seq(1,length(directories2),by=1)){

	data<-read.table(paste(resultsDirectory, directories2[i], "/GonadData.txt",sep=""), as.is=TRUE,header=FALSE,sep="\t");

	if(i==1){
		plot(data$V1+18.5, data$V5, xlab='Time (hours post-hatching)',type='l',main="20 hour adult cell cycle",
			 ylab="Number of proliferating cells",xlim=c(18,90),ylim=c(0,800),col=colors[i] )
	}else{
		points(data$V1+18.5, data$V5, type='l', col=colors[i] )
	}
}
legend("topleft",legend=c("death rate p = 0.025", "death rate p = 0.1", "death rate p = 0.2"),fill=colors, bty = "n")
abline(h=250, lty=2, col="red") #Experimental data for comparison
mtext("B)", 3, at=15, padj=-0.5)


# UNCOMMENT HERE FOR POSTSCRIPT OUTPUT------------------------------------------------------------------
#dev.off()
#embed_fonts("./deathRateVariationPrelim.eps", outfile = "./deathRateVariation.eps",options = "-dEPSCrop")
#-------------------------------------------------------------------------------------------------------
