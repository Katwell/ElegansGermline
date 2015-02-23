postscript("Rates.eps",width=6,height=6,horizontal=FALSE, pointsize=12, family="sans")

#Calculate rates
timeBetweenMeasurements = c(3.5, 4, 5, 4.5)
times = c(18.5,22,22,26,26,31,31,35.5)
totalLinearGrowth = c(61.7-32, 91.4-61.7, 198-91.4, 405-198)/timeBetweenMeasurements
radialGrowth = c(0, 6.08-5.48, 7.97-6.08, 11.3-7.97)/timeBetweenMeasurements
stretching = c(0, 0, 0, 126-53)/timeBetweenMeasurements

#Calculate combined standard deviation for each rate, based on individual measurement standard deviations
totalLinearGrowthSD = c(9.33,   sqrt((9.33^2)+(7.76^2)),      sqrt((23.3^2)+(7.76^2)),   	sqrt((23.3^2)+(46.7^2)))/timeBetweenMeasurements 
radialGrowthSD =      c(0, 		0, 							  sqrt((0.75^2)+(0.561^2)),	    sqrt((1.54^2)+(0.75^2)))/timeBetweenMeasurements
stretchingSD =        c(0, 		0,							  0,							sqrt((21.5^2)+(12.3^2)))/timeBetweenMeasurements
		


#plus and minus one standard deviation
totalLinearGrowthpSD = totalLinearGrowth + totalLinearGrowthSD 
radialGrowthpSD = radialGrowth + radialGrowthSD
stretchingpSD = stretching + stretchingSD                                                    															
totalLinearGrowthmSD = totalLinearGrowth - totalLinearGrowthSD 
radialGrowthmSD = radialGrowth - radialGrowthSD
stretchingmSD = stretching - stretchingSD
#limit ranges at zero
totalLinearGrowthmSD[totalLinearGrowthmSD<0]=0
radialGrowthmSD[radialGrowthmSD<0]=0
stretchingmSD[stretchingmSD<0]=0

 
#plotting 
par(lwd = 2)
plot(times,totalLinearGrowth[c(1,1,2,2,3,3,4,4)],type='l',col='firebrick2', xlab = 'Time (hours post-hatching)', ylab = expression(paste("Rate (", mu, "m/hr)")),ylim=c(0,70))
polygon(c(times,rev(times)), c(totalLinearGrowthpSD[c(1,1,2,2,3,3,4,4)],rev(totalLinearGrowthmSD)[c(1,1,2,2,3,3,4,4)]), col="lightcoral", lty=0)
polygon(c(times,rev(times)), c(radialGrowthpSD[c(1,1,2,2,3,3,4,4)],rev(radialGrowthmSD)[c(1,1,2,2,3,3,4,4)]), col="mediumspringgreen",lty=0)
polygon(c(times,rev(times)), c(stretchingpSD[c(1,1,2,2,3,3,4,4)],rev(stretchingmSD)[c(1,1,2,2,3,3,4,4)]), col="lightskyblue1",lty=0)
points(times,stretching[c(1,1,2,2,3,3,4,4)],type='l',col='dodgerblue2')
points(times,radialGrowth[c(1,1,2,2,3,3,4,4)],type='l',col='forestgreen')
points(times,totalLinearGrowth[c(1,1,2,2,3,3,4,4)],type='l',col='firebrick2')

legend("topleft",legend = c("Total linear growth","Radial growth", "Proximal arm growth"),fill=c("firebrick2","forestgreen","dodgerblue2"),bty="n")


#annotate borders of each growth phase
ytop = 55
ytext = 55
segments(18.5,0,18.5,ytop,lty = 3,lwd=1)
segments(22,0,22,ytop,lty = 3,lwd=1)
segments(26,0,26,ytop,lty = 3,lwd=1)
segments(31,0,31,ytop,lty = 3,lwd=1)
segments(35.5,0,35.5,ytop,lty = 3,lwd=1)
text(x=20.2,y=ytext,labels="Early L3")
text(x=23.9,y=ytext,labels="Late L3")
text(x=28.6,y=ytext,labels="Early L4")
text(x=33.4,y=ytext,labels="Late L4")

dev.off()