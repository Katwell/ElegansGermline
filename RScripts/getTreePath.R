# The function "getTreeOfPaths" defined here outputs a file suitable for opening in Paraview.
# The file contains the path taken by the cell with ID = startingCellID, as well as all its
# daughters, between startTime and endTime. Since cells are labelled in L3, when 16 germ
# cells are present, startingCellID ranges between 1 and 16.


getTreeOfPaths = function(startingCellID, startTime, endTime){

	# USER INPUT. Set the path to the directory containing your results folders:
	resultsDirectory = "../exampleOutput/"
	# USER INPUT. Name of the results folder you want to read from:
	baseFileName = "ForTracking"
	
	#READ IN A FILE TrackingData.txt CONTAINING CELL POSITION INFO
	pathDataFull<-read.table(paste(resultsDirectory, baseFileName, "/TrackingData.txt",sep=""),as.is=TRUE,header=FALSE,sep="\t");

	#READ IN A FILE DivisionData.txt CONTAINING CELL DIVISION INFO
	divisionData<-read.table(paste(resultsDirectory, baseFileName, "/DivisionData.txt",sep=""),as.is=TRUE,header=FALSE,sep="\t");


	# Crop data to time range of interest 
	croppedDivisionData = divisionData[(divisionData[,1]<=endTime),] 
	croppedPathData = pathDataFull[(pathDataFull[,1]<=endTime & pathDataFull[,1]>=startTime),]  
	measurementInterval = unique(croppedPathData[,1])[2]-unique(croppedPathData[,1])[1]



	# Find all the daughters of the selected cell
	searching = c(startingCellID);
	daughters = c(startingCellID);
	mothers =   c(startingCellID);
	iterations = 0
	while( length(searching)>0 & iterations<1000 ){
		
		# Look for daughters of the first cell in the searching list. Add those daughter to the search list.
		daughters = c(daughters, croppedDivisionData[ croppedDivisionData[,2]==searching[1],3] );
		searching = c(searching, croppedDivisionData[ croppedDivisionData[,2]==searching[1],3] );
		mothers   = c(mothers,   croppedDivisionData[ croppedDivisionData[,2]==searching[1],2] );
		
		# Remove front element of search list and iterate.
		if(length(searching)>1){
			searching = searching[2:length(searching)]
		}else{
			searching = c()
		}

		iterations = iterations+1
	}
	if(iterations >= 1000){
		print("WARNING: iteration limit exceeded while searching for cell daughters")
	}



	# Make a list of all the lines that should be drawn in the tree of paths:
	linesToDraw = data.frame();
	connections = c();
	lengthOfLines = c();
	timeData = c();
	directionData = data.frame();

	for(c in daughters){
		
		#Get the trail of this cell
		relevantPositions = croppedPathData[ croppedPathData[,2]==c, 3:5];
		nRows = dim(relevantPositions)[1] 										# <- number of datapoints we have for this cell
												
		currentRows = dim(linesToDraw)[1] 										# <- number of lines we have stored so far
		if(nRows > 1){ 															# <- we can only draw a line if we have more than one point on this cell's path
			
			linesToDraw = rbind(linesToDraw, relevantPositions[1:nRows,]);		# <- point on path
			connections = c(connections, (currentRows):(currentRows+nRows-1) ); # <- how to join them up
			lengthOfLines = c(lengthOfLines, nRows);
			timeData = c(timeData, croppedPathData[croppedPathData[,2]==c, 1]); # <- associate a time with each point on path
			directionData = rbind(directionData, relevantPositions[2:nRows,]-relevantPositions[1:(nRows-1),]); # <- associate a direction with each point

			# IF POSSIBLE, this block works out the cell's direction of motion at the point its path ends, by looking
			# at data recorded after the time period of interest
			times = croppedPathData[croppedPathData[,2]==c,1 ];
			nextTime = times[length(times)] + measurementInterval;
			nextPosition = pathDataFull[(pathDataFull[,2]==c & pathDataFull[,1]==nextTime) ,3:5];
			# If we were able to get a next position for the cell, after the end of the time interval
			if(dim(nextPosition)[1]>0){
				# Record cell's direction of travel
				directionData = rbind(directionData, nextPosition - relevantPositions[nRows,]); 
			}else{
				# Else record no direction
				directionData = rbind(directionData, c(0,0,0));
			}
		}
		

		# Adds in extra lines to represent division events, joining mother's position to daughter
		currentRows = dim(linesToDraw)[1]
		if(c!=startingCellID){ 																	 # <- For each daughter cell
		    							 
		    m = mothers[daughters==c];															 # <- Get mother's cell ID
			times = pathDataFull[pathDataFull[,2]==c,1];										 # <- Extract times at which position of the daughter cell is recorded
			birthTime = times[1];																 # <- First time recorded is closest to birth time
			p = croppedPathData[(croppedPathData[,2]==m & croppedPathData[,1]==birthTime),3:5];  # <- Position of the mother at birth time

			if(nRows>0 & dim(p)[1]>0){											# If we have both cell positions, add a line 
				linesToDraw = rbind(linesToDraw,p);								# From mother
				linesToDraw = rbind(linesToDraw,relevantPositions[1,]);			# To daughter
				connections = c(connections, (currentRows):(currentRows+1) );
				lengthOfLines = c(lengthOfLines, 2);							# Of length 2
				timeData = c(timeData, c(birthTime,birthTime));					# Created at time birth time.

				directionData = rbind(directionData, relevantPositions[1,]-p);  # Direction is mother to daughter
	
				# Again, IF POSSIBLE work out the direction the cell is travelling after it gets to its initial position  
				nextTime = birthTime + measurementInterval;
				nextPosition = pathDataFull[(pathDataFull[,2]==c & pathDataFull[,1]==birthTime) ,3:5];
				if(dim(nextPosition)[1]>0){
					directionData = rbind(directionData, nextPosition - relevantPositions[1,]);
				}else{
					directionData = rbind(directionData, c(0,0,0));
				}
			}	
		}
	}
	nPoints = dim(linesToDraw)[1];  # Total number of points, based on size of dataframe linesToDraw
	nLines = length(lengthOfLines); # Total number of lines



	# Write out the header for the vtk output file
	filename = paste("branchingPathsDataTime",startTime,"to",endTime,"Cell",startingCellID,".vtu",sep="")

	header = c("# vtk DataFile Version 2.0", 
		     paste("Branching path data time ",startTime," to ",endTime,", Cell ",startingCellID,sep=""),
		     "ASCII",
		     "DATASET POLYDATA",
		     paste("POINTS ", nPoints, " float",sep="")
		    )
    write(header, 
        file = filename,
		append = FALSE, 
		sep="\n"
	)

	# Write out the point data
	write(t(as.matrix(linesToDraw)), 
        file = filename,
		ncolumns = 3,
		append = TRUE, 
		sep="\t"
	)

    # Write out how to connect these points into lines
	bufSize = sum(lengthOfLines)+length(lengthOfLines)
	write(paste("LINES", nLines, bufSize, sep=" "),
		file = filename,
		append = TRUE, 
		sep="\n"
	)
	for(i in 1:length(lengthOfLines)){
		if(i==1){
			write(c(lengthOfLines[i], connections[1:sum(lengthOfLines[1:i])]), 
      		    file = filename,
				append = TRUE, 
				sep="\t"
			)
		}else{
			write(c(lengthOfLines[i], connections[(sum(lengthOfLines[1:i-1])+1):sum(lengthOfLines[1:i])]), 
      		    file = filename,
				append = TRUE, 
				sep="\t"
			)
		}
	}

	# Write out some extra data for visualisation (e.g. the time corresponding to each point, line directions)
	write(c(paste("POINT_DATA", nPoints, sep=" "),
                  "SCALARS time float 1",
		          "LOOKUP_TABLE pointsTable"), 
		file = filename,
		append = TRUE, 
		sep="\n"
	)
	write(timeData, 
		file = filename,
		append = TRUE, 
		sep="\n"
	)
	write(c("VECTORS directions float"), 
		file = filename,
		append = TRUE, 
		sep="\n"
	)
	write(t(as.matrix(directionData)),
		ncolumns = 3, 
		file = filename,
		append = TRUE, 
		sep="\t"
	)
	write(paste("LOOKUP_TABLE pointsTable ",nPoints,sep=""),
		file = filename,
		append = TRUE,
		sep = "\n"
	)
	write(t(as.matrix(linesToDraw)), 
        file = filename,
		ncolumns = 3,
		append = TRUE, 
		sep="\t"
	)
}