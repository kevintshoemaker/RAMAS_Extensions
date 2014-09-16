#########################
###  R Script: Run Maxent and generate HS surface for reading into RAMAS Metapop
###  K. Shoemaker 8/2014

##################
###    INPUTS

###  To use this function ("MakeMaxentHSMap"), these inputs are required:
###    (1) OccData: a dataframe or SpatialPointsVector object ('sp' package) of occurrence locations for the species of interest
###              if a data frame, first col. must be the X coords and second column must be Y coords.
###    (2) NBackground: number of background points, or pseudoabsences, to generate within the study region (region defined by non-NA values of the predictor maps)
###    (3) Predictors: a "RasterStack" object ('raster' package) containing one or more maps of potential environmental predictors of occurrence (e.g., climate, landcover) 
###    (4) ModelNum: a model number: integer value, allows you to specify and identify multiple Maxent models (e.g., with different predictors)
###    (5) Directory: a file path indicating where the habitat suitability map should be stored  

###################
### OUTPUTS

###  Raster (ESRI GRID format (.ASC)) of habitat suitability across your study region, derived from a Maxent model. This file can be read directly into RAMAS GIS
###     This file is writted to the specified directory, and has the name (e.g.) HabMap1.asc. The "1" in the file name corresponds to the ModelNum (see model inputs, above) 


###################
###  NOTES

### To run this script, you must 
###       (1) download Maxent (http://www.cs.princeton.edu/~schapire/maxent/)
###       (2) copy the maxent .JAR file into the Dismo package folder: Put the file 'maxent.jar' in the 'java' folder of the 'dismo' package. 
###               That is the folder returned by system.file("java", package="dismo"). 
###               You need MaxEnt version 3.3.3b or higher. Please note that this program (maxent.jar) cannot be redistributed or used for commercial or for-profit purposes.
###               See the 'dismo' package documentation for more detailed instruction. 

### This script runs MaxEnt with default values for regularization, etc.  You can add any MaxEnt arguments you want to the call to Maxent- in the 'RunMaxEnt' function
###   See the 'dismo' package documentation for more details

### If you have questions, please contact the author of this script: kevin@ramas.com

###################
###  LOAD REQUIRED PACKAGES

library(raster)
library(rJava)      # NOTE: your Java build must match with your R build: e.g., if your verstion of Java is 32-bit, you must run the 32-bit version of R.
library(dismo) 

###################
### START WITH A FRESH R WORKSPACE

rm(list=ls())


############################
#####    LOAD FUNCTIONS

RunMaxEnt <- function(ModelNum,Directory,Predictors,pr=pr,bg=bg,plot=T){    # takes a while to run.   this also plots out the importance values...
		
	if(plot){
		graphics.off()
		par(ask=T)
	}
		  # checking if the MaxEnt jar file is present. If not, skip this bit
	jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
	
	objName <- paste("SDM",ModelNum,sep="")

	   # RUN MAXENT
	if (file.exists(jar)) {
	  eval(parse(text=paste(objName,"<- maxent(Predictors, pr,a=bg,removeDuplicates=T)",sep="")))
	  if(plot) eval(parse(text=paste("plot(",objName,")",sep="")))
	}else {
	  cat('cannot run because maxent is not available')
	  graphics.off()
	}
	
	setwd(Directory)
	filename <- paste(objName,".RData",sep="") 
	eval(parse(text=paste("save(",objName,",file=filename)",sep="")))    # save MaxEnt object to file
	
	  if(plot){
		plot(1)
		graphics.off()
	}

}


###########################
#######   READ IN MAXENT MODEL FOR VISUALIZATION

ReadMaxEnt <- function(ModelNum=ModelNum,Directory){
	ObjName <- paste("SDM",ModelNum,sep="")
	filename <- paste(ObjName,".RData",sep="")
	setwd(Directory)
	load(filename)
	MaxEntObj <<- eval(parse(text=ObjName))
}

############################
######  EVALUATE PERFORMANCE (e.g., auc, false positive rate)

GeneratePerformanceMetrics <- function(ModelNum,Directory,Predictors,pr=OccData,bg=bg){
	ReadMaxEnt(ModelNum=ModelNum,Directory)
	evalObj <- evaluate(OccData, bg, MaxEntObj, Predictors)   
	evalObj    # print out performance statistics

	setwd(Directory)
	filename <- paste("SDMPerformance_Model_",ModelNum,".RData",sep="")
	save(evalObj,file=filename)
}


##############################
########   GENERATE RASTER MAP OF MAXENT PREDICTIONS

      ## evaluate performance-  [ultimately use rigorous cross-validation!]
	  
GenerateHabMap <- function(ModelNum,Directory,evalObj,Predictors){     # this function takes a while to run for the whole study area.  

    par(ask=T)
						# read in MaxEnt object
	setwd(Directory)	
	ReadMaxEnt(ModelNum=ModelNum,Directory)		
 
	HabMap <- predict(Predictors, MaxEntObj, progress='')   # generate habitat map: takes a long time
	par(mfrow=c(1,2))
	plot(HabMap, main='Maxent, raw values')
	tr <- threshold(evalObj, 'spec_sens')
	plot(HabMap > tr, main='presence/absence')
	points(OccData, pch='+')
	
							# SAVE RASTER OBJECT
	
	setwd(Directory)
	filename <- paste("HabMap",ModelNum,".asc",sep="")
	writeRaster(HabMap,filename=filename,format="ascii",overwrite=T)
	
	cat(paste("habitat suitability raster has been saved as: ",Directory,"/",filename,"\n",sep=""))

	graphics.off()
}

##################################
#####   GENERATE BACKGROUND POINTS (pseudoabsences)

GenerateBGPoints <- function(Predictors, Directory,OccData,seed=1977,nPoints=1000){

	graphics.off()
    par(ask=T)	   

	set.seed(seed)
	mask <- reclassify(Predictors[[1]],rcl=c(-Inf,Inf,1))
	bg <- randomPoints(mask, nPoints)    # random points
	  
	plot(!is.na(mask), legend=FALSE)
	points(bg, cex=0.5)
	points(OccData,pch=20,col="red")
	
	    # SAVE BACKGROUND POINTS
	setwd(Directory)
	write.table(bg,file="BGPoints.csv",row.names=F,col.names=T,sep=",")
	
	plot(1)
	graphics.off()
}


MakeMaxentHSMap <- function(OccData,NBackground,Predictors,ModelNum,Directory){
	
    graphics.off()
    par(ask=T)
	
	####################
	####  LOAD VARIABLES FROM HARD DISK   
				
	setwd(Directory)  

	###################
	###### VISUALIZE PREDICTORS

	predictors
	names(predictors) 
	plot(predictors)
    plot(1)
	#####################
	####  PREPARE BACKGROUND POINTS (if necessary)

	GenerateBGPoints(Predictors,Directory,OccData,seed=1977,nPoints=NBackground)                                   		# generate pseudoabsences (only run if necessary)

	####################
	#### LOAD BACKGROUND POINTS
	   
	bg <- read.csv("BGPoints.csv",header=T)                           		# Read in Background Points


	###################
	#### RUN MAXENT
	
	RunMaxEnt(ModelNum=ModelNum,Directory,Predictors,pr=OccData,bg=bg,plot=T)                  		# Run MaxEnt and store model (only run if necessary)

	###################
	#### LOAD MAXENT MODEL

	ReadMaxEnt(ModelNum=ModelNum,Directory)


	###################
	#### GENERATE PERFORMANCE METRICS

	GeneratePerformanceMetrics(ModelNum,Directory,Predictors,pr=OccData,bg=bg)                		# generate performance metrics                          

	###################
	##### LOAD PERFORMANCE METRICS
			   
	filename <- paste("SDMPerformance_Model_",ModelNum,".RData",sep="")
	load(filename) 
	evalObj   

	####################
	#####  GENERATE HABITAT SUITABILITY SURFACE      (takes a while...)            			

	GenerateHabMap(ModelNum,Directory,evalObj,Predictors)                       

}




###############  TEST OF SCRIPT

 ###MakeMaxentHSMap <- function(OccData,NBackground,Predictors,ModelNum,Directory){

  ### load predictors
setwd("\\\\AB-Server\\backup\\Eagle GIS\\EagleGISData_810m\\") 
load("PredictorsForSDM.RData")                                   		# load predictor data
    
setwd("\\\\AB-Server\\backup\\Eagle GIS\\EagleGISData_810m\\Other\\DigitizedEagleData\\All nests")
OccData <- read.csv("ThinnedNests.csv",header=T)[,c(2,3)]                      		# pr stands for presence locations. bg stands for background locations
  
Directory <- "C:\\Users\\Kevin\\Desktop\\TEST"

NBackground=5000

ModelNum = 1

MakeMaxentHSMap(OccData,NBackground,Predictors=predictors,ModelNum,Directory)





