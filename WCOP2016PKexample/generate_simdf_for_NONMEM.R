#WCOP2016 simulation PK example: Geberate simulation data frame for NONMEM
	rm(list=ls(all=TRUE))
	graphics.off()

#Generate NONMEM data frame for the 2comp-2 transit first-order absorption model
	#Generate data frame: (if you havea a NONMEM df ready then just read and apply function)
	#Set dose records:

	df <- NULL
	dosetimes <- c(seq(0,48,12))
	tlast <- 96

	#Now define finer sample times for after a dose to capture Cmax
		doseseq <- c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,7,8,9,10)

#Use the outer product but with addition to expand this doseseq for all dosetimes
	PKtimes <- outer(dosetimes,doseseq,"+")

#set number of subjects
	nsub <- 1000
	ID <- 1:nsub

#Make dataframe
	df <- expand.grid("CID"=ID,"TIME"=sort(unique(c(seq(0,tlast,1),PKtimes))),"AMT"=0,"MDV"=0,"DV"=NA,"CLCR"=120)
	df$CLCR[df$TIME <= 48] <- 90

#Dose: It can be any arbitrary dose
	doserows <- subset(df, TIME%in%dosetimes)
	doserows$AMT <- 500
	doserows$MDV <- 1

#Add back dose information
	df <- rbind(df,doserows)
	df <- df[order(df$CID,df$TIME,df$AMT),]       # arrange df by TIME (ascending) and by AMT (descending)
	df <- subset(df, (TIME==0 & AMT==0)==F) # remove the row that has a TIME=0 and AMT=0

	write.csv(df, file="TwoCompFirstOrderAbsTwoTransit.csv", row.names = F, quote= F, na=".")

