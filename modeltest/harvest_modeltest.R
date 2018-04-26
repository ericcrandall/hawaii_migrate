
harvest.modeltest<-function(path){
	
	modelparams<-c("bestmodel","fA","fC","fG","fT","kappa","titv","rAC","rAG","rAT","rCG","rCT","rGT","alpha_G") #the parameter names
	paramvalues<-matrix(nrow=length(dir(path,full.names=F,pattern=".",no..=T)),ncol=length(modelparams)) #initialize a matrix with the right number of rows (directories) and columns (parameters)
	n<-1 #initialize position marker


	for(i in dir(path,full.names=F,pattern=".",no..=T)){
	#loop over all the directories with modeltest output, go into each directory, read the output down to the line with lots of dashes, then read each parameter value, stick these into the empty matrix in the appropriate position
		modeltestout<-scan(file=paste(path,i,"modeltest_output.txt",sep="/"),skip=200,what="character") #make the appropriate path
		which(modeltestout=="-------------------------") #find the line names with this number of - on them
		pv<-as.numeric(modeltestout[seq(which(modeltestout=="-------------------------")[1]+2,which(modeltestout=="-------------------------")[2]-1,by=2)]) #create a vector of every other line between the two lines with ---- on them, and grab the corresponding values
	
		modl<-modeltestout[grep("selected:", modeltestout)+3] #this should find the line that comes after the last instance of BIC, which will be the best model
		pv<-c(modl,pv)
	
		paramvalues[n,]<-pv #stick them into the matrix at the appropriate position

		n<-n+1 #update the position marker by 1
}

	paramvalues<-data.frame(paramvalues, row.names=dir(path,full.names=F,pattern=".",no..=T)) #make a dataframe with the rows and columns named
	colnames(paramvalues)<-modelparams
	return(paramvalues)
}


 	
	









	


	


	








	
	


	