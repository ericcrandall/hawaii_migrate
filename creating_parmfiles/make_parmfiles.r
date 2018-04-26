#internal function ripped off from the phangorn package following Yang 1994
#discrete.gamma<-function (alpha, k) {
#    if (k == 1) 
#        return(1)
#    quants <- qgamma((1:(k - 1))/k, shape = alpha, rate = alpha)
#    diff(c(0, pgamma(quants * alpha, alpha + 1), 1)) * k
#}
#phangorn:::discrete.gama(0.5,4)
# Why Migrate doesn't do it this way according to Peter Beerli:Yang uses means or medians of categories that have equal probabilities, migrate (and phylip)
#calculate the more sensible Generalized Laguerre quadrature points  (Laguerre quadrature is optimal for the gamma distribution) this leads to means of categories that do have #different probabilities (you see that in the table migrate presents). For small to relatively large alphas migrate uses the generalized Laguerre quadrature and for very large alphas I #use Hermite’s quadrature.
#Felsenstein, J. 2001. Taking variation of evolutionary rates between sites into account in inferring phylogenies. Journal of Molecular Evolution 53: 447-455.


#############SETTINGS#############################

#path 
path<-"~/Dropbox/Crandall_tobo"

#Starting values - can also read in matrices
start.theta<-"0.01"
start.migration<-"100000"

#Priors
theta.prior<-" 0.000010 0.010000 10.000000 1.00000 "
mig.prior<-" 0.000100 100000.000000 1000000.000000 100000.000000 "

#Chain-Control
increment<-"500"
chainlength<-"50000"
burnin<-"0"
replicate<-"replicate=YES:3" #or "replicate=NO"

###########READ IN DATA#########################

#read in modeltest output for DNA sequence models
modeltest<-read.delim("~/Dropbox/Crandall_tobo/modeltest_output_table2.txt",as.is=T)

#read in sample size matrix, and edit it a bit to remove certain pops, and lump others
sampleN<-read.csv("~/Dropbox/Crandall_tobo/hawaii_sample_sizes.csv",stringsAsFactors=F)

sampleN[is.na(sampleN)]<-0 #replace NAs with zeros

sampleN<-cbind("filename"=sampleN[,1],"hawanui"=rowSums(sampleN[,2:6]),"mauinui"=rowSums(sampleN[7:9]),"oahu"=sampleN[,10],"kuainui"=rowSums(sampleN[,11:12]),sampleN[,13:31]) #lump big islands sites, Maui Nui = maui+lanaii+molokai, Kauai Nui = kuai+niihau

drops<-c("twin","broo","roga","rait","nort","pion","salm","pearmaze","john")
sampleN<-sampleN[,!(names(sampleN) %in% drops)]


#read in a stepping stone matrix to be edited
matrix<-as.matrix(read.delim(file="~/Dropbox/Crandall_tobo/stepping_stone_matrix.txt",header=F))


#make a matrix denoting which populations will be removed for each species and which will be treated as ghosts
#ghost pops where n=0 in MHI
#ghost pops where 0<n<5 in NWHI
#remove pops where n=0 in NWHI

ghost_test_mhi<-function(x){if(x < 5) "G" else "R"}
ghost_test_nwhi<-function(x){if(x == 0) "X" else if (x < 5) "G" else "R"} #R - population exists and was sampled, G - population exists but wasn't sampled, X - population doesn't exist
popstatus1<-data.frame(matrix(sapply(as.matrix(sampleN[,2:5]),ghost_test_mhi),ncol=4),row.names=sampleN[,1],stringsAsFactors=F)
popstatus2<-data.frame(matrix(sapply(as.matrix(sampleN[,6:15]),ghost_test_nwhi),ncol=10),row.names=sampleN[,1],stringsAsFactors=F)
popstatus<-cbind(popstatus1,popstatus2)
names(popstatus)<-names(sampleN[2:15])


#########MAKE MIGRATION MATRICES###############

for(i in rownames(popstatus)){ #i<-abuvai_CB
	
	print(paste("Now making matrices for",i))
	npops<-length(which(popstatus[i,]!="X")) #the number of extant populations for this species
	print(paste(npops,"extant populations"))
	nmatrix<-matrix[1:npops,1:npops] #crop the matrix to this size
		
	gmatrix<-nmatrix
	
	ghosts<-which(popstatus[i,which(popstatus[i,]!="X")]=="G")	#select the indices of the ghost populations, accounting for populations that have been removed as non-extant.
	
	print(paste("ghost pop at",names(popstatus)[which(popstatus[i,]=="G")]))
	
	if(length(ghosts)>0){
		for(j in seq(length(ghosts))){
	
		gmatrix[ghosts[j],]<-sub("\\*","c",gmatrix[ghosts[j],])
		gmatrix[,ghosts[j]]<-sub("\\*","c",gmatrix[,ghosts[j]])
	} #sub constant parameters in for * the rows and columns of the matrix that are ghost pops, if there are any ghosts pops
	}
	write.table(gmatrix,file=file.path(path,"migration_matrices",paste(i,"stepping.stone.txt",sep="_")),sep="\t",quote=F,row.names=F,col.names=F)
	
gmatrix[gmatrix=="*"]<-"m" #make the stepping-stone.1param matrix
	
	write.table(gmatrix,file=file.path(path,"migration_matrices",paste(i,"stepping.stone.1param.txt",sep="_")),sep="\t",quote=F,row.names=F,col.names=F)
	
if(length(ghosts)>0){
		for(j in seq(length(ghosts))){
	
		gmatrix[ghosts[j],]<-sub("0","c",gmatrix[ghosts[j],])
		gmatrix[,ghosts[j]]<-sub("0","c",gmatrix[,ghosts[j]])
	} #sub constant parameters in for * the rows and columns of the matrix that are ghost pops, if there are any ghosts pops
	}

gmatrix[gmatrix=="0"]<-"m" #make the n-island matrix


	write.table(gmatrix,file=file.path(path,"migration_matrices",paste(i,"n.island.txt",sep="_")),sep="\t",quote=F,row.names=F,col.names=F)

}


####################MAKE PARMFILES - Panmixia,N.island, Stepping_stone,Stepping_stone_1param#################

for(i in rownames(modeltest)) {#i<-"abuabd_CB_coleman_newH"
	dir.create(file.path(path,"Migrate_runs",i))
	setwd(file.path(path,"Migrate_runs",i))

	#Migration Matrices - read them from the migration matrices folder based on the first two sections of i. Add to these as necessary, and stick them into the for loop below
	stepping.stone<-as.matrix(read.delim(file=file.path("~/Dropbox/Crandall_tobo/migration_matrices",paste(strsplit(i,"_")[[1]][1],strsplit(i,"_")[[1]][2],"stepping.stone.txt",sep="_")),header=F))
	
	stepping.stone.1param<-as.matrix(read.delim(file=file.path("~/Dropbox/Crandall_tobo/migration_matrices",paste(strsplit(i,"_")[[1]][1],strsplit(i,"_")[[1]][2],"stepping.stone.1param.txt",sep="_")),header=F))
	
	n.island<-as.matrix(read.delim(file=file.path("~/Dropbox/Crandall_tobo/migration_matrices",paste(strsplit(i,"_")[[1]][1],strsplit(i,"_")[[1]][2],"n.island.txt",sep="_")),header=F))
	
	panmixia<-"*"
	
	#Population Lumping - add hypotheses as necessary and stick them into the if/then statements below
	npops<-length(which(popstatus[paste(strsplit(i,"_")[[1]][1],strsplit(i,"_")[[1]][2],sep="_"),]!="X")) #the number of extant populations for this species

	all_pop_relabel<-seq(npops)
	panmixia_pop_relabel<-rep(1,npops)
	
	
	for(migration.matrix in c("panmixia","stepping.stone.1param","stepping.stone","n.island")){#migration.matrix<-"stepping.stone"
	
		dir.create(file.path(path,"Migrate_runs",i,migration.matrix))
		setwd(file.path(path,"Migrate_runs",i,migration.matrix))
		
#####IF/THEN settings for writing various versions####
	
	# if the model includes a gamma parameter, then discretize the gamma shape parameter 
	model<-modeltest$bestmodel[which(rownames(modeltest)==i)]
	if(grepl("+G",model)){
		numrates<-4
		discretegammas<-c(modeltest$r1[which(rownames(modeltest)==i)],modeltest$r2[which(rownames(modeltest)==i)],modeltest$r3[which(rownames(modeltest)==i)],modeltest$r4[which(rownames(modeltest)==i)])
			
		discretegammaprobs<-c(modeltest$p1[which(rownames(modeltest)==i)],modeltest$p2[which(rownames(modeltest)==i)],modeltest$p3[which(rownames(modeltest)==i)],modeltest$p4[which(rownames(modeltest)==i)])

		} else {numrates<-1; discretegammas<-1;discretegammaprobs<-1}

	if(migration.matrix=="panmixia"){pop_relabel<-panmixia_pop_relabel} else {pop_relabel<-all_pop_relabel} 
	
	if(migration.matrix=="n.island"){increment<-"500";chainlength<-"10000"} else {increment<-"500"; chainlength<-"50000"}

#####Write the parmfiles################
	parmfile<-file("parmfile",open="w")
	
	#file header, plus some parameters that don't change for mtDNA
	cat("################################################################################",
"# Parmfile for Migrate 3.6.4 [do not remove these first TWO lines]",sep="\n",file=parmfile)
	
	#settings that will never change for sequence data
	cat("menu=NO", "nmlength=10","datatype=SequenceData",sep="\n",file=parmfile)
	
	#sequence model parameters for ttratio and base frequencies
	cat(paste("ttratio=",modeltest$titv[which(rownames(modeltest)==i)]),paste("freqs-from-data=NO:",modeltest$fA[which(rownames(modeltest)==i)],",",modeltest$fC[which(rownames(modeltest)==i)],",",modeltest$fG[which(rownames(modeltest)==i)],",",modeltest$fT[which(rownames(modeltest)==i)],sep=""),"seqerror-rate=0.0","categories=1",sep="\n",file=parmfile)
	
	#parameters for site rate variation (calculated above)
	cat(paste("rates=",numrates,":",paste(discretegammas,collapse=" "),sep=""),paste("prob-rates=",numrates,":",paste(discretegammaprobs,collapse=" "),sep=""),sep="\n",file=parmfile)
	
	#parameters that don't need to change among species (will mess with pop-relabel later)
	cat("autocorrelation=NO","weights=NO","interleaved=NO","fast-likelihood=NO","inheritance-scalars={1}",sep="\n",file=parmfile)
	
	#population re-labeling
	cat(paste("population-relabel={",paste(pop_relabel,collapse=" "),"}",sep=""),sep="\n",file=parmfile)
	
	#file names etc.
	cat("usertree=RANDOMTREE",paste("infile=../../../Migrate_datafiles_n14/",i,".mig",sep=""),"random-seed=AUTO",paste("title=",i,sep=""),"progress=YES","logfile=NO","print-data=NO",paste("outfile=",i,"_outfile.txt",sep=""),paste("pdf-outfile=",i,"_outfile.pdf",sep=""),sep="\n",file=parmfile)
	
	#more settings that don't need to change among species
	cat("use-M=YES","plot=NO","mathfile=mathfile","profile=ALL:QUICK","print-tree=NONE","write-summary=NO","aic-modeltest=NO","mig-histogram=NO","skyline=NO",sep="\n",file=parmfile)
	
	#theta start values. Can be a matrix of values input above.
	cat(paste("theta=own:{",paste(start.theta,sep=",",collapse=","),"}",sep=""),sep="\n",file=parmfile)
	
	#migration start values. Can be a matrix of values input above.
	cat(paste("migration=own:{",paste(start.migration,sep=",",collapse=","),"}",sep=""),sep="\n",file=parmfile)
	
	#more settings that don't need to change for species with msat only data
	cat("mutation=CONSTANT","fst-type=THETA",sep="\n",file=parmfile)
	
	#the migration matrix - input above
	cat(paste("custom-migration={",paste(get(migration.matrix),sep="",collapse=""),"}",sep=""),sep="\n",file=parmfile)
	
	#more settings that won't change among species
	cat("geo=NO","bayes-update=YES","bayes-updatefreq=0.500000","bayes-posteriorbins=500 500","bayes-posteriormaxtype=ALL","bayes-file=YES:bayesfile","bayes-allfile=YES:1:bayesallfile","bayes-proposals= THETA METROPOLIS Sampler","bayes-proposals= MIG SLICE Sampler",sep="\n",file=parmfile)
	
	#THETA PRIORS
	cat(paste("bayes-priors= THETA WEXPPRIOR:",paste(theta.prior,collapse=" "),sep=""),sep="\n",file=parmfile)
	
	#Migration PRIORS
	cat(paste("bayes-priors= MIG WEXPPRIOR:",paste(mig.prior,collapse=" "),sep=""),sep="\n",file=parmfile)
	
	#Chain Control - may move these settings up to "settings" later
	cat("long-chains=1",paste("long-inc=",increment,sep=""),paste("long-sample=",chainlength,sep=""),"burn-in=0",sep="\n",file=parmfile)
	
	#Heating
	cat("heating=YES:1:{1,1.5,3,100000}","heated-swap=YES",sep="\n",file=parmfile)

	#More settings having to do with convergence that won't change much, other than replicate, which should maybe be moved up
	cat("moving-steps=NO","long-chain-epsilon=INFINITY","gelman-convergence=No",replicate,"resistance=0.000100","end",sep="\n",file=parmfile)
	
	close(parmfile)
	}
}








####################MAKE PARMFILES - Regional Structure#################
#first we need to build the popstatus matrix a little differently, replacing ghost pops from the nwhi with regular pops because ghosts will screw things up down the line when they are lumped.
ghost_test_mhi<-function(x){if(x < 5) "G" else "R"}
ghost_test_nwhi<-function(x){if(x == 0) "X" else "R"} #R - population exists and was sampled, G - 
popstatus1<-data.frame(matrix(sapply(as.matrix(sampleN[,2:5]),ghost_test_mhi),ncol=4),row.names=sampleN[,1],stringsAsFactors=F)
popstatus2<-data.frame(matrix(sapply(as.matrix(sampleN[,6:15]),ghost_test_nwhi),ncol=10),row.names=sampleN[,1],stringsAsFactors=F)
popstatus<-cbind(popstatus1,popstatus2)
names(popstatus)<-names(sampleN[2:15])

#remove troublesome species
modeltest<-modeltest[-c(6,7,9,14,15),] #acapla, anon2,anon4, carign, carmel(still need to make empirical modle for this one)
############################################


for(i in rownames(modeltest)) {#i<-"abuabd_CB_coleman_newH"
	dir.create(file.path(path,"Migrate_runs_regional",i))
	setwd(file.path(path,"Migrate_runs_regional",i))
	
		
	#Population Lumping - add hypotheses as necessary and stick them into the if/then statements below
	structure<-read.table(file="~/Dropbox/Crandall_tobo/regional_structure.txt",header=T,sep="\t",row.names=1)
	
	extant<-which(popstatus[paste(strsplit(i,"_")[[1]][1],strsplit(i,"_")[[1]][2],sep="_"),]!="X") #identifiers of pops that are extant
	npops<-length(extant) #the number of extant populations for this species

	migration.matrix<-"stepping.stone"
	
	for(hypothesis in c("hi_lo","1current","2currents","structure_panmixia","structure_currents","empirical")){#"
	
		dir.create(file.path(path,"Migrate_runs_regional",i,hypothesis))
		setwd(file.path(path,"Migrate_runs_regional",i,hypothesis))

	#Migration Matrices - just the stepping-stone matrix will do here
	stepping.stone<-as.matrix(read.delim(file=file.path("~/Dropbox/Crandall_tobo/migration_matrices",paste(strsplit(i,"_")[[1]][1],strsplit(i,"_")[[1]][2],"stepping.stone.txt",sep="_")),header=F))
	
	#again, we need to remove constant parameters for ghost populations from the nwhi because they will screw with migration matrices for lumped populations
    if(npops>4){	
	stepping.stone[5:npops,][which(stepping.stone[5:npops,]=="c")]<-"*"
	stepping.stone[,5:npops][which(stepping.stone[,5:npops]=="c")]<-"*"
	}
#####IF/THEN settings for writing various versions####
	
	# if the model includes a gamma parameter, then discretize the gamma shape parameter 
	model<-modeltest$bestmodel[which(rownames(modeltest)==i)]
	if(grepl("+G",model)){
		numrates<-4
		discretegammas<-c(modeltest$r1[which(rownames(modeltest)==i)],modeltest$r2[which(rownames(modeltest)==i)],modeltest$r3[which(rownames(modeltest)==i)],modeltest$r4[which(rownames(modeltest)==i)])
			
		discretegammaprobs<-c(modeltest$p1[which(rownames(modeltest)==i)],modeltest$p2[which(rownames(modeltest)==i)],modeltest$p3[which(rownames(modeltest)==i)],modeltest$p4[which(rownames(modeltest)==i)])

		} else {numrates<-1; discretegammas<-1;discretegammaprobs<-1}

	#if(migration.matrix=="panmixia"){pop_relabel<-panmixia_pop_relabel} else {pop_relabel<-all_pop_relabel} 
	
	#if(migration.matrix=="n.island"){increment<-"500";chainlength<-"10000"} else {increment<-"500"; chainlength<-"50000"}
	
	clusters<-structure[hypothesis,extant]
	stepping.stone<-stepping.stone[1:max(clusters),1:max(clusters)]
#####Write the parmfiles################
	parmfile<-file("parmfile",open="w")
	
	#file header, plus some parameters that don't change for mtDNA
	cat("################################################################################",
"# Parmfile for Migrate 3.6.4 [do not remove these first TWO lines]",sep="\n",file=parmfile)
	
	#settings that will never change for sequence data
	cat("menu=NO", "nmlength=10","datatype=SequenceData",sep="\n",file=parmfile)
	
	#sequence model parameters for ttratio and base frequencies
	cat(paste("ttratio=",modeltest$titv[which(rownames(modeltest)==i)]),paste("freqs-from-data=NO:",modeltest$fA[which(rownames(modeltest)==i)],",",modeltest$fC[which(rownames(modeltest)==i)],",",modeltest$fG[which(rownames(modeltest)==i)],",",modeltest$fT[which(rownames(modeltest)==i)],sep=""),"seqerror-rate=0.0","categories=1",sep="\n",file=parmfile)
	
	#parameters for site rate variation (calculated above)
	cat(paste("rates=",numrates,":",paste(discretegammas,collapse=" "),sep=""),paste("prob-rates=",numrates,":",paste(discretegammaprobs,collapse=" "),sep=""),sep="\n",file=parmfile)
	
	#parameters that don't need to change among species (will mess with pop-relabel later)
	cat("autocorrelation=NO","weights=NO","interleaved=NO","fast-likelihood=NO","inheritance-scalars={1}",sep="\n",file=parmfile)
	
	#population re-labeling
	cat(paste("population-relabel={",paste(clusters,collapse=" "),"}",sep=""),sep="\n",file=parmfile)
	
	#file names etc.
	cat("usertree=RANDOMTREE",paste("infile=../../../Migrate_datafiles_n14/",i,".mig",sep=""),"random-seed=AUTO",paste("title=",i,sep=""),"progress=YES","logfile=NO","print-data=NO",paste("outfile=",i,"_outfile.txt",sep=""),paste("pdf-outfile=",i,"_outfile.pdf",sep=""),sep="\n",file=parmfile)
	
	#more settings that don't need to change among species
	cat("use-M=YES","plot=NO","mathfile=mathfile","profile=ALL:QUICK","print-tree=NONE","write-summary=NO","aic-modeltest=NO","mig-histogram=NO","skyline=NO",sep="\n",file=parmfile)
	
	#theta start values. Can be a matrix of values input above.
	cat(paste("theta=own:{",paste(start.theta,sep=",",collapse=","),"}",sep=""),sep="\n",file=parmfile)
	
	#migration start values. Can be a matrix of values input above.
	cat(paste("migration=own:{",paste(start.migration,sep=",",collapse=","),"}",sep=""),sep="\n",file=parmfile)
	
	#more settings that don't need to change for species with msat only data
	cat("mutation=CONSTANT","fst-type=THETA",sep="\n",file=parmfile)
	
	#the migration matrix - input above
	cat(paste("custom-migration={",paste(get(migration.matrix),sep="",collapse=""),"}",sep=""),sep="\n",file=parmfile)
	
	#more settings that won't change among species
	cat("geo=NO","bayes-update=YES","bayes-updatefreq=0.500000","bayes-posteriorbins=500 500","bayes-posteriormaxtype=ALL","bayes-file=YES:bayesfile","bayes-allfile=YES:1:bayesallfile","bayes-proposals= THETA METROPOLIS Sampler","bayes-proposals= MIG SLICE Sampler",sep="\n",file=parmfile)
	
	#THETA PRIORS
	cat(paste("bayes-priors= THETA WEXPPRIOR:",paste(theta.prior,collapse=" "),sep=""),sep="\n",file=parmfile)
	
	#Migration PRIORS
	cat(paste("bayes-priors= MIG WEXPPRIOR:",paste(mig.prior,collapse=" "),sep=""),sep="\n",file=parmfile)
	
	#Chain Control - may move these settings up to "settings" later
	cat("long-chains=1",paste("long-inc=",increment,sep=""),paste("long-sample=",chainlength,sep=""),"burn-in=0",sep="\n",file=parmfile)
	
	#Heating
	cat("heating=YES:1:{1,1.5,3,100000}","heated-swap=YES",sep="\n",file=parmfile)

	#More settings having to do with convergence that won't change much, other than replicate, which should maybe be moved up
	cat("moving-steps=NO","long-chain-epsilon=INFINITY","gelman-convergence=No",replicate,"resistance=0.000100","end",sep="\n",file=parmfile)
	
	close(parmfile)
	}
}

	
	
	
	



	
	
