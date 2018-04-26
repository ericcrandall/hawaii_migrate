# analyze bayesallfiles
# using coda and ggmcmc

# Load packages
library(coda)
library(ggplot2)
library(ggmcmc)
library(gridExtra)
library(reshape2)


#a function for creating Nm vectors out of m and Theta vectors.
migrants.per.gen<-function(x){
  #x<-x[[1]]
  m<-names(x)[which(grepl("M_",names(x)))] #names of m columns
  #theta<-names(x)[which(grepl("Theta_",names(x)))] #names of theta columns
  for(n in m){
    t<-paste("Theta",strsplit(n,split="_")[[1]][3],sep="_")
    x[,paste("Nm",strsplit(n,split="_")[[1]][2],strsplit(n,split="_")[[1]][3],sep="_")]<-  	x[,which(names(x)==n)]*x[,which(names(x)==t)] #this hairy little statement makes a new column named "Nm_X_Y" and then fills it by multiplying the M_X_Y column by the Theta_Y column	
  }
  return(x)
}

#a function to plot regressions from lm objects (from https://susanejohnston.wordpress.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/)
ggplot.lm <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}


# read in the best model table and other info
bestmodels<-read.csv("/Users/eric/Dropbox/Crandall_tobo/scripts/Final_BestModeltable_Using_Means_of_3reps.csv",stringsAsFactors=F)
#read in the population status table (R = regular, G = ghost (unsampled), X = nonextant)
popstatus<-read.table("/Users/eric/Dropbox/Crandall_tobo/scripts/popstatus.txt",header=T,sep="\t", stringsAsFactors = F)
#read in a stepping-stone graph
steporder<-read.csv("/Users/eric/Dropbox/Crandall_tobo/scripts/step_edges_migrate.csv",header=F)
# set the working directory
hawa_area<-read.csv(file="/Users/eric/Dropbox/Crandall_tobo/Hawaii_shallow_area.csv", stringsAsFactors = F)


wd<-"/Users/eric/Datasets/ToBo_Migrate_runs"
finaloutput<-"/Users/eric/Datasets/ToBo_Project/finalregressions"
setwd(wd)


#initialize an empty list
allstatstable<-list()

for(f in bestmodels$dataset){ #f<-"panpen_CI_Iacchei"
  
  #get species name and add underscores for spaces
  speciesname<-bestmodels[which(bestmodels$dataset==f),1]
  speciesname2<-gsub("\\s","_",speciesname,perl=T)
  
  #make an output directory
  outputdir<-file.path(finaloutput,speciesname2)
  dir.create(outputdir)
  
  #get the best model
  model<-bestmodels[which(bestmodels$dataset==f),3]
  print(paste("Now starting",f,";",speciesname,"; best model = ",model))
  
  data.list<-list() #initialize an empty data list
  
  for(r in 6:8){  #loop through the 3 replicate rounds (6,7,8)
    wd1<-file.path(wd,paste("round",r,sep=""),f,model)
    setwd(wd1)  #move into that species' and model directory
    
    print(paste("loading and processing bayesallfile for round", r))
    data<-read.table("bayesallfile", header=T) #this may take a minute or two
  
    # Split the whole list into the individual replicates, so that you will
    # have a list of data frames
    data.list.1<-split(data,data$Replicate)
    
    # Subset the parameters of interest, either by rows or by columns to remove burnin and uninteresting columns
    data.list.1<-lapply(data.list.1,subset,subset=Steps>(0.4*max(data.list.1[[1]]$Steps)))# by rows. Removing burnin. First 40% First 20,000 samples out of 50,000 or 10 million steps out of 25 million
    #data.list.1<-lapply(data.list.1,subset,select=c(4:(length(data.list.1[[1]])-6))) #by columns.
    
    # calculate Nm for each row
    data.list.2<-lapply(data.list.1,migrants.per.gen)
    
    # cat it onto the main list
    data.list<-c(data.list,data.list.2)
  }


  
  
  #convert each dataframe to an mcmc object and then convert the whole thing to an mcmc list
  data.list.mcmc<-mcmc.list(lapply(data.list,mcmc))
  #condense everything into a single mcmc replicate for the purposes of HPDinterval
  data.list.allinone<-mcmc(data=do.call("rbind",data.list)) 
  
  #calculate statistics
  print("calculating statistics")
  summ<-summary(data.list.mcmc)
  ess<-effectiveSize(data.list.mcmc)
  gelman<-gelman.diag(data.list.mcmc,multivariate=F)
  HPD<-HPDinterval(data.list.allinone)
  
  #cat the stats, man
  allstats<-cbind(summ$statistics[,1:2],HPD,ess,gelman$psrf)
  allstatstable[[f]]<-allstats
  
  #write the stats into the directory for future reference
  write.csv(allstats,file.path(outputdir,paste(speciesname2,model,"codastats.csv")))
  
  if(model=="panmixia"){next}
  
  #create plots of various diagnostics for visual evaluation, 
  print("plotting MCMC")
  
  data.list.plot<-lapply(data.list,subset,select=c(which(grepl("Theta_",names(data.list[[1]])) | grepl("M_",names(data.list[[1]])) | grepl("Nm_",names(data.list[[1]])) |grepl("lnPost",names(data.list[[1]])))))
  #log transform them, since they come from an exponential prior
  data.list.plot<-lapply(data.list.plot,log)
  data.list.plot<-mcmc.list(lapply(data.list.plot,mcmc))

  data.list.gg<-ggs(data.list.plot,burnin=F,description=paste(f,model,sep="_"))
  ggmcmc(data.list.gg,plot=c("ggs_traceplot","ggs_density","ggs_Rhat"),simplify_traceplot=0.25,file=paste("/Users/eric/Datasets/ToBo_Project/final_ggMCMC/",paste(speciesname2,model,sep="_"),".pdf",sep=""))
  
  
  
  cat("Plotting Nm Violin Plots") 
  #isolate Nm values for violin plots
  nmlist<-lapply(data.list,subset,select=c(which(grepl("Nm",names(data.list[[1]])))))
  nmlistallinone<-data.frame(data=do.call("rbind",nmlist))
  #get the first two parts of f into a string
  sl<-strsplit(f,"_")[[1]][1:2]
  sl<-paste(sl[1],sl[2],sep="_")
  #extract the population names that were sampled for this species
  pops<-popstatus[sl,]
  popnames<-names(pops[which(pops=="R")])
  #make a vector of the combos of names for Nm and m
  stepnames<-paste(popnames[steporder[0:(2*length(popnames)-2),1]],popnames[steporder[0:(2*length(popnames)-2),2]],sep=" to ")
  
  
  #melt nmlistallinone to a long data frame
  longNm <- melt(data = nmlistallinone,
                 measure.vars = names(nmlistallinone),
                 variable.name = "Parameter",
                 value.name = "Nm"
  )
  
  #make a vector to color the violins by direction of gene flow
  paramlength<-length(longNm[,1])/length(levels(longNm$Parameter)) #get the number of observations per parameter
  direction<-rep_len(c(rep_len("south",paramlength),rep_len("north",paramlength)),length(longNm[,2]))
  
  violin<-ggplot(longNm, (aes(x=Parameter, y=Nm, fill=direction))) + geom_violin(draw_quantiles=c(0.025,0.5,0.975))  + 
    scale_y_log10(breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000,10000,100000), labels=c("0.00001","0.0001","0.001","0.01","0.1","1","10","100","1,000","10,000","100,000"),limits=c(0.001,100000)) +
    scale_x_discrete(labels=stepnames) + scale_fill_manual(values=c("grey","white")) +
    labs(title=f, x="Parameter", y="Effective Female Migrants per Generation") + coord_flip() + 
    ggtitle(speciesname) + guides(fill=F)
  
  ggsave(file=file.path(outputdir,paste(speciesname2,"_",model,"_Nm_violin.pdf",sep="")),plot=violin,device="pdf",width=11,height=8.5,units="in")
  
  #pdf(file=file.path(outputdir,paste(speciesname2,"_",bestmodel,"_Nm_violin.pdf",sep="")), width=8.5, height=11, paper="letter") 
  #print(violin)
  #dev.off()
  
  cat("Plotting m/mu Violin Plots") 
  #isolate m values for violin plots
  mlist<-lapply(data.list,subset,select=c(which(grepl("M",names(data.list[[1]])))))
  mlistallinone<-data.frame(data=do.call("rbind",mlist))
  
  #melt mlistallinone to a long data frame
  longm <- melt(data = mlistallinone,
                measure.vars = names(mlistallinone),
                variable.name = "Parameter",
                value.name = "m"
  )
  
  #make a vector to color the violins by direction of gene flow
  paramlength<-length(longm[,1])/length(levels(longm$Parameter)) #get the number of observations per parameter
  direction<-rep_len(c(rep_len("south",paramlength),rep_len("north",paramlength)),length(longm[,2]))
  
  violinm<-ggplot(longm, (aes(x=Parameter, y=m, fill=direction))) + geom_violin(draw_quantiles=c(0.025,0.5,0.975))  + 
    scale_y_log10(breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000,10000,100000,1000000), labels=c("0.00001","0.0001","0.001","0.01","0.1","1","10","100","1,000","10,000","100,000","1,000,000"),limits=c(0.0001,1000000)) +
    scale_x_discrete(labels=stepnames) + scale_fill_manual(values=c("grey","white")) +
    labs(title=f, x="Parameter", y="proportion of migrants scaled by mutation rate") + coord_flip() + 
    ggtitle(speciesname) + guides(fill=F)
  
  ggsave(file=file.path(outputdir,paste(speciesname2,"_",model,"_m_violin.pdf",sep="")),plot=violinm,device="pdf",width=11,height=8.5,units="in")
  
  
  
  cat("Plotting Theta Violin Plots") 
  #isolate theta values for violin plots
  tlist<-lapply(data.list,subset,select=c(which(grepl("Theta",names(data.list[[1]])))))
  tlistallinone<-data.frame(data=do.call("rbind",tlist))
  
  #melt mlistallinone to a long data frame
  longt <- melt(data = tlistallinone,
                measure.vars = names(tlistallinone),
                variable.name = "Parameter",
                value.name = "Theta"
  )
  
  violint<-ggplot(longt, (aes(x=Parameter, y=Theta))) + geom_violin(draw_quantiles=c(0.025,0.5,0.975))  + 
    scale_y_log10(breaks=c(0.00001,0.0001,0.001,0.01,0.1), labels=c("0.00001","0.0001","0.001","0.01","0.1"),limits=c(0.00001,0.1)) +
    scale_x_discrete(labels=popnames) + scale_fill_manual(values=c("grey","white")) +
    labs(title=f, x="Parameter", y="Theta") + coord_flip() + 
    ggtitle(speciesname) + guides(fill=F)
  
  ggsave(file=file.path(outputdir,paste(speciesname2,"_",model,"_theta_violin.pdf",sep="")),plot=violint,device="pdf",width=11,height=8.5,units="in")
  
  cat("Plotting Stats as barplots")
  
  speciesNm<-as.data.frame(allstats[which(grepl("Nm",rownames(allstats))),])
  speciesTheta<-as.data.frame(allstats[which(grepl("Theta",rownames(allstats))),])
  
  if(length(speciesNm[1,])==1){cat("1 Nm parameter; skipping",speciesname, "\n"); next}
  if(length(speciesNm[,1])==0){cat("Panmixia; skipping",speciesname, "\n");next}
  
  pdf(file=paste("/Users/eric/Datasets/ToBo_Project/final_ggmcmc/",paste(speciesname2,model,"bars",sep="_"),".pdf",sep=""),width=8.5, height=11, paper="letter")
  
  
  nmstats<-ggplot(speciesNm) + 
    geom_crossbar( aes(x = rownames(speciesNm), y = Mean, ymin = lower, ymax=upper, color = speciesNm[,6] < 1.2, fill= ess > 200)) + 
    labs(title=speciesname, x="Parameter", y="Effective Female Migrants per Generation")  + 
    scale_y_log10(breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000,10000,100000), labels=c("0.00001","0.0001","0.001","0.01","0.1","1","10","100","1,000","10,000","100,000")) +
    scale_fill_manual(values=c("black","white"))
  
  thetastats<-ggplot(speciesTheta) + 
    geom_crossbar( aes(x = rownames(speciesTheta), y = Mean, ymin = lower, ymax=upper, color = speciesTheta[,6] < 1.2, fill= ess > 200)) + 
    labs(title=speciesname, x="Parameter", y="Effective Females x mutation rate")  + 
    scale_y_log10(breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,10), labels=c("0.00001","0.0001","0.001","0.01","0.1","1","10")) +
    scale_fill_manual(values=c("black","white"))
  
  grid.arrange(nmstats,thetastats,nrow=2)
  
  dev.off()
  

  cat("\n","done","next",sep="\n\n")
}

save(allstatstable,file=file.path(outputdir,"ToBo_bestmodel_stats_v.Rdata"))




# "Bayesian" regressions for questions about habitat area etc (new loop)

setwd(wd)
theta_area_plots<-list()
theta_area_pvalues<-list()

theta_census_plots<-list()
theta_census_pvalues<-list()

for(f in bestmodels$dataset){ #f<-"panpen_CI_Iacchei" f<-"abuvai_CB_coleman_newH"
  
  finaloutput<-"/Users/eric/Datasets/ToBo_Project/final_regressions/"
  
  #get species name and add underscores for spaces
  speciesname<-bestmodels[which(bestmodels$dataset==f),1]
  speciesname2<-gsub("\\s",".",speciesname,perl=T)
  
  #get the best model
  model<-bestmodels[which(bestmodels$dataset==f),3]
  print(paste("Now starting",f,";",speciesname,"; best model = ",model))
  if(model!="stepping.stone"){print("Not Stepping-Stone");next}
  
  data.list<-list() #initialize an empty data list
  
  for(r in 6:8){  #loop through the 3 replicate rounds (6,7,8)
    wd1<-file.path(wd,paste("round",r,sep=""),f,model)
    setwd(wd1)  #move into that species' and model directory
    
    print(paste("loading and processing bayesallfile for round", r))
    data<-read.table("bayesallfile", header=T) #this may take a minute or two
    
    # Split the whole list into the individual replicates, so that you will
    # have a list of data frames
    data.list.1<-split(data,data$Replicate)
    
    # Subset the parameters of interest, either by rows or by columns to remove burnin and uninteresting columns
    data.list.1<-lapply(data.list.1,subset,subset=Steps>(0.4*max(data.list.1[[1]]$Steps)))# by rows. Removing burnin. First 40% First 20,000 samples out of 50,000 or 10 million steps out of 25 million
    #data.list.1<-lapply(data.list.1,subset,select=c(4:(length(data.list.1[[1]])-6))) #by columns.
    
    # calculate Nm for each row
    data.list.2<-lapply(data.list.1,migrants.per.gen)
    
    # cat it onto the main list
    data.list<-c(data.list,data.list.2)
  }
  
  
  
  
  #convert each dataframe to an mcmc object and then convert the whole thing to an mcmc list
  #data.list.mcmc<-mcmc.list(lapply(data.list,mcmc))
  #condense everything into a single mcmc replicate for the purposes of HPDinterval
  dataframe<-do.call("rbind",data.list) # change the concatenated list to a data frame
  
  
  #get population names and migration vector names
  
  #get the first two parts of f into a string
  sl<-strsplit(f,"_")[[1]][1:2]
  sl<-paste(sl[1],sl[2],sep="_")
  #extract the population names that were sampled for this species
  pops<-popstatus[sl,]
  popnames<-names(pops[which(pops=="R")])
  #make a vector of the combos of names for Nm and m
  stepnames<-paste(popnames[steporder[0:(2*length(popnames)-2),1]],popnames[steporder[0:(2*length(popnames)-2),2]],sep=" to ")
  
  for(i in length(popnames):1){#count backwards to prevent substituting the 1 in '10'
    names(dataframe)<-gsub(paste("Theta",i,sep="_"),popnames[i],names(dataframe))
  }



## Theta~Area regressions
print("Bayesian regressions - Theta ~ Area")

#pick out the thetas and take a random sample of 10000 from them
randtheta <- sample(length(rownames(dataframe)), size=10000, replace=F)
subsample<-dataframe[randtheta,popnames]
rownames(subsample)<-NULL

areas<-hawa_area$Hardbottom_30m[which(hawa_area$island %in% popnames)] #lookup the areas of the islands that are in the dataset - the sub command strips "theta_" off the name


pearsonr.ta<-NULL
rsquared.ta<-NULL
error.ta<-NULL
slope.ta<-NULL
intercept.ta<-NULL
logl.ta<-NULL

subsample<-as.matrix(subsample)
#plot(log(subsample[1155,]) ~ areas)
for(n in 1:length(subsample[,1])){
  theta_area<-lm(formula = log10(subsample[n,]) ~ log10(areas) )
  
  #  abline(theta_area,col="grey",lty="dotted",lwd="0.2")
  
  rsquared.ta<-c(rsquared.ta,summary(theta_area)$r.squared)
  error.ta<-c(error.ta,summary(theta_area)$sigma)
  intercept.ta<-c(intercept.ta,as.numeric(coef(theta_area)[1]))
  slope.ta<-c(slope.ta,as.numeric(coef(theta_area)[2]))
  logl.ta<-c(logl.ta,logLik(theta_area))
  
  
}

pvalue<-length(which(slope.ta < 0))/length(slope.ta)

theta_area_pvalues[[speciesname2]]<-pvalue

theta_area_plots[[speciesname2]]<-ggplot(data.frame(slope.ta),aes(slope.ta)) + geom_histogram() + geom_vline(xintercept=0,color="Red") +
  labs(title = paste("Slopes for Theta ~ Area for:",speciesname,"p =",pvalue))


## Theta~CensusSize regressions
print("Bayesian regressions - Theta ~ CensusSize")

#pick out the thetas and take a random sample of 10000 from them
randtheta <- sample(length(rownames(dataframe)), size=10000, replace=F)
subsample<-dataframe[randtheta,popnames]
rownames(subsample)<-NULL

areas<-hawa_area$Hardbottom_30m[which(hawa_area$island %in% popnames)] #lookup the areas of the islands that are in the dataset - the sub command strips "theta_" off the name
density<-hawa_area[which(hawa_area$island %in% popnames),which(names(hawa_area)==speciesname2)]
if(length(density)==0){print("Not a fish, so can't do census size");next}
density[which(density==0)]<-0.000001
census<-areas*density*10000



pearsonr.ta<-NULL
rsquared.ta<-NULL
error.ta<-NULL
slope.ta<-NULL
intercept.ta<-NULL
logl.ta<-NULL

subsample<-as.matrix(subsample)
#plot(log(subsample[1155,]) ~ areas)
for(n in 1:length(subsample[,1])){
  theta_census<-lm(formula = log10(subsample[n,]) ~ log10(census) )
  
  #  abline(theta_area,col="grey",lty="dotted",lwd="0.2")
  
  rsquared.ta<-c(rsquared.ta,summary(theta_census)$r.squared)
  error.ta<-c(error.ta,summary(theta_census)$sigma)
  intercept.ta<-c(intercept.ta,as.numeric(coef(theta_census)[1]))
  slope.ta<-c(slope.ta,as.numeric(coef(theta_census)[2]))
  logl.ta<-c(logl.ta,logLik(theta_census))
  
  
}


pvalue<-length(which(slope.ta < 0))/length(slope.ta)

theta_census_pvalues[[speciesname2]]<-pvalue

theta_census_plots[[speciesname2]]<-ggplot(data.frame(slope.ta),aes(slope.ta)) + geom_histogram() + geom_vline(xintercept=0,color="Red") +
  labs(title = paste("Slopes for Theta ~ Census Size for:",speciesname,"p =",pvalue))
                                                                                                                                    
}

setwd(finaloutput)
ggsave("LogTheta_LogHardbottomArea_Plots.pdf",marrangeGrob(grobs = theta_area_plots, nrow=2, ncol=2))

ggsave("LogTheta_LogCensus_Plots.pdf",marrangeGrob(grobs = theta_census_plots, nrow=2, ncol=2))



# Read in the allstats and plot max and min values for each parameter for each species as bar charts
load("~/Google Drive/editing-reviewing/Hawaii_Migrate/ToBo_Project/final_violins/ToBo_bestmodel_stats_v.Rdata")

goodones<-bestmodels$dataset[which(bestmodels$ttest_pvalue<0.05 & bestmodels$bestmodel!="panmixia" | bestmodels$dataset=="acaoli_CB_gaither")]
goodones<-sort(goodones)

allNms<-data.frame(name=character(22),lower=numeric(22),upper=numeric(22),row.names = goodones,stringsAsFactors = F)
allThetas<-data.frame(name=character(22),lower=numeric(22),upper=numeric(22),row.names = goodones,stringsAsFactors = F)
allMs<-data.frame(name=character(22),lower=numeric(22),upper=numeric(22),row.names = goodones,stringsAsFactors = F)


for(f in goodones){ #f<-"panpen_CI_Iacchei" f<-"abuvai_CB_coleman_newH"
  print(f)
  
  #get species name 
  speciesname<-bestmodels[which(bestmodels$dataset==f),1]
  
  #pull out the stats
  speciesNm<-allstatstable[[f]][which(grepl("Nm",rownames(allstatstable[[f]]))),]
  speciesTheta<-allstatstable[[f]][which(grepl("Theta",rownames(allstatstable[[f]]))),]
  speciesM<-allstatstable[[f]][which(grepl("M",rownames(allstatstable[[f]]))),]
  
  #deal with species that only have one value
  if(is.null(dim(speciesNm))){
    allNms[f,]<-c(speciesname,min(speciesNm[3]),max(speciesNm[4]))
    allThetas[f,]<-c(speciesname,min(speciesTheta[3]),max(speciesTheta[4]))
    allMs[f,]<-c(speciesname,min(speciesM[3]),max(speciesM[4]))
  }
  else{
  allNms[f,]<-c(speciesname,min(speciesNm[,3]),max(speciesNm[,4]))
  allThetas[f,]<-c(speciesname,min(speciesTheta[,3]),max(speciesTheta[,4]))
  allMs[f,]<-c(speciesname,min(speciesM[,3]),max(speciesM[,4]))
  }
   
} 

allNms<-as.data.frame(allNms)
allNms$lower<-as.numeric(allNms$lower)
allNms$upper<-as.numeric(allNms$upper)

allThetas<-as.data.frame(allThetas)
allThetas$lower<-as.numeric(allThetas$lower)
allThetas$upper<-as.numeric(allThetas$upper)

allMs<-as.data.frame(allMs)
allMs$lower<-as.numeric(allMs$lower)
allMs$upper<-as.numeric(allMs$upper)


allNms$type<-"Nm"
allThetas$type<-"Theta"
allMs$type<-"m"

all<-rbind(allNms,allThetas,allMs)

all$type<-factor(all$type,levels=c("Theta","m","Nm"))

dodge<-position_dodge(width=0.9)

bars<-ggplot(data=all) + 
  geom_linerange(data=all, aes(x = name, ymin = lower, ymax=upper, color=type), position=dodge) +
  scale_y_log10(breaks=c(0.000001,0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000,10000,100000,1000000), labels=c("0.000001","0.00001","0.0001","0.001","0.01","0.1","1","10","100","1,000","10,000","100,000","1,000,000")) +
  scale_color_manual(values=c("red","blue","purple"), labels=c(expression(Theta),expression(m/mu),expression(paste(N[e],m)))) +
  theme(axis.text.x = element_text(angle=90,hjust=1), axis.title.x = element_blank(), axis.title.y=element_blank(), legend.title=element_blank(), legend.text=element_text(size=12))
