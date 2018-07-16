migrate_harvester<-function(wd,n=3,models,multilocus=F,quiet=F){ #wd<-"/Users/eric/Datasets/simulations/ibdsim/migrate_pareto2_mtdna_results"; n=3; models<-c("panmixia","5island","5stepping.stone","5stepping.stone1","10island","10stepping.stone","10stepping.stone1")
  
  setwd(wd)
  likelists<-list() #initialize an empty list
  
  for(r in list.files(wd)){
    wd1<-file.path(wd,r)
    setwd(wd1)
    marglike<-data.frame(model=character(0),thermodynamic=numeric(0),bezier.corrected=numeric(0),harmonic.mean=numeric(0),stringsAsFactors=F) #initialize a data frame to take the values
    l=1 #initialize l
    for(m in models){ #i<-"stepping.stone"
      wd2<-file.path(wd1,m)
      if(quiet==F) {print(wd2)}
      if(!file.exists(wd2)){next}
      setwd(wd2)
      outfile<-scan(file="outfile",what="character",sep="\n",quiet=T) #scan in the outfile, separating at each newline
      
      if(multilocus==F){
      #get the result from thermodynamic integration
      thermoline<-grep("(1a)",outfile,value=T) #find the line with the thermodynamic likelihood on it
      if(length(thermoline)==0){next}
      thermoline<-strsplit(thermoline,split="=",fixed=T) #split it up
      thermo<-as.numeric(substr(thermoline[[1]][2],start=1,stop=12)) #grab the thermodynamic likelihood
      bezistart<-grep("\\(",strsplit(thermoline[[1]][2],split="")[[1]])+1
      bezier<-as.numeric(substr(thermoline[[1]][2],start=bezistart,stop=bezistart+11)) #and the bezier-corrected value
      #get the harmonic mean
      harmoline<-grep("\\(2\\) H",outfile,value=T) #find the line with harmonic mean likelihood on it
      harmoline<-strsplit(harmoline,split="=",fixed=T) #grab the harmonic mean
      harmo<-as.numeric(substr(harmoline[[1]][2],start=1,stop=12))
      marglike[l,]<-c(m,thermo,bezier,harmo) #add this as a row to the data frame
      l=l+1
      }
      
      if(multilocus==T){
        all.locus.line<-outfile[grep("\\[Scaling",outfile,value=F)-1] #find the line with all three values on it (it comes right before the line that has the scaling factor on it)
        all.locus.line<-strsplit(all.locus.line," +")
        if(length(all.locus.line)==0){next}
        thermo<-as.numeric(all.locus.line[[1]][3])
        bezier<-as.numeric(all.locus.line[[1]][4])
        harmo<-as.numeric(all.locus.line[[1]][5])
        marglike[l,]<-c(m,thermo,bezier,harmo) #add this as a row to the data frame
        l=l+1
      }
    }
    likelists[[r]]<-marglike
  }
  setwd(wd)
  return(likelists)
}


bfcalcs<-function(df,ml="bezier.corrected"){
  df$thermodynamic<-as.numeric(df$thermodynamic)
  df$bezier.corrected<-as.numeric(df$bezier.corrected)
  df$harmonic.mean<-as.numeric(df$harmonic.mean)
  mlcol<-df[,ml] 
  bmvalue<-mlcol[which.max(mlcol)]
  lbf<-2*(mlcol-bmvalue)
  choice<-rank(-mlcol)
  modelprob<-exp(lbf/2)/sum(exp(lbf/2))
  dfall<-cbind(df,lbf,choice,modelprob)
  return(dfall)
}


getmodels<-function(dfr){
  model1<-dfr$model[which(dfr$choice==1)]
  model2<-dfr$model[which(dfr$choice==2)]
  modelprob1<-dfr$modelprob[which(dfr$choice==1)]
  modelprob2<-dfr$modelprob[which(dfr$choice==2)]
  
  print(model1)
  print(model2)
  print(modelprob1)
  c(model1,modelprob1,model2,modelprob2)
}

bfcalcs_reps<-function(modeltables,models,ml_type="bezier.corrected", reps=3){
  #modeltables must be a list of bfcalcs output
  #ml_type = "bezier.corrected","thermodynamic", or "harmonic.mean"
  
  require(perm)
  output<-list()
  
  #create an empty data frame of reps*length(models) rows and cbind each list element into it.
  likes<-data.frame(matrix(nrow=reps*length(models),ncol=4))
  rowindex<-1
  for(r in 1:reps){
   likes[rowindex:(r*length(models)),]<-cbind(modeltables[[r]][,1:4])
   rowindex<-rowindex+length(models)
  }
  colnames(likes)<-c("model","thermodynamic","bezier.corrected","harmonic.mean")


  likes$model<-factor(likes$model, models) #factor to turn missing models into NA
  likes<-likes[!(is.na(likes$model)),] #remove NAs 
  likes$model<-factor(likes$model, models) #and refactor
  mean.ml<-as.vector(by(likes[[ml_type]],likes$model,mean)) #take the means for each model
  sd.ml<-as.vector(by(likes[[ml_type]],likes$model,sd))
  sem.ml<-as.vector(by(likes[[ml_type]],likes$model,sem))
    
  bmvalue<-mean.ml[which.max(mean.ml)] #get the best mean marg like
  lbf<-2*(mean.ml-bmvalue) #get the log bayes factors
  choice<-rank(-mean.ml) # make a choice column
  modelprob<-exp(lbf/2)/sum(exp(lbf/2)) # relative model probability
 
  dfall<-data.frame(models,mean.ml,sd.ml,sem.ml,lbf,choice,modelprob) #put it all together
  
  #do permutation tests
  kmeans.p<-permKS(x = likes[[ml_type]], g = likes$model, method="exact.mc")$p.value
  ttest.p<-permTS(x=likes[[ml_type]][which(likes$model==dfall$models[which(dfall$choice==1)])], 
                y=likes[[ml_type]][which(likes$model==dfall$models[which(dfall$choice==2)])], 
                alternative="greater", method="exact.mc")$p.value #non-parametric t-test
  
  output[["table"]]<-dfall
  output["perm_anova.p"]<-kmeans.p
  output["perm_t.test.p"]<-ttest.p
  return(output)
  
}


bfplot<-function(modeltables,models,ml_type="bezier.corrected", reps=3,title="Title"){

  pdf(file = "thermodynamic_marginal_likelihoods_final.pdf",width=8.5,height=3)
  
  means<-list() #initialize means
  
  for(dataset in names(modeltables)){
    
    likes<-data.frame(matrix(nrow=reps*length(models),ncol=4))
    rowindex<-1
    for(r in 1:reps){
      likes[rowindex:(r*length(models)),]<-cbind(modeltables[[r]])
      rowindex<-rowindex+length(models)
    }
    colnames(likes)<-c("model","thermodynamic","bezier.corrected","harmonic.mean")
    likes$model<-factor(likes$model, models)
    likes<-likes[!(is.na(likes$model)),]
    
    #likes<-likes[which(likes[ml_type] > max(likes$bezier.corrected)-100),]
    y.mean<-as.vector(by(likes[[ml_type]],likes$model,mean))
    y.sd<-as.vector(by(likes[[ml_type]],likes$model,sd))
    y.min<-y.mean-((y.sd/sqrt(3))*4.303)
    y.max<-y.mean+((y.sd/sqrt(3))*4.303)
    
    
    likes.mean<-data.frame(model=factor(models,models),y.mean,y.min,y.max,y.sd)
    means[[dataset]]<-likes.mean
    
    #l<-ggplot(data=likes, aes(x=model,y=bezier.corrected,colour=factor(rep), 
    #                                 shape=factor(rep), size=20 ))
    l<-ggplot(data=likes, aes(x=model,y=bezier.corrected))
    
    l<-l+geom_point(colour="blue", size=3)+
      geom_pointrange(data=likes.mean,y=y.mean,ymin=y.min,ymax=y.max, size=0.5)+
      scale_x_discrete(drop=FALSE)+
      theme(axis.text.y = element_text(size=16),legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_blank(),plot.title=element_text(size=20))+
      ggtitle(title)+ylab("Marginal Log-Likelihood")+
      coord_fixed(0.1)+ coord_flip()
    print(l)
    #  plots<-c(plots,l)
  }
  dev.off()
}

# a ggplot wrapper to plot regressions from https://susanejohnston.wordpress.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

# this one is specific to my needs for IBD plots
ggplotRegression2 <- function (fit, title) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = title) + 
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())
}


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

NeGivenMu<-function(x,mu){
  t<-names(x)[which(grepl("Theta_",names(x)))]
  for(n in t){
    x[,paste("Ne",strsplit(n,split="_")[[1]][2],sep="_")]<-x[,which(names(x)==n)]/mu #this hairy little statement makes a new column named "Ne_X" and then fills it by dividing the ThetaX column by mu
  }
  return(x)
}

#Equation2 from Waples and Gaggioti (Cockerham and Weir, 1987,1993)
expectedTheta<-function(N,m,mu,n){
  #N = Ne, m = fraction of migrants,mu = mutation rate,n = number of demes 
  expected<- 1/(1 + 4*N*mu + (4*N*m*n/(n-1)))
}

sem<-function(x){
  sd(x)/sqrt(length(x))
}

