library(strataG)

setwd("~/Dropbox/Crandall_tobo/Final_FASTA/")

summstats<-data.frame()

for(fastafile in list.files("~/Dropbox/Crandall_tobo/Final_FASTA/")){
  
  print(fastafile)
  
  f<-read.fasta(fastafile)
  pops<-sub(pattern=".+population:([A-Za-z]+).+",replacement="\\1",x=labels(f))
  g<-DNAbin2genind(f, pop=pops)
  h<-genind2gtypes(g)
  phist<-statPhist(h, nrep=1000)
  phist_value<-phist$result[1]
  phist_pvalue<-phist$result[2]
  tajd<-tajimasD(f)
  tajd_value<-tajd[1]
  tajd_pvalue<-tajd[2]
  k<-labelHaplotypes(f, prefix="Hap")
  haps<-length(labels(k$hap.seqs))
  hapdiv<-swfscMisc::diversity(k$haps)
  
  summstats<-rbind(summstats,data.frame(fastafile, phist_value, phist_pvalue, tajd_value, tajd_pvalue,haps,hapdiv))
}