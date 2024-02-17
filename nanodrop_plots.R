##dna isolation quality summary 

setwd("D:/jayasuriya_documents/dna_isolations/")

library(stringr)

ss1<-read.csv("qc1.csv")

ss2<-read.csv("qc2.csv")

dna_qual<-rbind(ss1,ss2)

dna_qual$Sample.Name<-paste0("GRS_T2D","_",str_pad(9:104,width = 3,pad = "0"))

########################################
plot(9:104,dna_qual$Nucleic.Acid.ng.uL.)

plot(9:104,dna_qual$A260.A230)

plot(dna_qual$A260.A230,dna_qual$A260.A280,col=c("green","red")[(1*(dna_qual$Nucleic.Acid.ng.uL.<60))+1])

text(dna_qual$A260.A230,dna_qual$A260.A280,labels = as.character(9:104),cex=.6)
abline(h=c(1.85,1.88),lty=c(2,2),col="blue")
abline(v=c(2.3,2.8),lty=c(2,2),col="red")
abline(v=1.8,col="green",lty=2)
######################################

edge_samples<-as.character(c(1:8,89:96,9,17,25,33,41,49,57,65,73,81,16,24,32,40,48,56,64,72,80,88))

overall_18<-mean(1*(dna_qual$A260.A280<1.8&dna_qual$A260.A230<1.8))

edge_18<-mean(1*(dna_qual$A260.A280<1.8&dna_qual$A260.A230<1.8&(row.names(dna_qual)%in%edge_samples)))

nonedge_78<-mean(1*(dna_qual$A260.A280<1.8&dna_qual$A260.A230<1.8&!(row.names(dna_qual)%in%edge_samples)))

write.csv(dna_qual,"dna_qual.csv")

devtools::install_github("https://github.com/omegahat/RDCOMClient.git")

##################################
###### dna isolation slides 
library(DescTools)

GetCurrPP()

PpText("hello world", x = 100, y = 200, height = 50, width = 100, fontname = "Calibri", fontsize = 18,
       bold = FALSE, italic = FALSE, col = "black", bg = "grey" ,
       hasFrame = F, pp = DescToolsOptions("lastPP"))

txt_v<-apply(dna_qual, 1,function(x){
  
  c<-paste0(x[2],"\n")
  q1<-paste0(x[3],"\n")
  q2<-paste0(x[4],"\n" )
  paste0(c,q1,q2)
} )

lapply(txt_v[71:96], function(x){PpText(x, x = 1, y = 1, height = 50, width = 60, fontname = "Calibri", fontsize = 9,
                                 bold = FALSE, italic = FALSE, col = "black", bg = "grey" ,
                                 hasFrame = F, pp = DescToolsOptions("lastPP"))})


###################################

###confounding simulations 
###loading dna methylation data

BiocManager::install("WGCNA")

library(WGCNA)

##filtering probes in promoter regions 

library(recountmethylation)

#reading the manifest file 

library(data.table)

library(dplyr)

BiocManager::install("plyranges")
library(GenomicRanges)
library(plyranges)

library(ENmix)

annotation<-readmanifest("D:/jayasuriya_documents/r_projects/beta_values/manifest_file/MethylationEPIC_v-1-0_B4.csv")


load("D:/jayasuriya_documents/r_projects/beta_values/GenoSTAN_promoters.rda")


grep("E033",GenoSTAN_promoters@elementMetadata[,1])


GenoSTAN_promoters@ranges[GenoSTAN_promoters@ranges@NAMES=="E033"|GenoSTAN_promoters@ranges@NAMES=="E031"] %>% length()


tcells<-GenoSTAN_promoters[which(grepl(pattern="E033",GenoSTAN_promoters[,1]$name))]

bcells<-GenoSTAN_promoters[which(grepl(pattern="E031",GenoSTAN_promoters[,1]$name))]

####################################################
index_prom<-function(x,y,anno){
  pos<-anno[x,]$pos
  chr<-anno[x,]$chr
  seqname_chr<-paste0("chr",chr)
  print(seqname_chr)
  ys<-y %>% filter(seqnames==seqname_chr)
  k<-which(data.table::between(pos,start(ys),end(ys)))
  if(length(k)!=0){
    return(x)
  }else{return(NULL)}
}

index_promoter<-sapply(1:nrow(promoter),index_prom,y=tcells,anno=promoter)

index_promoter_b<-sapply(1:nrow(promoter),index_prom,y=bcells,anno=promoter)
####################################################

anno_subset<-annotation[[1]] %>% dplyr::select(c('chr','pos','Name','Strand'))

anno_subset<-anno_subset %>% mutate(end=pos,chr=paste0("chr",chr),Strand=ifelse(Strand=="F","+",ifelse(Strand=="R","-","*"))) %>% 
  rename(start=pos) %>% select(c(chr,start,end,Strand,Name))
  
cpgs<-makeGRangesFromDataFrame(anno_subset,na.rm = T,keep.extra.columns = T) 

tc_ov<-subsetByOverlaps(cpgs,tcells,type="within")
bc_ov<-subsetByOverlaps(cpgs,bcells,type="within")

c(tc_ov,bc_ov)[which(duplicated(c(tc_ov,bc_ov)))]


promoter_unique<-unique(c(tc_ov,bc_ov))

bt_specific<-unlist(mcols(promoter_unique),use.names = F)

###########################################
library(recountmethylation)
library(ExperimentHub)

recountmethylation::data_mdpost("remethdb2.h5")


library(HDF5Array)


test<-loadHDF5SummarizedExperiment("hdf5")

#######################################

test[row.names(test)%in%bt_specific,]

getBeta(test[row.names(test)%in%bt_specific,])

powers1=seq(1,20,by=1) 

RpowerTable=pickSoftThreshold(betas_n, powerVector=powers1)[[2]]

gc()

cex1=0.7 

par(mfrow=c(1,2)) 

plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],xlab="soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n") 


text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],labels=powers1,cex=cex1,col="red") 
# this line corresponds to using an R^2 cut-off of h 

abline(h=0.95,col="red") 

plot(RpowerTable[,1], RpowerTable[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n") 

text(RpowerTable[,1], RpowerTable[,5], labels=powers1, cex=cex1,col="red") 

Connectivity=softConnectivity(betas_n,power=10)-1

par(mfrow=c(1,1)) 

scaleFreePlot(Connectivity, main=paste("soft threshold, power=",10), truncated=F)

ConnectivityCut = 260 # number of most connected genes that will be considered 
# Incidentally, in the paper by Mischel et al (2005) we considered all 3600 #genes. 

ConnectivityRank = rank(-Connectivity) 

restConnectivity = ConnectivityRank <= ConnectivityCut 

# thus our module detection uses the following number of genes 
sum(restConnectivity) 

dissTOM=TOMdist(ADJ) 

ADJ= adjacency(betas_n[,restConnectivity],power=6) 

hierTOM = hclust(as.dist(dissTOM),method="average"); 

par(mfrow=c(1,1)) 

plot(hierTOM,labels=F) 

colorh1= cutreeStaticColor(hierTOM,cutHeight = 0.94, minSize = 5)

# The above should be identical to colorh1=datSummary$color1[restConnectivity] 
par(mfrow=c(2,1),mar=c(2,4,1,1)) 

plot(hierTOM, main="Cluster Dendrogram", labels=F, xlab="", sub=""); 

plotColorUnderTree(hierTOM,colors=data.frame(module=colorh1)) 

title("Module (branch) color") 

TOMplot(dissTOM , hierTOM, colorh1)

fwrite(list(colnames(betas_n)[which(colorh1=="blue")]),"gs.txt")

unique(colorh1)

means_cpg<-apply(betas_n,2,function(x){return(mean(as.numeric(x),na.rm=T))})
var_cpg<-apply(betas_n,2,function(x){return(var(as.numeric(x),na.rm=T))})
plot(unname(means_cpg),var_cpg,col=colorh1)




