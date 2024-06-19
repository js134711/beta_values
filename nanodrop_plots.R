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

#

library(HardyWeinberg)
library(bigsnpr)
library(bigstatsr)
library(dplyr)

gwas_samples<-snp_attach("/home/storage/user14/files/t2d_gsa/flp_allelq/maf_allele_recode_samples.rds")

grl<-gwas_samples$fam %>% group_by(family.ID) 

str_l<-group_rows(grl)


names(str_l)<-unlist(group_keys(grl),use.names = F)

hwe_max<-function(ind,genot,grow,all_count_thres){
  
  gindex<-genot[,ind]
  
  chisq_list_ind<-lapply(grow, function(x,gs = gindex,thr = all_count_thres){
    
    kl<-apply(gs[x,],2,FUN = function(y){
      
      count_t<-table(y);
      
      y_cr_chisq<-(-99L)
      
      if(length(count_t) == 3){
        a_allele_count<-((count_t[1]*2+count_t[2])<thr);
        b_allele_count<-((count_t[3]*2+count_t[2])<thr);
        if(a_allele_count|b_allele_count){
          y_cr_chisq<-(-9L)
        }else{
          count_t<-as.vector(count_t);
          names(count_t)<-c("AA","AB","BB");
          y_cr_chisq<-HardyWeinberg::HWChisq(count_t,verbose = F)$chisq;
        }
      }
      
      if(length(count_t)==2){
         
          n_cou<-names(count_t)
          rep_name<-dplyr::case_match(n_cou,"0"~"AA","1"~"AB","2"~ "BB")
          count_t<-as.vector(count_t)
          names(count_t)<-rep_name
          len_chk<-dplyr::case_when(is.na(count_t["AA"])~any(c(count_t["AB"]+2*count_t["BB"],count_t["AB"])<thr),
                    is.na(count_t["BB"])~any(c(count_t["AB"]+2*count_t["AA"],count_t["AB"] )<thr),
                    is.na(count_t["AB"])~any(c(count_t["AA"],count_t["BB"])<thr))
          
          count_t<-dplyr::case_when(is.na(count_t["AA"])~c("AA"=0,count_t),
                           is.na(count_t["BB"])~c("BB"=0,count_t),
                           is.na(count_t["AB"])~c("AB"=0,count_t))
          
          if(len_chk){
            y_cr_chisq<-(-8L)
          }else{
            y_cr_chisq<-HardyWeinberg::HWChisq(count_t,verbose = F)$chisq;
          }
      }
      
       return(list(counts1 = count_t, hwe_chisq = y_cr_chisq))})
    #print(kl)
    return((kl))
  })
  
  return(do.call(cbind,chisq_list_ind))
}

options(bigstatsr.check.parallel.blas = FALSE)

res_hwe<-bigstatsr::big_apply(X=gwas_samples$genotypes,
                                a.FUN = hwe_max,
                                all_count_thres = 5,
                                ind = cols_along(gwas_samples$genotypes),
                                grow = str_l,
                                ncores = 15,
                                a.combine = "rbind")

save(res_hwe,file = "hwe_geno_maf_fil.RData")

dev_hwe<-function(x2_groups,number_group,min_x2,atleast_minx2){
  
  hwe_diseq<-apply(x2_groups, 1, function(x){
    x<-sapply(x,function(y){return(y[[2]])})
    sum_group<-sum(x>=min_x2)
    max_group<-sum(x>=atleast_minx2)
    if(sum_group>=number_group&max_group>=1){
      return(T)
    }else{
      return(F)
    }
  })
  return(which(hwe_diseq))
}


diseq<-dev_hwe(x2_groups = res_hwe,number_group = 4, min_x2 = 4, atleast_minx2 = 8)

apply(res_hwe[diseq,],1,unlist,recursive=T,simplify = T)

View(res_hwe[diseq,])


snp_stratfrq<-function(gen,ind,fam,map){
    print(dim(gen[,ind]))
    mer_gen<-cbind(fam,gen[,ind])
    print(dim(mer_gen))
    sum_frq<-mer_gen%>%group_by(family.ID)%>%summarise(across(as.character(1):tail(colnames(mer_gen),1),~mean(.x,na.rm=T)/2))
    print(dim(sum_frq))
    all_frq_t<-as_tibble(cbind(nms = names(sum_frq), t(sum_frq)))[,-1]
    colnames(all_frq_t)<-all_frq_t[1,]
    all_frq_t<-all_frq_t[-1,]
    loci_mer_atgc<-cbind(map[ind,],all_frq_t)
    return(loci_mer_atgc)
}


ran_all_frq<-big_apply(gwas_samples$genotypes,a.FUN = snp_stratfrq,
    ind  = which(apply(res_hwe,1,function(x){if(any(x==-99)){return(T)}else{FALSE}})),
    map = gwas_samples$map,
     fam = gwas_samples$fam)
################################################################################################
#merge r2 value from vcf to plink result
library(bigstatsr)
library(bigreadr)
library(foreach)
library(stringr)
library(dplyr)
vcf_list<-list.files('~/T2D_GWAS',pattern = 'vcf.gz')

lapply(vcf_list,function(x){
    ofile<- str_replace(x,pattern = ".*(chr\\d+).*",replacement = "\\1")
    system(glue::glue('bcftools query --format "%CHROM\\t%POS\\t%INFO/R2\\n" {x} >  {ofile}.r2'))
})
setwd("~/T2D_GWAS/r2_val")
r2_file<-list.files(pattern = ".r2")

r2_fbm<-function(x1){
  fln<-paste0(sample(x = c(letters),size = 8,replace = F),collapse = "")
  y<-as_FBM(x = x1,backingfile = paste0("./",fln),type = "double")
  return(y)
}

fbm_tr_fun<-function(d){
  chr_num<-str_replace(d,"chr(\\d+).r2","\\1")
  fbm_files<-big_fread1(file = d,every_nlines = 1e5,.transform = r2_fbm,select = 2:3)
  return(list(cn = chr_num, fbm = fbm_files))
}

r2_bind_com<-function(number,input_file,r2_tr){
  fbm_bds<-r2_tr[[which(vapply(r2_tr, function(x){x[[1]]==number}, FUN.VALUE = logical(1)))]][[2]]
  input_file$R2<-NA_real_
  input_file<-input_file[CHR == number]
  s_r2<-r2_bind_fun(index = 1,z = fbm_bds, res = input_file)
  return(s_r2)
}

r2_bind_fun<-function(index,z,res){
            r2_table<-as.data.table(z[[index]][,])            
            names(r2_table)<-c("POS","R2")  
            r2_table<-unique(r2_table,by="POS")
            res<-res%>%rows_update(r2_table,by = "POS",unmatched = "ignore")          
            if(!any(is.na(res$R2))){
                return(res)
            }else{
                 res<-try(r2_bind_fun(index+1,z = z,res=res))
                 return(res)
            }
                }

start1<-Sys.time()
r2_fbm_trsfmd<-foreach(d = r2_file)%do%fbm_tr_fun(d)
r2_res_fl<-foreach(a = 1:22,.combine = 'rbind')%do%r2_bind_com(number = a, input_file = res_fl,r2_tr = r2_fbm_trsfmd)
end1<-Sys.time()

end1-start1 #22.00813



























