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
inp_c_file<-read_csv("C://Users/Jayasuriya/Desktop/dna_2806_gel.csv",col_types = c("c","d"))

#function for calculation of NFW, DNA and dilution steps needed 
#no dilution single dilution stepwise dilution 

#concentration is greater than 40 but less than 75 one microliter is taken 

#concentration is greater than 75 but less than 90 two microliters are taken and diluted to 50ng/ul

#concentration is greater than 90 but less than 150 diluted to half and one microliter is taken 

#greater than 150, serial diluted to half and two microliters is normalised to 50 ng 
sym(names(inp_c_file)[2])->conc_var

mt_inp_c_file<-inp_c_file |> mutate(c_in = case_when(!!conc_var<20~1,
                                                     !!conc_var>=20&!!conc_var<=40~2,
                                                     !!conc_var>40&!!conc_var<=75~3,
                                                     !!conc_var>75~4),
                                    rn = row_number()) |> 
  group_by(c_in) |> arrange(c_in,!!conc_var) |> filter(c_in!=1) |> ungroup()


dlk_cn<-function(c,concn){
  
  if(c == 2){
    sl<-(-9)
    dil_cn<-c(concn,0,2,3,1,concn*2,sl)
  }
  if(c==3){
    sl<-(-9)
    dil_cn<-c(concn,0,1,4,1,concn,sl)
  }
  if(c==4&concn<90){
    sl<-(-9)
    nvnv<-(concn/25)-2
    dil_cn<-c(concn,nvnv,2,4,1,50,sl)
  }
  if(c==4&concn>90&concn<150){
    sl<-(-9)
    dil_cn<-c(concn,1,1,4,1,concn/2,sl)
  }
  if(c==4&concn>=150){
    concn1<-concn/2
    sl<-99
    dil_cn_sec<-dlk_cn(4,concn1)[2:6]
    dil_cn<-c(concn,dil_cn_sec,sl)
  }
  return(dil_cn)
}

gel_md_array<-array(data = rep(0,26*8*number_of_gel),
                    dim = c(2,13,8,number_of_gel),
                    dimnames = list(x = c("R1","R2"), y = paste0("sample",1:13),
                                    z = c("SAMPLE_ID","NDC","DIL_VL","DGL_VL","N_VL","CL_VL","FIN_VL","SRL_DIL"),
                                    k = paste0("GL_NO:",1:number_of_gel))
)

for(i in rows_along(mt_inp_c_file)){
  
  sample_id<-mt_inp_c_file[i,4,drop=T];
  
  concn<-mt_inp_c_file[i,2,drop=T];
  
  indn<-mt_inp_c_file[i,3,drop=T];
  
  rnum<-ceiling(i/13);
  
  gelnum<-floor(i/25)+1;
  
  #print(gelnum)
  
  if(i>13){
    i<-i-13
  }
  
  gel_md_array[rnum,i,,gelnum]<-c(sample_id,dlk_cn(indn,concn));
  
  print(gel_md_array[rnum,i,,gelnum])
}


labels = apply(apply(aperm(gel_md_array,perm = c(2,1,3,4)) , 3, function(x)x),1,function(x){
  if(x[1]!=0){
    print(x[1])
    x1<-inp_c_file[as.integer(x[1]),,drop = T]
    print(x[2:8])
    y1<-paste0(c(x1,x[3:8]),"\n",collapse ="")
    
  }else{
    
    y1<-""
    
  }
  
})

xyplot(y~x,cast,
       scales = list(x = list(at = NULL), y= list(at = NULL)),
       panel = function(x,y,...){panel.rect(x = x, y = 20, width = 1, height = 5);
         panel.rect(x = x, y = 6, width = 1, height = 5)
         panel.text(x = c(x),y = rep(c(20),each=13),labels = c("lambda",labels[1:12]),fontsize = 8)
         panel.text(x = c(x),y = rep(c(6),each=13),labels = labels[13:26],fontsize = 8)
       }
)


library(bigsnpr)
library(purrr)
library(data.table)
library(foreach)
library(doParallel)
library(glue)
library(stringr)
library(lattice)
library(rtracklayer)       
library(GenomicRanges, help, pos = 2, lib.loc = NULL)


res_list<-function(path,c){fread(glue::glue("{path}/chr{c}.csv.gz"))}

index_get_meta<-function(chr_mg,ss_lst){              
                                        res_smch<-chr_mg[!duplicated(chr_mg,by="pos38",fromLast=T),.(chr = Chr,
                                                pos = pos37,
                                                a0=purrr::pmap_chr(.l=.(x=REF,y=ALT,z=A1),.f=function(x,y,z){ifelse(x==z,y,x)}),
                                                a1=A1)]
                                        tmp1<-imap(ss_lst,~merge_sst(res_t2d = res_smch,summ_stats= .x, prefix = .y))
                                        tmp11<-keep(tmp1,is.list)
                                        if(length(tmp11)>0){
                                            tmp2<-Reduce(function(x,y){merge(x,y,by="Position",all=T)},tmp11)
                                            return(merge(cbind(chr_mg[!duplicated(chr_mg,by="pos38",fromLast=T),!c("REF","ALT","A1")],"A0"= res_smch[,a0],"A1"= res_smch[,a1]),tmp2,by.x = "pos37",by.y="Position",all=T))
                                        
                                        }else{
                                            return(tmp11)
                                            }                  
                                }

merge_sst<-function(res_t2d,summ_stats,prefix){
    if(prefix == "mult_ans"){clnamestp<- c("Position","SE","PvalAssociation","PvalAncestryCorrelatedHeterogeneity","Ncases","Ncontrols","Neff")}else{clnamestp<-c("Position","SE","Pval","Ncases","Ncontrols","Neff")}
    merged_22<-summ_stats[,
                            {                tmp1=tryCatch(expr={snp_match(sumstats = .(chr =Chromsome,pos=Position,a0=NonEffectAllele,a1=EffectAllele,beta=Beta),info_snp = res_t2d,strand_flip=F,return_flip_and_rev = T,match.min.prop = 0)},error=function(e){return(NA)})
                                             if(is.list(tmp1[1])){tmp2=list(
                                                    beta=tmp1[,5],
                                                    #res_index = tmp1[,9],
                                                    eaf= .SD[tmp1[,6]][,ifelse(tmp1[,8],1-EAF,EAF)],
                                                    .SD[tmp1[,6]][,clnames,env=list(clnames=as.list(clnamestp))]) }else{tmp2=NA}                                    
  
                                       }
                                                ]
          if(is.list(merged_22)){clnm<-names(merged_22) 
               clnm<-setdiff(clnm,"Position")                                          
               setnames(merged_22,clnm,paste0(prefix,"_",tolower(clnm)),T)
               return(merged_22)}else{return(NA)}}


wrap_meta<-function(c){
     files_chr_wise<-lapply(c(eur = "./eur_files",sas = "./sas_files", mult_ans = "multians"),res_list,c=c)
     res_t2d<-sig[Chr==c]
     return(index_get_meta(chr_mg=res_t2d,ss_lst=files_chr_wise))
}                        
setwd("/home/storage/user14/files/t2d_gsa/summary_stat")
library(readxl)
loci_multans<-read_xlsx("41586_2024_7019_MOESM3_ESM.xlsx",sheet=25,range="A5:M640",col_names=F)
loci<-readRDS("loci_boundary.rds")
names(loci_multans)<-c(names(loci),paste0(rep(c("lead_snp","pos","pval"),times=3),rep(c("_mv","_dia","_glob"),each=3)))
setDT(loci_multans)
study_wise_index<-melt(loci_multans, id.vars = c("Locus", "Chr"),measure.vars = patterns("lead_snp_","pos_","pval_"),variable.name="study",value.name=c("rsid","pos","pval"))

get_pos<-function(res_rds=NULL,res_dt,match_vector,froms){
    #rds path contain the file path for rds file 
    #postions in match vector to be matched 
    if(!is.null(res_rds)){
        res_file<-if(froms){readRDS(res_rds)}else{res_rds}
        nsnps<-nrow(res_file[[1]])
        snps_pos<-res_file[[2]][5*((1:nsnps)-1)+1]
        matched_indx<-which(snps_pos%in%match_vector[-c(1)])
        pvalue_from_res<-res_file[[2]][(5*(matched_indx-1)+1),,drop=F]
        attr(pvalue_from_res,"dimnames")<-list(NULL,c("pos","wgen_or","wgen_beta_se","wgen_zstat","wgen_pval"))
        res_loci<-as.data.table(cbind(res_file[[1]][matched_indx,],pvalue_from_res))     
     
    }else{
        cs<-match_vector[1]
        pos_s<-match_vector[-c(1)]
        res_loci<-res_dt[`#CHROM`==cs&POS%in%pos_s]
        res_loci[,`:=`(TEST=NULL,OBS_CT=NULL)]
        res_loci[,`:=`(L95=NULL,U95=NULL)]
        setnames(res_loci,names(res_loci),c("Chr","pos38","id","REF","ALT","A1","wgen_or","wgen_beta_se","wgen_zstat","wgen_pval"))


    }
     return(res_loci)

}
pos_list_chrwise<-study_wise_index[,.(index_variant_pos=list(c(chr = unique(Chr) , pos = unique(na.omit(pos38))))),by=Chr]
pos_list_chrwise[,res:=list(purrr::map(index_variant_pos,function(x){get_pos(glue::glue("/home/storage/user14/files/backing_file/chr{x[1]}_2.TRAIT.glm.logistic.hybrid.rds"),x,T)}))]
index_res_table_locusk<-study_wise_index[!is.na(rsid)][!duplicated(pos38)][pos_list_chrwise[,res[[1]],by=Chr],on=.(Chr=Chr,pos38=pos)]
setnames(index_res_table_locusk,"pos","pos37",T)
setkey(index_res_table_locusk,Locus,pos38,pos37)

clust <- makeCluster(8) 
clusterExport(cl=clust,varlist=c("merge_sst","res_list","sig","index_get_meta"),envir=.GlobalEnv)
doParallel::registerDoParallel(clust)

foreach(i = 1:22,
        .packages = c("data.table","purrr","bigsnpr"),
        .verbose=T,
        .errorhandling = "pass") %dopar% {wrap_meta(c=i)} -> res_meta_index1

stopCluster(clust)
fwrite(res_meta_index1,"res_meta.csv")


#ld var 0.6
lst_index<-study_wise_index[,.(index_variants=list(.SD)),by=Locus]

ldsnps_topld<-function(gene_name,index_variant_dt,indictr){
    paste0("~/user14/files/t2d_gsa/summary_stat/ld_index/",gsub(" ","",gene_name),c("info","ld"),".txt")->lociname
    if(indictr){
            tmpfile<-tempfile(tmpdir="~/user14/files/t2d_gsa/summary_stat/rtemp_files",fileext=".txt")
            writeLines(unique(index_variant_dt[!is.na(rsid),paste0("chr",Chr,":",pos38)]),tmpfile,sep="\n")
            system(glue::glue("~/user14/softwares/topld_api/topld_api"," -inFile {tmpfile}", " -outputInfo {lociname[1]}"," -outputLD {lociname[2]}", " -pop EUR", " -thres 0.6"))
    }
    ldsnps_06<-fread(lociname[2],select=c(2:5,7,9,11,14,12),col.names=c("query_pos","target_pos","r2","dprime","target_rsid","target_maf","target_ref","target_alt","sign"))
    ldsnps_06_list<-split(ldsnps_06,by="query_pos")
    return(ldsnps_06_list)
}

lst_index[,ldvar:=list(purrr::pmap(.l=.(x=Locus,y=index_variants),function(x,y)ldsnps_topld(x,y,F)))]

rds_list<-lapply(1:22,function(x){k<-readRDS(glue::glue("/home/storage/user14/files/backing_file/chr{x}_2.TRAIT.glm.logistic.hybrid.rds"));return(k)})

lst_index_ld_res<-lst_index[,.(chr_pos=list(c(
              Chr=unique(na.omit(index_variants[[1]][,Chr])),
              pos=unique(na.omit(index_variants[[1]][,pos38])))),ldvar),by=Locus][,list(purrr::pmap(.(x=chr_pos,y=ldvar),function(x,y,rds=rds_list){
                ld_list_res<-lapply(y,function(k,chr=x[1],rds1=rds[[x[1]]]){get_pos(NULL,resf,c(chr,k[[2]]),F)})
                #names(ld_list_res)<-as.character(x[-c(1)])
                return(ld_list_res)
                 }
                )),by=Locus]


lst_index<-lst_index[lst_index_ld_res,on="Locus"]
#lst_index<-lst_index[loci_multans[,c(1:2)],on="Locus"]
setnames(lst_index,"V1","ldvar_res",T)
rtracklayer::import.chain('hg38ToHg19.over.chain')->C

for(i in 1:22){
    files_chr_wise<-lapply(c(eur = "./eur_files",sas = "./sas_files", mult_ans = "multians"),res_list,c=i)
    lst_index[Chr==i,`:=`(chr_mrg_res=purrr::pmap(.(x=Chr,y=ldvar_res),function(x,y,ss_lst=files_chr_wise){                    
                    ld_list_res<-lapply(y,function(k,chr=x){
                        while(nrow(k)>0){
                                chr_mg<-as.data.table(merge(x=list(Chr=chr),y=k))
                                lfted<-chr_mg[,as.data.table(rtracklayer::liftOver(GRanges(Rle(paste0("chr",Chr)),IRanges(start=pos38,width=1),score = pos38),C))][,.(pos37=start,pos38=score)]
                                chr_mg<-lfted[chr_mg,nomatch=0,on=.(pos38=pos38)]   
                                res_smch<-chr_mg[!duplicated(chr_mg,by="pos38",fromLast=T),.(chr = Chr,
                                        pos = pos37,
                                        a0=purrr::pmap_chr(.l=.(x=REF,y=ALT,z=A1),.f=function(x,y,z){ifelse(x==z,y,x)}),
                                        a1=A1)]
                                tmp1<-imap(ss_lst,~merge_sst(res_t2d = res_smch,summ_stats= .x, prefix = .y))
                                tmp11<-keep(tmp1,is.list)
                                if(length(tmp11)>0){
                                    tmp2<-Reduce(function(x,y){merge(x,y,by="Position")},tmp11)
                                    return(merge(cbind(chr_mg[!duplicated(chr_mg,by="pos38",fromLast=T),!c("REF","ALT","A1")],"A0"= res_smch[,a0],"A1"= res_smch[,a1]),tmp2,by.x = "pos37",by.y="Position",all=T))
                                
                                }else{
                                    return(tmp11)
                                    }                  
                            }
                        })
                    return(ld_list_res)
                    }
                    )),by=Locus]
                }
lst_index[,.I[nrow(rbindlist(chr_mrg_res[[1]]))==0],by=Locus][!is.na(V1),V1]->locus_no_ld_res
loci_merged<-rbindlist(lapply(1:594, function(i)lst_index[-c(locus_no_ld_res),][i,rbindlist(chr_mrg_res[[1]]),by=Locus]))
loci_merged<-loci_merged[!duplicated(loci_merged,by=c("Chr","pos38"))]

saveRDS(lst_index,"sas_lst_res.rds")
setkey(loci_merged,Locus,Chr,pos38,pos37)

#######################################
lst_ld<-lst_index[,.(merged_ld_res=list(imap(ldvar[[1]],~merge(x = .x,y = ldvar_res[[1]][[.y]],by.x =c("target_pos","target_rsid"),by.y=c("pos","ID"))))),by=.(Chr,Locus)]
ld_res_table_locusk<-lst_ld[,rbindlist(merged_ld_res[[1]]),by=.(Chr,Locus)]
setnames(ld_res_table_locusk,"target_pos","pos38",T)
setkey(ld_res_table_locusk,Locus,Chr,query_pos,pos38)


library(purrr, help, pos = 2, lib.loc = NULL)
find_comp<-function(x){
  c("A","T","G","C")->d1
  c("T","A","C","G")->d2
  y<-d2[which(d1%in%x)]
  return(y)  
}
std_f<-function(a1,a2,testa,testb,x){
  kl<-fcase(testa==a1&testb==a2,x,
                testa==a2&testb==a1,-x,
                find_comp(testa) == a1&find_comp(testb)==a2,x,
                find_comp(testb)==a2&find_comp(testa)==a1,-x,
                default=NA_real_)
  return(kl)
}
##index_var_table
##1050 var with res
index_res_table_locusk<-study_wise_index[!is.na(rsid)][!duplicated(pos38)][pos_list_chrwise[,res[[1]],by=Chr],on=.(Chr=Chr,pos38=pos38)]

study_wise_index[!is.na(rsid)][!duplicated(pos38)][,.N,by=.(study)]

setnames(index_res_table_locusk,"pos","pos37",T)

setkey(index_res_table_locusk,Locus,pos38,pos37)

study_wise_index[study=="global_t2d"][!is.na(rsid)][pos_list_chrwise[,res[[1]],by=Chr],nomatch=NULL,on=.(Chr=Chr,pos38=pos38)]

##ld_indexvar_table
#0.6 ld with res
lst_index[,.(list(names(ldvar[[1]])),list(names(ldvar_res[[1]]))),Locus][,all.equal(V1,V2)]
lst_ld<-lst_index[,.(merged_ld_res=list(imap(ldvar[[1]],~merge(x = .x,y = ldvar_res[[1]][[.y]],by.x =c("target_pos","target_rsid"),by.y=c("pos","ID"))))),by=.(Chr,Locus)]
ld_res_table_locusk<-lst_ld[,rbindlist(merged_ld_res[[1]]),by=.(Chr,Locus)]
setnames(ld_res_table_locusk,"target_pos","pos38",T)
setkey(ld_res_table_locusk,Locus,Chr,query_pos,pos38)
#number check
ld_res_table_locusk[,.N]
lst_index[,rbindlist(ldvar_res[[1]]),by=Locus][,.N]

index_res_table_locusk[Chr==1]
##
#res with index and meta
clust <- makeCluster(8) 
clusterExport(cl=clust,varlist=c("merge_sst","res_list","index_res_table_locusk","index_get_meta"),envir=.GlobalEnv)
doParallel::registerDoParallel(clust)

foreach(i = 1:22,
        .packages = c("data.table","purrr","bigsnpr"),
        .verbose=T,
        .errorhandling = "pass") %dopar% {wrap_meta(c=i)} -> res_meta_index1

stopCluster(clust)


##
res_meta_index1<-rbindlist(res_meta_index1)
setkey(res_meta_index1,Locus,pos38,pos37)
ld_res_table_locusk[res_meta_index1[res_meta_index1[,.I[is.na(sas_pval)]],.(Locus,Chr,pos38)]]

ld_res_table_locusk[,`:=`(a0=purrr::pmap_chr(.l=.(x=REF,y=ALT,z=A1),.f=function(x,y,z){ifelse(x==z,y,x)}),REF=NULL,ALT=NULL,A1=NULL,a1=A1)]

ld_res_table_locusk[,`:=`(r2_signed=pmap_dbl(.l=.(a0,a1,target_ref,target_alt,ifelse(sign=="+",r2,-r2)),std_f))]

res_meta_index1[duplicated(res_meta_index1,by="pos37")]
#res_meta_index1~index variants and their results ~ key(Locus,pos38,pos37)
#ld_res_table_locusk ~ ld variants r2 and res merged ~ key(Locus,Chr,query_pos,pos38)
#loci_merged ~ ld variants merged with results and meta ~ key(Locus,Chr,pos38,pos37)
#index, ld with minimum p value in wellgen, ld difference between eur and sas 

ld_res_table_locusk[res_meta_index1[.("TCF7L2",unique(pos38),114758349),.(Locus,Chr,pos38),nomatch=0]][order(pos38),levelplot(r2_signed~as.factor(pos38)*as.factor(pos38),col.regions=topo.colors(5),
                                                                                                       shrink=c(0.5,1),scales = list(x= list(rot=90)),at=pretty(c(0.6,1),4),
                                                                                                       colkey=list(at = pretty(c(0.6,1),4)))]

loci_merged[ld_res_table_locusk[res_meta_index1[61,.(Locus,Chr,pos38),nomatch=NULL]][,.(Locus,Chr,pos38)],nomatch=0][,
    .(sapply(.SD,function(x){-log10(as.numeric(x))}),pos38),.SDcols = c("wgen_pval","sas_pval","eur_pval")][,
    xyplot(wgen_pval+sas_pval+as.numeric(eur_pval)~pos38,cex=3,type=c("g","p"),pch=c(18,19,20),key=list(points = list(pch=c(18,19,20)),columns=3))]


ld_res_table_locusk[res_meta_index1[.("TCF7L2",unique(pos38),114758349),.(Locus,Chr,pos38),nomatch=0]][order(wgen_pval),.SD[1],by=.(r2_signed>0.9,r2_signed>0.8)]

ld_res_table_locusk[res_meta_index1[61,.(Locus,Chr,pos38),nomatch=NULL]][
    order(pos38),
    xyplot(-log10(wgen_pval)+I(max(-log10(wgen_pval))+r2_signed)~pos38-112998590,
    cex=3,type=c("p","p","g"),distribute.type=T,pch=c(18,19),
    par.settings = list(reference.line = list(col="darkgrey",alpha=1,lwd=1)))]

cls<-c("Chr","ID","A1_FREQ","wgen_or","wgen_pval","sas_eaf","sas_beta","sas_pval","eur_eaf","eur_beta","eur_pval","mult_ans_pvalancestrycorrelatedheterogeneity")

loci_merged[ld_res_table_locusk[res_meta_index1[17,.(Locus,Chr,pos38)],nomatch=NULL][,.(Locus,Chr,pos38)],nomatch=NULL][c(unique(sapply(list(wgen_pval,sas_pval,eur_pval),which.min))),.SD,.SDcols=cls]

fwrite(res_meta_index1,"res_meta.csv")



library(bigsnpr)
library(bigreadr)
library(data.table)
library(Matrix)
library(doParallel)
library(tracklayer)
library(collapse)
install.packages("logr")
library(logr)
install.packages("collapse")
chr5<-fread(input="./multians/chr5.csv.gz",nThread = 8)
#chr,pos37,pos38,a0,a1,a1frq,beta,se,ncase,ncontrol,neff

res_chr5<-resf[`#CHROM`==5]
rs5lfted<-res_chr5[,as.data.table(rtracklayer::liftOver(GRanges(Rle(paste0("chr",`#CHROM`)),IRanges(start=POS,width=1),score = POS),C))][,.(pos37=start,pos38=score)]
res_chr5<-merge(res_chr5,rs5lfted,by.x = "POS",by.y="pos38",all.x=T)
mrgd_res<-res_chr5[chr5,on=.(pos37=Position),nomatch = NULL]#384890


setwd("C:/user14/t2d/bfiles/")

geno_hg38<-snp_attach("sorted_hg38_bf.rds")

imput_mis<-snp_fastImputeSimple(geno_hg38$genotypes,method = "mean0")

ind.rows1k<-sample.int( n = 5877,size = 1000)
ind.cols<-sample.int(n=dim(imput_mis)[2],size = 1e5)
ind.rows1k<-ind.rows1k[order(ind.rows1k)]
ind.cols<-ind.cols[order(ind.cols)]


chr1_ind<-which(geno_hg38$map$chromosome==1)

ind_s_chr1<-ind.cols[(ind.cols%in%chr1_ind)]


m1<-snp_cor(imput_mis,
        ind.row = ind.rows1k,
        ind.col = ind_s_chr1,
        size    = 500,
        infos.pos = geno_hg38$map$physical.pos[ind_s_chr1]
        )

#sparsity
length(m@x)/m@Dim[1]/m@Dim[2]

setorder(ldmatrix,SNP1,SNP2)

m <- Matrix(nrow = length(keyld), ncol = length(keyld), data = 0, sparse = TRUE)

#m <- as(m, "CsparseMatrix")
m <- as(m, "symmetricMatrix")

keyld<-ldmatrix[,{
  tmp1<-unique(c(SNP1,SNP2))
  tmp2<-tmp1[order(tmp1)]
}]

ldmatrix[,`:=`(key1=match(SNP1,keyld),key2=match(SNP2,keyld))]
cl <- makeCluster(detectCores() - 1)  
registerDoParallel(cl)

for(i in seq_len(dim(ldmatrix)[1])){
  ind1<-ldmatrix[i,key1];
  ind2<-ldmatrix[i,key2];
  m[ind1,ind2]<-ldmatrix[i,R2]
}
#https://link.springer.com/article/10.1186/s12859-018-2289-9
corr<-as_SFBM(spmat = m, "./test_1e4")

ldmatrix[duplicated(ldmatrix,by=c("SNP1","SNP2"))|duplicated(ldmatrix,by=c("SNP1","SNP2"),fromLast=T)]

info<-fread("c://user14//ldfiles//sas//infochr1.csv.gz",nThread = 8,select = 1:5)
summ_stat_chr1<-fread("./sas_fls/chr1.csv.gz")
setkeyv(summ_stat_chr1,"Position")
# saveRDS(resf,"t2d_res_new.rds")
# rm(resf)
#########################################################
      #                                           #
      #               LIFTOVER START              #        
c1938<-import.chain("./GRCh37_to_GRCh38.chain")

lf<-log_open("chr1_liftover")

lfted<-summ_stat_chr1[,as.data.table(rtracklayer::liftOver(GRanges(Rle(paste0("chr",Chromsome)),IRanges(start=Position,width=1),score = Position),C1))][,.(pos38=start,pos37=score)]
summ_stat_chr1<-merge(summ_stat_chr1,lfted,by.x = "Position",by.y="pos37",all.x=T)
summ_stat_chr1[is.na(pos38),.N]
summ_stat_chr1[(duplicated(pos38,incomparables = NA)|duplicated(pos38,incomparables = NA,fromLast=T))]


function(chr_num,to_build,input_genome_build,input_dt,lfted_dt,posColName){
  
  s0<-paste0("Chromosome = ", chr_num)
  s1<-paste0("Number of positions given = ", nrow(input_dt))
  dp_rows_count<-input_dt[duplicated(pos),env = list(pos = posColName)]
  s2<-paste0("No of duplicate rows in input = ", dp_rows_count)

  
}

#                                           #
#               LIFTOVER END                # 

##########################################################

indx_ldinfo<-match(summ_stat_chr1$pos38,info$Position)


mtched_pos_ld<-info[na.omit(indx_ldinfo),Position]

ldmatrix<-fread("c://user14//ldfiles//sas_ld_chr1.csv.gz",select=1:5,nThread = 8)


fmatch(summ_stat_chr1$pos38,ldmatrix$SNP1)->s1m
fmatch(summ_stat_chr1$pos38,ldmatrix$SNP2)->s2m
intersect(na.omit(s1m),na.omit(s2m))->ld_match_matrix
ss(x = ldmatrix,i=ld_match_matrix)->ld_matrix_subset

#farthest non zero value in matrix/data.frame

res_chr1<-fread("./chr1.T2D.glm.logistic.hybrid",nThread = 8,select = c(1,2,3,4,5,7,9,12,16,17,19))
setnames(res_chr1,names(res_chr1)[1],"CHR",T)

summ_stat_chr1<-summ_stat_chr1[!is.na(pos38)]
names(summ_stat_chr1)

allele_match<-snp_match(sumstats = summ_stat_chr1[,.(chr=Chromsome,pos=pos38,a0=NonEffectAllele,a1=EffectAllele,beta=Beta,beta_se=SE,n_eff=Neff)],
          info_snp = res_chr1[,.(chr=CHR,pos=POS,a0=fcase(A1==REF,ALT,is.na(A1),NA,default=REF),a1=A1)],
          strand_flip = F,
          return_flip_and_rev = T,
          match.min.prop = 0,
          remove_dups = T)

fmatch(allele_match$pos,ldmatrix$SNP1)->s1m
fmatch(allele_match$pos,ldmatrix$SNP2)->s2m
intersect(na.omit(s1m),na.omit(s2m))->ld_match_matrix
ss(x = ldmatrix,i=ld_match_matrix)->ld_matrix_subset



setorder(ld_matrix_subset,SNP1,SNP2)


m <- Matrix(nrow = length(allele_match$pos), ncol = length(allele_match$pos), data = 0, sparse = TRUE)

#m <- as(m, "CsparseMatrix")
m <- as(m, "symmetricMatrix")

ld_matrix_subset[,`:=`(key1=match(SNP1,allele_match$pos),key2=match(SNP2,allele_match$pos))]

for(i in seq_len(dim(ld_matrix_subset)[1])){
  ind1<-ld_matrix_subset[i,key1];
  ind2<-ld_matrix_subset[i,key2];
  m[ind1,ind2]<-ld_matrix_subset[i,R2]
}


corr_chr1<-as_SFBM(m,"chr1",compact = T)



df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]


ldsc <- snp_ldsc(colSums(m), 
                    nrow(m)-1e6, 
                    chi2 = (allele_match$beta / allele_match$beta_se)^2,
                    sample_size = allele_match$n_eff, 
                    blocks = NULL)

h2_est <- ldsc[["h2"]]

h2_est

p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
h2_seq <- round(2.5 * c(0.7, 1, 1.4), 4)
grid.param <-
  expand.grid(p = p_seq,
              h2 = h2_seq,
              sparse = c(FALSE, TRUE))
# Get adjusted beta from grid model
beta_grid <-
  snp_ldpred2_grid(corr_chr1, allele_match, grid.param)
  











































