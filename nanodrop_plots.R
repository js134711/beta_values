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
































