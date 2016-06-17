library(seqc)
library(edgeR)
source("./Real_Data_Analysis.R")

############################# Data 1 #############################

dat3= ILM_refseq_gene_CNL[,-c(2,3,4)]
rownames(dat3)=dat3[,1]
dat3=dat3[,-1]
dat3=data.matrix(dat3)
dat2=dat3[, grepl("A_|B_", colnames(dat3))]
dat=dat2[, grep("_L01_FlowCellA", colnames(dat2))]
group1= length(grep("A_", colnames(dat)))
group2= length(grep("B_", colnames(dat)))
result <- Real_Data_Analysis(count_matrix=dat,length_group1=group1, length_group2=group2)
save(result, file= "./outputDirectory/result.RData")

############################# Data 2 #############################

katz_mouse= load(url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/katz_mouse_eset.RData"))
library(Biobase)
dat2=exprs(katz.mouse.eset)
dat=dat2[,c(2,4,1,3)]
group1= 2
group2= 2

result.katz <- Real_Data_Analysis(count_matrix=dat,length_group1=group1, length_group2=group2)
save(result.katz, file= "./outputDirectory/result.katz.RData")

