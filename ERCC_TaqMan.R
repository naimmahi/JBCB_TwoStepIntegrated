
############### ERCC spike-in data validation ###########################

source("./Real_Data_Analysis.R")
path_to_write= "./outputDirectory/"
kNumOfMethods=6
method=c("edgeR", "DESeq2", "DSS")
ercc= read.table("./inputDirectory/ERCC_Controls_Analysis.txt", stringsAsFactors=FALSE, header=T, sep="\t", row.names=2)
load("./outputDirectory/result.RData")

deseq.nodg.ercc <- merge(result$res.deseq.nodg, ercc, by.x='row.names', by.y='row.names')
deseq.LR.ercc <- merge(result$res.LR.deseq, ercc, by.x='row.names', by.y='row.names')
edger.nodg.ercc <- merge(result$res.edger.nodg, ercc, by.x='row.names', by.y='row.names')
edger.LR.ercc <- merge(result$res.LR.edger, ercc, by.x='row.names', by.y='row.names')
dss.nodg.ercc <- merge(result$res.dss.nodg, ercc, by.x='row.names', by.y='row.names')
dss.LR.ercc <- merge(result$res.LR.dss, ercc, by.x='row.names', by.y='row.names')

plot.dat <- list()
plot.dat["Regular_DESeq2"] <- list(deseq.nodg.ercc)
plot.dat["Integrated_DESeq2"] <- list(deseq.LR.ercc)
plot.dat["Regular_edgeR"] <- list(edger.nodg.ercc)
plot.dat["Integrated_edgeR"] <- list(edger.LR.ercc)
plot.dat["Regular_DSS"] <- list(dss.nodg.ercc)
plot.dat["Integrated_DSS"] <- list(dss.LR.ercc)

qval.index <- list(Regular_DESeq2=7, Integrated_DESeq2=4, Regular_edgeR=3, Integrated_edgeR=4, Regular_DSS=3, Integrated_DSS=4) #col number of the padj in res
for(i in 1:length(method))
{
	setEPS()
	postscript(paste(path_to_write,"/ercc_roc_",method[i],".eps", sep=""),
	pointsize=20, fonts=c("serif","Helvetica"), horizontal=FALSE)
	plot_ercc(method[i])
	dev.off()
}

########################################## taqman RT-PCR data validation #################################################

load("./inputDirectory/TaqManData.Rdata")
taq <- taq.dat
deseq.nodg.taq <- merge(result$res.deseq.nodg, taq, by.x='row.names', by.y='row.names')
deseq.LR.taq <- merge(result$res.LR.deseq, taq, by.x='row.names', by.y='row.names')
edger.nodg.taq <- merge(result$res.edger.nodg, taq, by.x='row.names', by.y='row.names')
edger.LR.taq <- merge(result$res.LR.edger, taq, by.x='row.names', by.y='row.names')
dss.nodg.taq <- merge(result$res.dss.nodg, taq, by.x='row.names', by.y='row.names')
dss.LR.taq <- merge(result$res.LR.dss, taq, by.x='row.names', by.y='row.names')

plot.dat <- list()
plot.dat["Regular_DESeq2"] <- list(deseq.nodg.taq)
plot.dat["Integrated_DESeq2"] <- list(deseq.LR.taq)
plot.dat["Regular_edgeR"] <- list(edger.nodg.taq)
plot.dat["Integrated_edgeR"] <- list(edger.LR.taq)
plot.dat["Regular_DSS"] <- list(dss.nodg.taq)
plot.dat["Integrated_DSS"] <- list(dss.LR.taq)
qval.index <- list(Regular_DESeq2=7, Integrated_DESeq2=4, Regular_edgeR=3, Integrated_edgeR=4, Regular_DSS=3, Integrated_DSS=4) #col number of the padj in res

kLog2Cutoff <- 0.5 
for(i in 1:length(method))
{
	setEPS()
	postscript(paste(path_to_write,"/taq_AUC_",method[i],".eps", sep=""),
	pointsize=20, fonts=c("serif","Helvetica"), horizontal=FALSE)
	plot_AUC_taq(method[i])
	dev.off()
}

