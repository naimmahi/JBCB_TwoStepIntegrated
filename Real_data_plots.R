################################### Venn Diagram ###################################
source("./Real_Data_Analysis.R")
load("./outputDirectory/result.katz.RData")
path_to_write= "./outputDirectory/"
method=c("edgeR","DESeq2","DSS")
fdr=0.05
lfc=1	#log2
kNumOfMethods=6

plot.dat <- list()
plot.dat["Regular_edgeR"] <- list(result.katz$res.edger[which(result.katz$res.edger$padj<fdr & abs(result.katz$res.edger$log.fold.change)>lfc),])
plot.dat["Integrated_edgeR"] <- list(result.katz$res.integ.edger[which(result.katz$res.integ.edger$padj<fdr & abs(result.katz$res.integ.edger$log.fold.change)>lfc),])
plot.dat["Regular_DESeq2"] <- list(result.katz$res.deseq.odg[which(result.katz$res.deseq.odg$padj<fdr & abs(result.katz$res.deseq.odg$log.fold.change)>lfc),])
plot.dat["Integrated_DESeq2"] <- list(result.katz$res.integ.deseq[which(result.katz$res.integ.deseq$padj<fdr & abs(result.katz$res.integ.deseq$log.fold.change)>lfc),])
plot.dat["Regular_DSS"] <- list(result.katz$res.dss[which(result.katz$res.dss$padj<fdr & abs(result.katz$res.dss$log.fold.change)>lfc),])
plot.dat["Integrated_DSS"] <- list(result.katz$res.integ.dss[which(result.katz$res.integ.dss$padj<fdr & abs(result.katz$res.integ.dss$log.fold.change)>lfc),])
qval.index <- list(Regular_edgeR=3, Integrated_edgeR=3, Regular_DESeq2=3, Integrated_DESeq2=3, Regular_DSS=3, Integrated_DSS=3) #col number of the padj in res
lfc.index <- list(Regular_edgeR=1, Integrated_edgeR=1, Regular_DESeq2=1, Integrated_DESeq2=1, Regular_DSS=1, Integrated_DSS=1)

for(i in 1:length(method))
{
	setEPS()
	postscript(paste(path_to_write, method[i],"_venn_katz_",exp(log(2))*lfc,"fc.eps", sep=""),fonts=c("serif","Helvetica"), 
	pointsize=20, width=4, height=2.5, horizontal=FALSE)
	venn_diagram(vector1=rownames(plot.dat[[paste("Integrated_", method[i], sep="")]]),vector2=rownames(plot.dat[[paste("Regular_", method[i], sep="")]]), method[i])
	dev.off()
}


################################### ECDF plots ###################################

plot.dat <- list()
plot.dat["Regular_edgeR"] <- list(result.katz$res.edger.nodg["padj"])
plot.dat["Integrated_edgeR"] <- list(result.katz$res.LR.edger["padj"])
plot.dat["Regular_DESeq2"] <- list(result.katz$res.deseq.nodg["padj"])
plot.dat["Integrated_DESeq2"] <- list(result.katz$res.LR.deseq["padj"])
plot.dat["Regular_DSS"] <- list(result.katz$res.dss.nodg["padj"])
plot.dat["Integrated_DSS"] <- list(result.katz$res.LR.dss["padj"])
qval.index <- list(Regular_edgeR=1, Integrated_edgeR=1, Regular_DESeq2=1, Integrated_DESeq2=1, Regular_DSS=1, Integrated_DSS=1) #col number of the padj in res

for(i in 1:length(method))
{
	setEPS()
	postscript(paste(path_to_write, method[i],"_ecdf_NOD.eps", sep=""),fonts=c("serif","Helvetica"), 
	pointsize=20, width=8, height=6, horizontal=FALSE)
	Plot_ECDF(method[i], dat1=plot.dat[[paste("Regular_", method[i], sep="")]],dat2=plot.dat[[paste("Integrated_", method[i], sep="")]], qval.index)
	dev.off()
}

################################### GO enrichment ###################################

plot.dat <- list()
plot.dat["Regular_edgeR"] <- list(result.katz$res.edger)
plot.dat["Integrated_edgeR"] <- list(result.katz$res.integ.edger)
plot.dat["Regular_DESeq2"] <- list(result.katz$res.deseq)
plot.dat["Integrated_DESeq2"] <- list(result.katz$res.integ.deseq)
plot.dat["Regular_DSS"] <- list(result.katz$res.dss)
plot.dat["Integrated_DSS"] <- list(result.katz$res.integ.dss)
qval.index <- list(Regular_edgeR=2, Integrated_edgeR=2, Regular_DESeq2=2, Integrated_DESeq2=2, Regular_DSS=2, Integrated_DSS=2) #col number of the padj in res
lfc.index <- list(Regular_edgeR=1, Integrated_edgeR=1, Regular_DESeq2=1, Integrated_DESeq2=1, Regular_DSS=1, Integrated_DSS=1)

goseq= vector("list", kNumOfMethods)
for(i in 1: kNumOfMethods)
{
	print(i)
	goseq[[i]]= GO(plot.dat[[i]], qval.index[[i]], lfc.index[[i]])
	names(goseq)[i]= names(plot.dat[i])
}

for(i in 1:kNumOfMethods) 
{
	write.csv(goseq[[i]], file=paste(path_to_write,names(goseq[i]),".csv", sep="")) 
}

############# create a table of the above results #############

qval.index <- list(Regular_edgeR=2, Integrated_edgeR=3, Regular_DESeq2=6, Integrated_DESeq2=3, Regular_DSS=2, Integrated_DSS=3) #col number of the padj in res
res= vector("list", length(method))
for(i in 1:length(method))
{
	res[[i]] <- table_go(goseq[[grep(method[i], names(qval.index))[1]]],
	goseq[[grep(method[i], names(qval.index))[2]]])
}
GO_table <- rbind(res[[1]], res[[2]], res[[3]])
rownames(GO_table) = names(plot.dat)
GO_table

###################### Venn diagrams for GO ######################

setwd(path_to_write)
list.filenames = list.files(pattern="*.csv")
list.data<-list()
for (i in 1:length(list.filenames))
{
	list.data[[i]]<-read.csv(list.filenames[i])
}
names(list.data)<-list.filenames

plot.dat <- list()
plot.dat["Regular_edgeR"] <- list(list.data[[3]])
plot.dat["Integrated_edgeR"] <- list(list.data[[6]])
plot.dat["Regular_DESeq2"] <- list(list.data[[1]])
plot.dat["Integrated_DESeq2"] <- list(list.data[[4]])
plot.dat["Regular_DSS"] <- list(list.data[[5]])
plot.dat["Integrated_DSS"] <- list(list.data[[2]])

setdiff.reg= vector("list", length(method))
setdiff.integ= vector("list", length(method))
for(i in 1:length(method))
{
	setdiff.reg[[i]]= plot.dat[[paste("Regular_",method[i], sep="")]][match(setdiff(plot.dat[[paste("Regular_",method[i], sep="")]][,2], 
	plot.dat[[paste("Integrated_",method[i], sep="")]][,2]),as.character(plot.dat[[paste("Regular_",method[i], sep="")]][,2])),-1]

	setdiff.integ[[i]]= plot.dat[[paste("Integrated_",method[i], sep="")]][match(setdiff(plot.dat[[paste("Integrated_",method[i], sep="")]][,2],
	plot.dat[[paste("Regular_",method[i], sep="")]][,2]),as.character(plot.dat[[paste("Integrated_",method[i], sep="")]][,2])),-1]
	names(setdiff.reg)[i]= paste("Regular_", method[i], sep="")
	names(setdiff.integ)[i]= paste("Integrated_", method[i], sep="")
}
for(i in 1:length(method)) 
{
	write.csv(setdiff.reg[[i]], file=paste(path_to_write,names(setdiff.reg[i]),"_uni.csv", sep=""), row.names=FALSE) 
	write.csv(setdiff.integ[[i]], file=paste(path_to_write,names(setdiff.integ[i]),"_uni.csv", sep=""), row.names=FALSE) 
}

for(i in 1:length(method))
{
	setEPS()
	postscript(paste(path_to_write, method[i],"_venn_GO.eps", sep=""),fonts=c("serif","Helvetica"), 
	pointsize=20, width=4, height=2.5, horizontal=FALSE)
	venn_diagram(vector1=plot.dat[[paste("Integrated_", method[i], sep="")]][,2],vector2=plot.dat[[paste("Regular_", method[i], sep="")]][,2], method=NULL)
	dev.off()
}


