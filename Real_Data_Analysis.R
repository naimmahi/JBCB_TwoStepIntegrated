###################################################
###### Real Data analysis functions #########
###################################################

Real_Data_Analysis <- function (count_matrix, length_group1, length_group2, result_directory) 
{

	library(edgeR)
	library(DESeq2)
	library(DSS)

	source(paste(directory,"OD.test.R", sep=""))
	source(paste(directory,"LR.test.R", sep=""))

############################### OD test ###############################
## We assume data set is already filtered 
	o <- order(rowSums(count_matrix))
	reads <- count_matrix[o,]
	keep <- rowSums(cpm(reads)>1) >= min(length_group1,length_group2)
	counts <- data.matrix(reads[keep,])

	col=c(rep(1, length_group1), rep(0, length_group2))
	alt.groups=factor(col)
	nod.dgelist = DGEList(counts = counts, group = factor(col))
	res.nod.odt = OD.test(data = counts, alt.groups, lib.size = nod.dgelist$samples$lib.size)
	nodg=res.nod.odt$nodg	#set of NOD genes
	odg=res.nod.odt$odg	#set of OD genes

############################### Regular edgeR ###############################

	data= counts
	d.reg= DGEList(counts=data, group=col)
	d.reg= calcNormFactors(d.reg)
	d.reg= estimateCommonDisp(d.reg)
	d.reg= estimateTagwiseDisp(d.reg)
	result.edger= exactTest(d.reg)
	fdr.edger= p.adjust(result.edger$table$PValue, method = "BH", n = length(result.edger$table$PValue))
	res.edger= data.frame(log.fold.change=result.edger$table$logFC, pvalues= result.edger$table$PValue, padj= fdr.edger)
	rownames(res.edger)=rownames(data)

############### OD edgeR ###################

	data.odg= odg
	d.odg= DGEList(data.odg, group=col)
	d.odg= calcNormFactors(d.odg)
	d.odg= estimateCommonDisp(d.odg)
	d.odg= estimateTagwiseDisp(d.odg)
	result.edger.odg= exactTest(d.odg)
	fdr.edger.odg= p.adjust(result.edger.odg$table$PValue, method = "BH", n = length(result.edger.odg$table$PValue))
	res.edger.odg= data.frame(log.fold.change=result.edger.odg$table$logFC, pvalues=result.edger.odg$table$PValue, padj= fdr.edger.odg)
	rownames(res.edger.odg)=rownames(odg)
	colnames(res.edger.odg)=c("log.fold.change","pvalues","padj")

############### NOD edgeR ###################

	data.nodg= nodg
	d.nodg= DGEList(data.nodg, group=col)
	d.nodg= calcNormFactors(d.nodg)
	d.nodg= estimateCommonDisp(d.nodg)
	d.nodg= estimateTagwiseDisp(d.nodg)
	result.edger.nodg= exactTest(d.nodg)
	fdr.edger.nodg= p.adjust(result.edger.nodg$table$PValue, method = "BH", n = length(result.edger.nodg$table$PValue))
	res.edger.nodg= data.frame(log.fold.change=result.edger.nodg$table$logFC,pvalues=result.edger.nodg$table$PValue, padj= fdr.edger.nodg)
	rownames(res.edger.nodg)=rownames(nodg)

############ LR test in edgeR ################

	alt.groups=factor(col)
	null.groups <- rep(1, times=length(counts[1,]))
	nod.dgelist = DGEList(counts = nodg, group = factor(col))
	norm.dgelist = calcNormFactors(nod.dgelist, method = "TMM")
	lib.sizes.edger = as.vector(norm.dgelist$samples$norm.factors) *as.vector(norm.dgelist$samples$lib.size)
	res.LR.edger = LR.test(data = nodg, alt.groups,null.groups, lib.size = lib.sizes.edger)
	res.LR.edger =as.data.frame(res.LR.edger )
	rownames(res.LR.edger)=rownames(nodg)

############ Combine OD&NOD genes for edgeR ################

	res.integ.edger.pre <- rbind(res.LR.edger,res.edger.odg )
	res.integ.edger <- res.integ.edger.pre[rownames(res.edger),]

############################### Regular DESeq2 ###############################

	data=counts
	colData <- data.frame(condition=factor(col))
	dds <- DESeqDataSetFromMatrix(data, colData, formula(~ condition))
	dds <- DESeq(dds)
	result.deseq <- as.data.frame(results(dds))
	fdr.deseq= p.adjust(result.deseq$pvalue, method = "BH", n = length(result.deseq$pvalue))
	res.deseq= result.deseq[,c(2,5)]
	res.deseq$padj= fdr.deseq
	colnames(res.deseq)=c("log.fold.change","pvalues","padj")

############### OD DESeq2 ###################

	dds.odg <- DESeqDataSetFromMatrix(odg, colData, formula(~ condition))
	dds.odg <- DESeq(dds.odg)
	result.deseq.odg <- as.data.frame(results(dds.odg))
	fdr.deseq.odg= p.adjust(result.deseq.odg$pvalue, method = "BH", n = length(result.deseq.odg$pvalue))
	res.deseq.odg= result.deseq.odg[,c(2,5)]
	res.deseq.odg$padj= fdr.deseq.odg
	colnames(res.deseq.odg)=c("log.fold.change","pvalues","padj")

############### NOD DESeq2 ###################

	colData <- data.frame(condition=factor(col))
	dds.nodg <- DESeqDataSetFromMatrix(nodg, colData, formula(~ condition))
	dds.nodg <- DESeq(dds.nodg)
	result.deseq.nodg <- as.data.frame(results(dds.nodg))
	fdr.deseq.nodg= p.adjust(result.deseq.nodg$pvalue, method = "BH", n = length(result.deseq.nodg$pvalue))
	res.deseq.nodg= result.deseq.nodg[,c(2,5)]
	res.deseq.nodg$padj= fdr.deseq.nodg
	colnames(res.deseq.nodg)=c("log.fold.change","pvalues","padj")

############ LR test in DESeq2 ###############

	alt.groups=factor(col)
	null.groups <- rep(1, times=length(counts[1,]))
	nod.dgelist = DGEList(counts = nodg, group = factor(col))
	colData <- data.frame(condition=factor(col))
	norm.deseq <- DESeqDataSetFromMatrix(nodg, colData, formula(~ condition))
	norm.deseq2=estimateSizeFactors(norm.deseq)
	norm.deseq3 =sizeFactors(norm.deseq2)
	lib.sizes.deseq = as.vector(norm.deseq3)*as.vector(nod.dgelist$samples$lib.size)
	res.LR.deseq = LR.test(data = nodg, alt.groups, null.groups, lib.size = lib.sizes.deseq)
	res.LR.deseq =as.data.frame(res.LR.deseq )
	rownames(res.LR.deseq)=rownames(nodg)

############ Combine OD&NOD genes for DESeq2 ################

	res.integ.deseq.pre <- rbind(res.LR.deseq,res.deseq.odg )
	res.integ.deseq <- res.integ.deseq.pre[rownames(res.deseq),]

############################### Regular dss ###############################

	data= counts
	colnames(data)= c(1:dim(data)[2])
	d=newSeqCountSet(data, col)
	dge <- estNormFactors(d)
	dge2 <- estDispersion(dge)
	result.dss <- waldTest(dge2,0,1)
	fdr.dss= p.adjust(result.dss$pval, method = "BH", n = length(result.dss$pval))
	res.dss= data.frame(log.fold.change=result.dss$lfc, pvalues=result.dss$pval, padj= fdr.dss)
	rownames(res.dss)=rownames(data)

################ OD dss ###################

	data.nb= data.matrix(odg)
	colnames(data.nb)= c(1:dim(data.nb)[2])
	d.nb= newSeqCountSet(counts=data.nb, col)
	d.nb= estNormFactors(d.nb)
	dge.nb <- estDispersion(d.nb)
	result.dss.odg <- waldTest(dge.nb,0,1)
	fdr.dss.odg= p.adjust(result.dss.odg$pval, method = "BH", n = length(result.dss.odg$pval))
	res.dss.odg= data.frame(logfc=result.dss.odg$lfc, pvalues=result.dss.odg$pval, padj= fdr.dss.odg)
	rownames(res.dss.odg)=rownames(odg)
	colnames(res.dss.odg)=c("log.fold.change","pvalues","padj")

################ NOD dss ###################

	data.nodg= data.matrix(nodg)
	colnames(data.nodg)= c(1:dim(data.nodg)[2])
	d.nodg= newSeqCountSet(counts=data.nodg, col)
	d.nodg= estNormFactors(d.nodg)
	dge.nodg <- estDispersion(d.nodg)
	result.dss.nodg <- waldTest(dge.nodg,0,1)
	fdr.dss.nodg= p.adjust(result.dss.nodg$pval, method = "BH", n = length(result.dss.nodg$pval))
	res.dss.nodg= data.frame(log.fold.change=result.dss.nodg$lfc, pvalues=result.dss.nodg$pval,padj= fdr.dss.nodg)
	rownames(res.dss.nodg)=rownames(nodg)

############### LR test in DSS ###########

	data.nodg= data.matrix(nodg)
	colnames(data.nodg)= c(1:dim(data.nodg)[2])
	nod.dgelist = newSeqCountSet(data.nodg, col)
	norm.factor = estNormFactors(nod.dgelist)
	norm.factor=normalizationFactor(norm.factor)
	lib.sizes.dss = as.vector(norm.factor) *as.vector(colSums(data.nodg))
	res.LR.dss = LR.test(data = data.nodg, alt.groups,null.groups, lib.size = lib.sizes.dss)
	res.LR.dss =as.data.frame(res.LR.dss )
	rownames(res.LR.dss)=rownames(nodg)

############ Combine OD&NOD genes for DSS ################

	res.integ.dss.pre <- rbind(res.LR.dss,res.dss.odg )
	res.integ.dss <- res.integ.dss.pre[rownames(res.dss),]

#################### result list #######################

	result <- list(od.test=res.nod.odt,
	res.deseq=res.deseq,res.integ.deseq=res.integ.deseq, res.deseq.odg=res.deseq.odg, res.deseq.nodg=res.deseq.nodg, res.LR.deseq=res.LR.deseq,
	res.edger=res.edger,res.integ.edger=res.integ.edger, res.edger.odg=res.edger.odg, res.edger.nodg=res.edger.nodg, res.LR.edger=res.LR.edger,
	res.dss=res.dss, res.integ.dss=res.integ.dss, res.dss.odg=res.dss.odg, res.dss.nodg=res.dss.nodg, res.LR.dss=res.LR.dss)

	return(result)
}



###########################################################################
###### Sensitivity and specificity with spike-ins and taqman data #########
###########################################################################

############### ERCC spike-in data validation ###########################

PlotRocs_ERCC <- function(i, dat, qval.index, logmix.index, main)
{
  outcome= rep(1, dim(dat)[1])
  outcome[dat[,logmix.index] == 0] =0
  #ylab=NULL
  if(i==3){ylab="Sensitivity"}
  else {ylab=NA}
  if(i==1|i==3|i==5){
    roc <- plot.roc(outcome, dat[,qval.index],lty=2,lwd=4.5, ylab=ylab, main=main) #main="ROC of ERCC spike-in data"

  }else{
		roc <- lines.roc(outcome, dat[,qval.index], add=TRUE, lwd=4.5,lty=1)
	}
  return(roc)
}
plot_ercc <- function(method)
{
	require(pROC)
	#par(mgp=c(1.4,1.5,0),mar=c(5.1,5,4,2))
	res <- lapply(seq(kNumOfMethods)[grep(method, names(qval.index))], function(i) PlotRocs_ERCC(i, plot.dat[[i]],qval.index[[i]],
	dim(plot.dat[[i]])[2], main=method))
	names(res) <- names(plot.dat)[grep(method, names(qval.index))]
	legends <- lapply(c(1,2), function(i) paste(names(res)[i], "AUC =", format(res[[i]]$auc, digits=3), sep=' '))
	legend("bottomright", legend=legends, lty=c(2,1), lwd=2,  cex=.75)
}

#################### TaqMan data validation ###########################(unused)

PlotRocs_Taqman <- function(i, dat, qval.index, logFC.index)
{
  require(pROC)
  outcome= rep(1, dim(dat)[1])
  outcome[abs(dat[,logFC.index]) <= kLog2Cutoff] =0
  if(i==1|i==3|i==5)
  {
    roc <- plot.roc(outcome, dat[,qval.index], lty=2,lwd=4.5, main="ROC of TaqMan data", ylim=c(0,1.05))
    mtext(paste("logFC cutoff= ", kLog2Cutoff, sep=''), side=3, padj=-1.75, cex=.8)

  }else{
    roc <- lines.roc(outcome, dat[,qval.index], lty=1,lwd=4.5,add=TRUE)
  }
  return(roc)
}
plot_taq <- function(method)
{
	require(pROC)
	#par(mgp=c(1.4,1.5,0),mar=c(5.1,5,4,2))
	res <- lapply(seq(kNumOfMethods)[grep(method, names(qval.index))], function(i) PlotRocs_Taqman(i, plot.dat[[i]],qval.index[[i]],dim(plot.dat[[i]])[2]))
	names(res) <- names(plot.dat)[grep(method, names(qval.index))]
	legends <- lapply(c(1,2), function(i) paste(names(res)[i], "AUC =", format(res[[i]]$auc, digits=3), sep=' '))
	legend("bottomright", legend=legends, lty=c(2,1), lwd=2,  cex=.75)
}

########## Calculate AUCs by changing log2 cutoff ##############

Calculate_AUC <- function(i, dat, qval.index, logFC.index)
{
  ## calculate ROC
  ## return AUC vector for a range of logFC cutoffs
  require(pROC)
  auc.res <- matrix(nrow=length(seq(0.5,2,0.1)), ncol=1)

  ## logFC cutoff range
  cutoff <- seq(0.5,2,0.1)

  for(i in seq(1:length(cutoff)))
  {
    outcome <- rep(1, dim(dat)[1])
    outcome[abs(dat[,logFC.index]) <= cutoff[i]] =0
    auc.res[i] <- roc(outcome, dat[,qval.index])$auc[[1]]
  }
  return(auc.res)
}

plot_AUC_taq <- function(method)
{
	auc.res <- sapply(seq(kNumOfMethods)[grep(method, names(qval.index))], function(i) Calculate_AUC(i, plot.dat[[i]], qval.index[[i]], dim(plot.dat[[i]])[2])) 
	colnames(auc.res) <- names(plot.dat)[grep(method, names(qval.index))]
	  if(method=="edgeR"){ylab="AUC"}
	  else {ylab=NA}
	  par(mar=c(5.1, 4.1, 2.1, 2.1))
	plot(seq(0.5,2,0.1), auc.res[,1], type='n', xlab="logFC cutoff values", ylab=ylab, main=method, ylim=c(0.4,1)) #c(min(auc.res)-0.1, max(auc.res)+0.1))
	lines(seq(0.5,2,0.1), auc.res[,1],lwd=4.5, lty=2, ylim=c(0.45,0.85))
	lines(seq(0.5,2,0.1), auc.res[,2],lwd=4.5, lty=1, ylim=c(0.45,0.85))
	legend("topright", legend=colnames(auc.res), lty=c(2,1), lwd=2,  cex=.75)
}

############################ ECDF plots ##############################

Plot_ECDF <- function (method, dat1,dat2, qval.index) 
{
	require(sSeq)
	d= data.frame(dat1[,as.numeric(qval.index[grep(method, names(qval.index))[1]])], dat2[,as.numeric(qval.index[grep(method, names(qval.index))[2]])])
	colnames(d)= c(names(qval.index[grep(method[1], names(qval.index))[1]]), names(qval.index[grep(method[1], names(qval.index))[2]]))

	  if(method=="edgeR"){ylab="ECDF"}
	  else {ylab=NA}
	par(mgp=c(1.4,1.5,0),mar=c(5.1,5,4,2))  
	ECDF <- ecdfAUC(d
	, xlab="p-value",main="",ylab=ylab,lineType=c(5,1),addLeg=FALSE, drawRef=TRUE,
	cex.axis=1.5, cex.lab=1.7,col.line= c("black","black"),
	axis.padj=c(-0.5, 0.5), lab.padj=c(-3, 3), lwd=c(3.5,3.5), box.lwd=1.2)
	legends <- lapply(c(1,2), function(i) paste(colnames(d)[i], "AUC= ", format(as.numeric(ECDF[i]), digits=3), sep=' '))
	legend("bottomright", legend=legends, lty=c(2,1), lwd=2.5,  cex=0.75)
	return(ECDF)
}


####################################### Enrichment Analysis ############################################

GO <- function(dat, qval.index, lfc.index) 
{

	require(goseq)
	library("org.Mm.eg.db")
	library(GO.db)

	genes <- ifelse(dat[, qval.index]< 0.05 & abs(dat[, lfc.index])>0.5, 1, 0)
	names(genes)=row.names(dat)
	pwf=nullp(genes,"mm10","ensGene")
	GO.wall=goseq(pwf,"mm10","ensGene")
	enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue,method="BH")<.05]
	sig.GO = GO.wall[match(enriched.GO,GO.wall$category),c(1,2,6,7)]
	return(sig.GO)
}

table_go <- function (data1, data2)
{
	regular_uniq= setdiff(data1[,1],data2[,1])
	integrated_uniq= setdiff(data2[,1],data1[,1])
	common= length(intersect(data1[,1],data2[,1]))
	Regular = c( dim(data1)[1], length(regular_uniq), common)
	Integrtaed = c(dim(data2)[1],length(integrated_uniq), common)
	final = rbind(Regular, Integrtaed)
	colnames(final)= c("Total categories", "Unique categories", "Common categories")
	return(final)
}

venn_diagram <- function(vector1, vector2, method=NULL)
{
	library(VennDiagram)
	venn= venn.diagram(list(vector1,vector2),cex=1.2,lwd=1,NULL, main=method,main.cex=1.2,main.pos= c(0.5, 1.05), category.names=c(" ", " ")) 
	return(grid.draw(venn))
}
