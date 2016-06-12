LR.test <- function(data, alt.groups, null.groups, lib.size){
g <- dim(data)[1]
X2 <- NULL
log.fold.change <- NULL
for(i in 1:g){
		alt.model <- glm(as.numeric(data[i,]) ~ alt.groups, offset=log(lib.size), family=poisson)
		null.model <- glm(as.numeric(data[i,]) ~ null.groups, offset=log(lib.size), family=poisson)
		X2[i] <- deviance(null.model)-deviance(alt.model)
		log.fold.change[i] <- -alt.model$coef[2]
	}
chisq <- qchisq(df=1, (1:g-0.5)/g)
pval <- pchisq(X2, df=1, lower.tail=FALSE)
adj.pval <- p.adjust(pval, method="BH")

### Output ###
list(log.fold.change=log.fold.change, pvalues=pval, padj=adj.pval)
}
