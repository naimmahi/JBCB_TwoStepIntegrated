OD.test <- function(data, alt.groups, lib.size){
g <- dim(data)[1]
score.test <- NULL
for(i in 1:g){
		alt.model <- glm(as.numeric(data[i,]) ~ alt.groups, offset=log(lib.size), family=poisson)
		hat.values <- hatvalues(alt.model)
		score.test[i] <- 1/(2*length(data[i,])) * sum(residuals(alt.model, type="pearson")^2 - ((data[i,] - hat.values*alt.model$fitted.values)/alt.model$fitted.values))^2
	}
chisq <- qchisq(df=1, (1:g-0.5)/g)
MSE <- 2
UL <- NULL
WH=0.05

#### Obtain the upper boundary of the WH bands #######################################
xbar <- mean(chisq)
bottom <- sum((chisq-xbar)^2)
top <- (chisq-xbar)^2
s <- sqrt(MSE*(1/g) + (top/bottom))
W <- sqrt(2*qf(df1=1, df2=g-1, p=1-(WH/g)))
UL <- pmax(chisq + W*s,1)

###### Obtain the indices of the over-dispersed and not-over-dispersed genes, respectively ##########

boundary <- min(which(sort(score.test)-UL > 0))
out <- boundary-1 + seq(boundary:length(score.test))
over.disp <- which(score.test %in% sort(score.test)[out])
not.over.disp <- setdiff(1:length(score.test), over.disp)

###### get splitted data ######
nodg=data[not.over.disp,]
odg=data[over.disp,]

### Output ###
list(index.over.disp=over.disp, index.not.over.disp=not.over.disp, test.statistic=score.test, chisq_quantile= chisq, nodg=nodg, odg=odg)
}
