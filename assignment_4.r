#BINF702 SPRING 2015 ASSIGNMENT #4
#Name - Sithalechumi Narayanan
#Date - April 13th 2015
#Name of file - assignment_4.r
#Based in part of Dr. Solka's Assignment #1 key)
#Based on both Dr. Solka's Krijnen Chapter 7 HW Solutions and Krijnen Chapter 7 solutions )
#Merge the most relevant data found in the 3 tables (golub.gnames, golub, and golub.cl) that make-up the golub data in the library(multtest) into one data.frame with the following properties: 
  
#  Name: golub.df 
#Dimensions: patient rows and named gene columns, and an additional named column for the cancer classifications 
#Column Names: use the gene name (column 2) from golub.gnames and "classification" 
#Classification Column: use a factor column in golub.df that uses "ALL" and "AML" as the classifications 

data(golub, package = "multtest")

dim(golub)
#this returns 3051 by 38

#let's transpose the matrix
golubt = t(golub)

dim(golubt)

#this returns 38 by 3051

#this creates a data frame filled with NAs
golub.df = as.data.frame(matrix(nrow=38,ncol=3052))

#this fills the first 3051 columns of the data frame with the 
#transposed Golub data
golub.df[1:38,1:3051]=golubt

#here are the factors of classifications calculated on the 
#class file
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))

golub.df[,3052]=gol.fac

#here we create the column names
mycolnames = c(golub.gnames[1:3051,2],"classification")

#here we set the column names
colnames(golub.df) = mycolnames

dim(golub.df)
#this returns 38 by 3052 as expected

#Chapter 7 exercises.
#1. 

data <- data.frame(golub.df[,2124])
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
stripchart(golub.df[,2124]~gol.fac, pch=as.numeric(gol.fac))
plot(hclust(dist(data,method="euclidian"),method="single"))
initial <- matrix(tapply(golub.df[,2124],gol.fac,mean), nrow = 2, ncol=1) 
cl<- kmeans(data, initial, nstart = 10)
table(cl$cluster,gol.fac)
n <- length(data); nboot<-1000
boot.cl <- matrix(0,nrow=nboot,ncol = 2)
for (i in 1:nboot){
dat.star <- data[sample(1:n,replace=TRUE)]
cl <- kmeans(dat.star, initial, nstart = 10)
boot.cl[i,] <- c(cl$centers[1,],cl$centers[2,])
}
quantile(boot.cl[,1],c(0.025,0.975))
quantile(boot.cl[,2],c(0.025,0.975))

#2.
source("http://bioconductor.org/biocLite.R")
biocLite("genefilter")
library("genefilter"); data(golub, package = "multtest")
closeg <- genefinder(golub, 1042, 10, method = "euc", scale = "none")
golub.dataframe <- as.data.frame(x,row.names = golub.gnames, optional = FALSE)
#converted golub data matrix into a dataframe.
golub.gnames[closeg[[1]][[1]],2]
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
boxplot(golub.df[,394] ~gol.fac)

#3. 
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
x <- golub.df[,2289] 
y <- golub.df[,2430]
plot(x,y)
which.min(y) 
cor.test(x[-21],y[-21])
nboot <- 1000; boot.cor <- matrix(0,nrow=nboot,ncol = 1)
data <- matrix(c(x[-21],y[-21]),ncol=2,byrow=FALSE)
nboot <- 1000; boot.cor <- matrix(0,nrow=nboot,ncol = 1)
data <- matrix(c(x[-21],y[-21]),ncol=2,byrow=FALSE)
for (i in 1:nboot){
dat.star <- data[sample(1:nrow(data),replace=TRUE),]
boot.cor[i,] <- cor(dat.star)[2,1]}
mean(boot.cor)
quantile(boot.cor[,1],c(0.025,0.975))

#4. 
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
o1 <- grep("oncogene",golub.gnames[,2])
plot(hclust(dist(golub.df[,o1],method="euclidian"),method="single"))
o2 <- agrep("^antigen",golub.gnames[,2])
plot(hclust(dist(golub.df[,o2],method="euclidian"),method="single"))
o3 <- grep("receptor",golub.gnames[,2])
plot(hclust(dist(golub.df[,o3],method="euclidian"),method="single"))

#5.
library(ALL); data(ALL)
ALLB <- ALL[,ALL$BT %in% c("B1","B2","B3")]
ALLB.df <- as.data.frame(x,row.names = NULL, optional = FALSE) 
#converted into a dataframe
panova <- apply(exprs(ALLB), 1, function(x) anova(lm(x ~ ALLB$BT))$Pr[1])
ALLBsp <- ALLB[panova<0.001,]
dim(exprs(ALLBsp))
min(cor(exprs(ALLBsp)))
eigen(cor(exprs(ALLBsp)))$values[1:5]
data <- exprs(ALLBsp); p <- ncol(data); n <- nrow(data) ; nboot<-1000
eigenvalues <- array(dim=c(nboot,p))
for (i in 1:nboot){dat.star <- data[sample(1:n,replace=TRUE),]
eigenvalues[i,] <- eigen(cor(dat.star))$values}
for (j in 1:p) print(quantile(eigenvalues[,j],c(0.025,0.975)))
biplot(princomp(data,cor=TRUE),pc.biplot=T,cex=0.5,expand=0.8)