#BINF702 SPRING 2015 ASSIGNMENT #5
#Name - Sithalechumi Narayanan
#Date - April 27th 2015
#Name of file - assignment_5.r
#Part 1 (Based in part of Dr. Solka's Assignment #1 key)
#Part 2 (Based on both Dr. Solka's Chapter 8 Solutions and Krijnen Chapter 8 solutions )
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

#1.
library(multtest);data(golub); 

#recursive partitioning
library (rpart)

#to find the optimal gene with respect to Golub data
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
maxgol <- apply(golub.df[gol.fac=="ALL",!names(golub.df)%in%c("classification")], 2, function(x) max(x))
mingol <- apply(golub.df[gol.fac=="AML",!names(golub.df)%in%c("classification")], 2, function(x) min(x))
sum(maxgol < mingol)
which.min(maxgol - mingol)
golub.gnames[2124,]
boxplot(golub.df[,2124] ~gol.fac)

#Classification tree for the optimal gene
gol.rp <- rpart(gol.fac ~ golub.df[,2124], method="class", cp=0.001)
plot(gol.rp, branch=0,margin=0.1); text(gol.rp, digits=3, use.n=TRUE)

#To find gene Gdf5
grep("Gdf5",golub.gnames[,2])

#Classification tree for gene Gdf5
gol.rp <- rpart(gol.fac ~ golub.df[,2058], method="class", cp=0.001)
plot(gol.rp, branch=0,margin=0.1); text(gol.rp, digits=3, use.n=TRUE)
gol.rp <- rpart(gol.fac ~., data.frame(t(golub)), method="class", cp=0.001)
plot(gol.rp, branch=0,margin=0.1); text(gol.rp, digits=3, use.n=TRUE)

#2. Sensitivity vs specificity plot for gene expression CCND3 Cyclin D3
golub.clchanged <- -golub.cl -1
pred <- prediction(golub.df[,1042], golub.clchanged)
perf <- performance(pred, "sens", "spec")
plot(perf)

#b)It resembles the same function as figure 8.2

#c) Area under the curve sensitivity vs specificity curve
pred <- prediction(golub.df[,1042], golub.clchanged)
perf <- performance (pred, "auc")
perf

#4. Prediction of achieved remission for ALL data
#Constructing an expression set containing the patients with values
#on the phenotypical variable remission and the gene expressions
#with a significant p-value on the t-test with the patient groups CR
#or REF.
library(ALL); library(hgu95av2.db); library(rpart); data(ALL)
ALLrem <- ALL[,which(pData(ALL)$remission %in% c("CR","REF"))]
remfac <-factor(pData(ALLrem)$remission)
pano <- apply(exprs(ALLrem),1,function(x) t.test(x ~ remfac)$p.value)
names <- featureNames(ALLrem)[pano<.001]
ALLremsel<- ALLrem[names,]
data <- data.frame(t(exprs(ALLremsel)))
dim(data)
rownames(data)
colnames(data)
all.rp <- rpart(remfac ~., data, method="class", cp=0.001)
plot(all.rp, branch=0,margin=0.1); text(all.rp, digits=3, use.n=TRUE)
rpart.pred <- predict(all.rp, type="class")
table(rpart.pred,remfac)
remfac
rpart.pred 
7/(93+1+6+14)
mget(c("1840_g_at","36769_at","1472_g_at","854_at"), env = hgu95av2GENENAME)

#5. Gene selection by area under the curve 
library(ROCR); data(golub, package = "multtest")
gol.true <- factor(golub.cl,levels=0:1,labels= c("TRUE","FALSE"))
auc.values <- apply(golub.df,2,
function(x) performance(prediction(x, gol.true),"auc")@y.values[[1]])
o <- order(auc.values,decreasing=TRUE)
#To get the first 25 gene names
golub.gnames[o[1:25],2]
#yes CCND3 Cyclin D3 gene is one among them