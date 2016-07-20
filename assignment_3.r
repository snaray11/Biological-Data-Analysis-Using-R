#BINF702 SPRING 2015 ASSIGNMENT #3
#Name - Sithalechumi Narayanan
#Date - March 30th 2015
#Name of file - assignment_3.r
#Part 1 (Based in part of Dr. Solka's Assignment #1 key)
#Part 2 (Based on both Dr. Solka's Krijnen Chapter 3 and 4 HW Solutions and Krijnen Chapter 3 and 4 solutions )
#Merge the most relevant data found in the 3 tables (golub.gnames, golub, and golub.cl) that make-up the golub data in the library(multtest) into one data.frame with the following properties: 
  
#  Name: golub.df 
#Dimensions: patient rows and named gene columns, and an additional named column for the cancer classifications 
#Column Names: use the gene name (column 2) from golub.gnames and "classification" 
#Classification Column: use a factor column in golub.df that uses "ALL" and "AML" as the classifications 
#Part 1

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

#Part 2

#9. a) Three genes with the largest absolute t-values
data(golub, package="multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
tval <- apply(golub.df[gol.fac=="ALL",!names(golub.df)%in%c("classification")],2,function(x) sqrt(27) * mean(x)/sd(x))
o <- order(tval,decreasing=TRUE)
tval[o[1:3]]
golub.gnames[o[1:3],2]

# b) 
sdall <- apply(golub.df[gol.fac=="ALL",!names(golub.df)%in%c("classification")],2, sd)
sdaml <- apply(golub.df[gol.fac=="AML",!names(golub.df)%in%c("classification")],2, sd)
sdratio <- sdall/sdaml
sum( sdratio > 0.5 & sdratio < 1.5)

#8. a)

agrep("^Zyxin",golub.gnames[,2])
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
shapiro.test(golub.df[gol.fac=="ALL",2124])
#Since the p-value is smaller than the significance level, the null
#hypothesis is rejected and the Zyxin gene expression for ALL patients 
#is not normally distributed.
#to confirm
qqnorm(golub.df[gol.fac=="ALL",2124]) 
qqline(golub.df[gol.fac=="ALL",2124]) 

#b)

hist(golub.df[gol.fac=="ALL",2124],xlim = c(-2,2), ylim = c(0,10),xlab = 'Zyxin gene exp. ALL', ylab = 'Density', main = 'Data model for Zyxin gene exp. ALL patients')
curve(dnorm(x, 0.3, 0.75), add = TRUE)
#In the curve the mean value seems to lie at about 0.3, so the histogram 
#and the curve are somewhat similar.

#c)

gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
sigma <- 0.75; n <- 2; mu <- 0.3
x <- golub.df[gol.fac=="ALL",2124]
z.value <- sqrt(n)*(mean(x) - mu)/sigma

2*pnorm(-abs(z.value))
#Since the p-value is clearly larger than 0.05, we conclude that the null hypothesis
#is not rejected (accepted) and the Zyxin gene expression for ALL patients with 
#N(0.3, 0.75^2) distribution is different when compared to the Zyxin gene expression 
#for ALL patients itself. Therefore this distribution is perhaps based on N(0.3, 0.75^2)
#distribution


#Part 3

#1. a)
agrep("^CD33",golub.gnames[,2])
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
shapiro.test(golub.df[gol.fac=="ALL",808])
shapiro.test(golub.df[gol.fac=="AML",808])
#p-value for ALL = 0.592 and p-value for AML
#= 0.2583. Hence, for normality is accepted.

# b) 
var.test(golub.df[,808] ~ gol.fac) 
#gives p-value = 0.1095 so equality of variances is accepted.

# c) 
t.test(golub.df[,808] ~ gol.fac, var.equal = TRUE)
#gives p-value = 1.773e-09, so equality of means is rejected.

# d) Yes, t = -7.9813 is quite extreme.

#3. a)HOXA9. 

# a) 
shapiro.test(golub.df[gol.fac=="ALL",1391]) 
#gives p-value = 1.318e-07, so that normality is rejected.

# b) 
wilcox.test(golub.df[,1391] ~ gol.fac) 
#gives p-value = 7.923e-05, so that equality 
#of means is rejected. Note that the p-value from
#Grubbs test of the ALL expression values is 0.00519, so the null
#hypothesis of no outliers is rejected. Nevertheless the Welch two-
#sample T-test is also rejects the null-hypothesis of equal means.
#Its t-value equals -4.3026 and is quite large.

#6. 

gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
pt <- apply(golub.df[!names(golub.df)%in%c("classification")], 2, function(x) t.test(x ~ gol.fac)$p.value)
index <-agrep("^antigen",golub.gnames[,2])
golub.index<-golub[index,]
pt.index<-pt[index]
golub.gnames.index<-golub.gnames[index,]
golub.gnames.index[order(pt.index)[1:length(index)],2]

#8.

all66 <- golub.df[gol.fac=="ALL",66]
all790 <- golub.df[gol.fac=="ALL",790]
boxplot(all66,all790)
mean(all66);mean(all790)
median(all66);median(all790)
sd(all66);sd(all790)
IQR(all66)/1.349 ;IQR(all790)/1.349
mean(all66);mean(all790)
mad(all66);mad(all790)
shapiro.test(all66);shapiro.test(all790)

#10.

# a)
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
pt <- apply(golub.df[!names(golub.df)%in%c("classification")], 2, function(x) t.test(x ~ gol.fac)$p.value)
pw <- apply(golub.df[!names(golub.df)%in%c("classification")], 2, function(x) wilcox.test(x ~ gol.fac)$p.value)
o <- order(pt,decreasing=FALSE)
golub.gnames[o[1:10],2]

# b)
o <- order(pw,decreasing=FALSE)
golub.gnames[o[1:10],2]