#BINF702 SPRING 2015 ASSIGNMENT #2
#Name - Sithalechumi Narayanan
#Date - March 2nd 2015
#Name of file - assignment2.r
#Part 1 (Based in part of Dr. Solka's Assignment #1 key)
#Part 2 (Based on both Dr. Solka's Krijnen Chapter 2 HW Solutions and Krijnen Chapter 2 solutions )
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

#Part 2
#2 Comparing two genes.
#a) Boxplot and outliers
i <- 66
boxplot(golub.df[,i])
i <- 790
boxplot(golub.df[,i])
#Column 790 has 3 outliers while column 66 doesn't have any

#b) QQ plot and formulating a hypothesis
qqnorm(golub.df[gol.fac=="ALL",i]) 
qqline(golub.df[gol.fac=="ALL",i]) 
#almost all values of 66 lie on or around the line, while the 3 outliers of 790
#are away from the normality line. Hypothesis:The expression values of 66 are 
#normally distributed, but those of column 790 are not.

#c) Mean and median for both the genes.
i <- 790
mean(golub.df[gol.fac=="ALL",i]) 
median(golub.df[gol.fac=="ALL",i])
#Mean (-1.174024) is larger than the median (-1.28137) due to
#outliers on the right hand side. 
i <- 66
mean(golub.df[gol.fac=="ALL",i]) 
median(golub.df[gol.fac=="ALL",i])
#For the gene in column 66 the mean is 1.182503 and 
#the median 1.23023. The differences are smaller.

#3. Effect size.
#a) Five genes with the largest effect size.
efs <- apply(golub.df[(gol.fac=="ALL"),!names(golub.df)%in%c("classification")],2,function(x) mean(x)/sd(x))
o <- order(efs,decreasing=TRUE)
efs[o[1:5]]
golub.gnames[o[1:5],2]

#b) Robust variant of the effect size
refs <- apply(golub.df[gol.fac=="ALL",!names(golub.df)%in%c("classification")],2,function(x) median(x)/mad(x))
o <- order(refs,decreasing=TRUE)
refs[o[1:5]]
golub.gnames[o[1:5],2]

#4. Gene expression for CCND3 Cyclin D3
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
stripchart(golub.df[,1042] ~ gol.fac,method="jitter")
stripchart(golub.df[,1042] ~ gol.fac,method="jitter",vertical = TRUE)
stripchart(golub.df[,1042] ~ gol.fac,method="jitter",col=c("red", "blue"),
vertical = TRUE)
stripchart(golub.df[,1042] ~ gol.fac,method="jitter",col=c("red", "blue"),
pch="*",vertical = TRUE)
title("CCND3 Cyclin D3 expression value for ALL and AMl patients")

#5. Box and Whisker plot for CCND3 Cyclin D3
x11()
data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
boxplot.stats(golub.df[(gol.fac=="ALL"),1042], coef = 1.5, do.conf = TRUE, do.out = TRUE) #finds values
boxplot(golub.df[(gol.fac=="ALL"),1042], xlim=c(0,4))
locator()
arrows(2.0, 1.93, 1.24, 1.93); text(2.5,1.93,"Median")
arrows(2.0, 1.1, 1.24, 1.1) ; text(2.5,1.1,"Outlier")
arrows(2.0, 1.79, 1.24, 1.79); text(2.5,1.79,"first quartile")
arrows(2.0, 2.17, 1.24, 2.17); text(2.5,2.17,"third quartile")
arrows(2.0, 1.27, 1.24, 1.27); text(2.5,1.27,"lower whisker")
arrows(2.0, 2.59, 1.24, 2.59); text(2.5,2.59,"upper whisker")
dev.copy2eps(device=x11,file="BoxplotWithExplanation.eps")

# 6. Box-and-whiskers plot of persons of Golub et al. (1999) data
#a) The medians are all around zero, the inter quartile range differ
#only slightly, the minimal values are all around minus 1.5. All
#persons have outliers near three.
boxplot(data.frame(golub.df))

#b) The means are very close to zero. The medians are all between
#(-0.15383, 0.06922), so these are also close to zero.
personmean <- apply(golub.df[!names(golub.df)%in%c("classification")],1,mean)
personmean
personmedian <- apply(golub.df[!names(golub.df)%in%c("classification")],1,median)
personmedian

#c) The data seem preprocessed to have standard deviation equal to
#one. The re-scaled IQR and MAD have slightly larger range.
range(apply(golub.df[!names(golub.df)%in%c("classification")],1,sd))
range(apply(golub.df[!names(golub.df)%in%c("classification")],1,function(x) IQR(x)/1.349))
range(apply(golub.df[!names(golub.df)%in%c("classification")],1,mad))

#7. Oncogenes of Golub et al. (1999) data.

#a) 
data(golub, package="multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
rowindex <- agrep("^oncogene",golub.gnames[,2])
oncogol <- golub[rowindex,]
oncogolub.gnames <- golub.gnames[rowindex,]
row.names(oncogol) <- oncogolub.gnames[,3]
boxplot(data.frame(t(oncogol[,gol.fac=="ALL"])))

#b) The plot gives a nice overview of the distributions of the gene
#expressions values of the onco gene separately for the ALL and
#the AML patients. Several genes behave similarly for ALL and
#AML. Some are clearly distributed around zero, but others not.
#Also, some have a small inter quartile ranges, while for others
#this is large. A similar statement holds for outliers, some do not
#have outliers, but others certainly have. Some gene show distinct
#distributions between patient groups. For instance, the sixth has
#ALL expressions around zero, but those for AML are larger than
#zero.
par(mfrow=c(2,1))
boxplot(data.frame(t(oncogol[,gol.fac=="ALL"])))
title("box-and-whiskers plot for oncogenes of ALL patients ")
boxplot(data.frame(t(oncogol[,gol.fac=="AML"])))
title("box-and-whiskers plot for oncogenes of AML patients ")
par(mfrow=c(1,1))

#8. Descriptive statistics for the ALL gene expression values of the Golub
#et al. (1999) data.

#a) The ranges indicate strong difference in means. The range of the
#mean and of the median are similar. The bulk of the data seems symmetric.
range(apply(golub.df[gol.fac=="ALL",!names(golub.df)%in%c("classification")],2,mean))
range(apply(golub.df[gol.fac=="ALL",!names(golub.df)%in%c("classification")],2,median))

#b) The range of the standard deviation is somewhat smaller than of
#the re-scaled IQR and MAD.
range(apply(golub.df[gol.fac=="ALL",!names(golub.df)%in%c("classification")],2,sd))
range(apply(golub.df[gol.fac=="ALL",!names(golub.df)%in%c("classification")],2,function(x) IQR(x)/1.349))
range(apply(golub.df[gol.fac=="ALL",!names(golub.df)%in%c("classification")],2,mad))

