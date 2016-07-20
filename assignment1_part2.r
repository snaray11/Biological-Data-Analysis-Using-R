#Name: Sithalechumi Narayanan
#Name of file: assignment1_part2.r
#Date: 02/09/2015


library (multtest)
data (golub)
golub.df <- data.frame(golub.gnames,golub) #merge golub.gnames and golub into  a data frame. Not merging golub.cl here
#because in the next step transpose will convert all the rows to columns and so we will then have to remove the golub.cl row and
# add it as a column. To save this step, merging only 2 datasets here.
View(golub.df)
tgolub.df <- t(golub.df) #transposes
View(tgolub.df)
rmrowsgolub.df <- tgolub.df [-c(1,3),] #removing rows 1 and 3
View(rmrowsgolub.df)

addcolgolub.df <- cbind(rmrowsgolub.df, golub.cl) #adding a column for cancer classification - golub.cl merge
View(addcolgolub.df)
addcolgolub.df [,3052]      #lets check what's found in column 3052
genenames.fac <- (addcolgolub.df[1,])    #create a factor for all the gene names
rownames (addcolgolub.df) <- NULL    #nullify the row and column names
colnames (addcolgolub.df) <- NULL
colnames (addcolgolub.df) <- genenames.fac   #have the gene names as column names
colnames (addcolgolub.df)[3052] <- "Classifications"   #have "classifications" as the last column name
View(addcolgolub.df)

#3. a) 
meangolub <- apply (rmrowsgolub.df,2,mean) #mean gene expression value
meangolub

#3. b)
o <- order(meangolub, decreasing=TRUE)
addcolgolub.df[o,]

#3. d)
golub.gnames [o[1:3],2] #biological names

#4. a)
sdgolub <- apply(rmrowsgolub.df,2,sd)

#4. b)
genesdgolub <- rmrowsgolub.df[sdgolub>0.5,]

#4. c)
sum (sdgolub>0.5) #gives the number of genes whose SD is more than 0.5
#1498

#5. a)
length(agrep("^oncogene", golub.gnames[,2]))
#agrep is used for string matching

#5.b)
rowindex <- agrep("^oncogene",golub.gnames[,2]) 
oncogolub <- golub[rowindex,]
oncogolub.gnames <- golub.gnames[rowindex,]
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
meangol <- apply(oncogolub[,gol.fac=="ALL"],1,mean)
o <- order(meangol,decreasing=TRUE)
oncogolub.gnames[o[1:3],2]

#5. c)
meangol <- apply(oncogol[,gol.fac=="AML"],1,mean)
o <- order(meangol,decreasing=TRUE)
oncogolub.gnames[o[1:3],2]

