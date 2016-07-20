#Name: Sithalechumi Narayanan
#Name of file: assignment1_part1.r
#Date: 02/09/2015

library (multtest)
data (golub)
golub.df <- data.frame (cbind(golub.gnames,golub)) #merge golub.gnames and golub into  a data frame. Not merging golub.cl here
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
addcolgolub.df[1,3052]   #lets check if the last column name has changed

#addcolgolub.df <- cbind(addcolgolub.df, gol.fac)
#levels(0)=="TRUE" <- "ALL"  #check the values and input strings as either "ALL" or "AML"
#levels(1)=="TRUE" <- "AML"
#levels(0)=="FALSE"<- "N/A"
#levels(1)=="FALSE"<- "N/A"

