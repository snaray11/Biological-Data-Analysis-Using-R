#BINF702 SPRING 2015 Project R Script
#Name - Sithalechumi Narayanan
#Date - May 4th 2015
#Name of file - epilepsy_dataset.r
#Based on "A handbook of Statistical Analyses Using R" - Chapter 6, 
#Chapter 11 and Chapter 13 - Brian S. Everitt and Torsten Hothorn 

#R Documentation
data("epilepsy", package = "HSAUR")
library(lattice)
str(epilepsy)
dotplot(I(seizure.rate / base) ~ period | subject, data = epilepsy,
subset = treatment == "Progabide")

#Finding mean and variance
itp <- interaction(epilepsy$treatment, epilepsy$period)
tapply(epilepsy$seizure.rate, itp, mean)
tapply(epilepsy$seizure.rate, itp, var)

install.packages("geepack")

#Boxplot for Placebo and Progabide with the number of seizures
layout(matrix(1:2, nrow = 1))
ylim <- range(epilepsy$seizure.rate)
placebo <- subset(epilepsy, treatment == "placebo")
progabide <- subset(epilepsy, treatment == "Progabide")
boxplot(seizure.rate ~ period, data = placebo,
ylab = "Number of seizures",
xlab = "Period", ylim = ylim, main = "Placebo")
boxplot(seizure.rate ~ period, data = progabide,
main = "Progabide", ylab = "Number of seizures",
xlab = "Period", ylim = ylim)

#To find Generalized Linear Models GLM
library(gee) 
library(MASS) 
data(OME)
per <- rep(log(2),nrow(epilepsy))
epilepsy$period <- as.numeric(epilepsy$period)
fm <- seizure.rate ~ base + age + treatment + offset(per)
epilepsy_glm <- glm(fm, data = epilepsy, family = "poisson")

#To find Generalized Estimating Equation GEE
epilepsy_gee1 <- gee(fm, data = epilepsy, family = "poisson",
id = subject, corstr = "independence", scale.fix = TRUE,
scale.value = 1)
epilepsy_gee2 <- gee(fm, data = epilepsy, family = "poisson",
id = subject, corstr = "exchangeable", scale.fix = TRUE,
scale.value = 1)
epilepsy_gee3 <- gee(fm, data = epilepsy, family = "poisson",
id = subject, corstr = "exchangeable", scale.fix = FALSE,
scale.value = 1)

#Boxplot for Placebo and Progabide with Log number of seizures
layout(matrix(1:2, nrow = 1))
ylim <- range(log(epilepsy$seizure.rate + 1))
boxplot(log(seizure.rate + 1) ~ period, data = placebo,
main = "Placebo", ylab = "Log number of seizures",
xlab = "Period", ylim = ylim)
boxplot(log(seizure.rate + 1) ~ period, data = progabide,
main = "Progabide", ylab = "Log number of seizures",
xlab = "Period", ylim = ylim)

#Summary details for all the estimate values
summary(epilepsy_glm)
summary(epilepsy_gee1)
summary(epilepsy_gee2)
summary(epilepsy_gee3)
