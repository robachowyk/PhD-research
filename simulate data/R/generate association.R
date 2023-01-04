setwd("C:\\Users\\VNOB-0958\\Documents\\GitHub\\PhD-Research\\simulate data")

DF <- read.csv("DF_N=47075_2022-12-28.csv")
DF <- DF[!duplicated(DF),] # remove duplicates
DF <- na.omit(DF)
head(DF)

# generate treatment:
DF$treatment <- rbinom(nrow(DF),1,0.4)

# generate covariates
DF$X1 <- 2020 - DF$birth.year # include age

loc <- sample.int(10,1)
scale <- runif(1,0,2)
DF$X2 <- rnorm(nrow(DF),loc,scale)

loc <- sample.int(10,1)
scale <- runif(1,0,2)
DF$X3 <- rnorm(nrow(DF),loc,scale)

loc <- sample.int(5,1)+5
scale <- runif(1,0,0.5)
DF$X4 <- rnorm(nrow(DF),loc,scale)

loc <- sample.int(10,1)
scale <- runif(1,0,2)
DF$X5 <- rnorm(nrow(DF),loc,scale)

residual_errors <- rnorm(nrow(DF),0,1)
a <- sample.int(20,1)/10
b <- sample.int(30,1)/10
c <- sample.int(40,1)/10
d <- sample.int(10,1)/10

DF$Y <- - 10 + a*DF$treatment*DF$X1 + b*log(DF$X4) + c*DF$X3*DF$X2 + d*DF$X5

write.csv(DF, sprintf("DF_associations_N=%s_%s.csv", nrow(DF), Sys.Date()))
