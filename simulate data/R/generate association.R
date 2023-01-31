
# install.packages("reticulate")
# used for: pandas read pickle file
library(reticulate)

name_DF <- 'DF_N=4401_2023-01-16.csv'
DF <- read.csv(file.path("..", "..", "datasets", name_DF))
DF <- DF[!duplicated(DF),] # remove duplicates
# remove NaN values or empty elements:
DF <- DF[DF$name!="",]
DF <- DF[DF$family_name!="",]
DF <- DF[DF$was_assigned_female!="",]
DF <- DF[DF$country!="",]
DF <- DF[DF$birth_year!="",]
DF <- na.omit(DF)
head(DF)

Sys.which('python')
py_install("pandas") # say no to miniconda
use_virtualenv("r-reticulate")
source_python("pickle_reader.py")
name_dict <- 'country_proba_names_3226names.pkl'
country_proba_names <- read_pickle_file(file.path("..", "dictionaries", name_dict))

# generate treatment:
DF$treatment <- rbinom(nrow(DF),1,0.4)

# generate covariates
DF$X1 <- 2020 - DF$birth_year # include age

loc <- 2.5
scale <- 1
DF$X2 <- rnorm(nrow(DF),loc,scale)

loc <- 1
scale <- 1
DF$X3 <- rnorm(nrow(DF),loc,scale)

loc <- 1
scale <- 1
DF$X4 <- rnorm(nrow(DF),loc,scale)

loc <- 1
scale <- 1
DF$X5 <- rnorm(nrow(DF),loc,scale)

residual_errors <- rnorm(nrow(DF),0,1)
a <- 5.5
b <- 0.01
c <- 0.08
d <- 0.7

# for non-fixed treatment effect:
DF$Y <- - 10 + a*DF$treatment*DF$X2 + b*exp(DF$X4) + c*DF$X3*DF$X1 + d*DF$X5

head(DF)

write.csv(DF, file.path("..", "..", "datasets", sprintf("DF_associations_N=%s_%s.csv", nrow(DF), Sys.Date())), row.names=FALSE)
