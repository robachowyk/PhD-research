setwd("C:\\Users\\VNOB-0958\\Documents\\code")

# install.packages("comparator")
library(comparator)
# install.packages("dplyr")
library(dplyr)
# install.packages("plyr")
library(plyr)
# install.packages("NCmisc")
library(NCmisc)
# install.packages("sqldf")
library(sqldf)

DF <- read.csv("DF_associations_N=9397_2022-12-30.csv", header=TRUE, row.names=1)

identifiers <- list("family.name" = "jaro-winkler",
                    "gender" = "strict",
                    "country" = "strict",
                    "birth.year" = "large")

overlap <- DF[sample(nrow(DF),150),]

A <- rbind(DF[sample(nrow(DF),450),], overlap)
A <- A[,!names(A) %in% c("Y")]
rownames(A) <- NULL
head(A)

B <- rbind(DF[sample(nrow(DF),350),], overlap)[,append(names(identifiers),"Y")]
rownames(B) <- NULL
head(B)

levenshtein_similarity <- function(a,b){
  if (1 - elementwise(Levenshtein(),a,b)/max(length(a),length(b)) >= 0.95){return(1)} else{return(0)}
}

jaro_winkler_similarity <- function(a,b){
  if (elementwise(JaroWinkler(),a,b) >= 0.95){return(1)} else{return(0)}
}

strict_equality <- function(a,b){
  return(a==b)
}

large_equality <- function(a,b){
  return(substring(a,1,nchar(a)-1)==substring(b,1,nchar(b)-1)) # a revoir
}

comparison_vector <- function(A_record, B, identifiers) {
  methods <- list('jaro-winkler' = jaro_winkler_similarity, 'levenshtein' = levenshtein_similarity, 'strict' = strict_equality, 'large' = large_equality)
  comparisons <- list()
  for (linking_var in names(identifiers)) {
    method <- methods[[identifiers[[linking_var]]]]
    comparisons[[linking_var]] <- as.array(apply(B, 1, function(row) method(A_record[[linking_var]], row[[linking_var]])))
  }
  return(do.call(cbind, comparisons))
}
A_record <- A[1,]
comparison_vector(A_record, B, identifiers)

##############
# almost_exact_matches <- function(A, B, identifiers){
#   
#   A$source <- 'A'
#   B$source <- 'B'
#   A$true.index <- rownames(A)
#   B$true.index <- rownames(B)
#   # we remove A and B duplicates to ensure returning 1-2-1 true matches:
#   df <- rbind.fill(A[!duplicated(A[names(identifiers)], keep=FALSE),], B[!duplicated(B[names(identifiers)], keep=FALSE),])
#   #df <- merge(, B[!duplicated(B[names(identifiers)], keep=FALSE),], all=TRUE, join='inner')
#   duplicata <- df[duplicated(df[names(identifiers)], keep=FALSE),]
#   
#   
#   return( list( 'A' = as.array( to_vec(for(idx in duplicata) idx[1]) ), 'B' = as.array( to_vec(for(idx in duplicata) idx[2]) ) ) )
# }
# almost_exact_matches(A, B, identifiers)
# 
# 
# A$source <- 'A'
# B$source <- 'B'
# A$true.index <- rownames(A)
# B$true.index <- rownames(B)
# a <- A[!duplicated(A[names(identifiers)], keep=FALSE),append(names(identifiers),c("true.index","source"))]
# b <- B[!duplicated(B[names(identifiers)], keep=FALSE),append(names(identifiers),c("true.index","source"))]
# sqldf('SELECT * FROM a INTERSECT SELECT * FROM b')
# 
# semi_join(a, b, by = names(identifiers), copy = TRUE)
#           
# append(names(identifiers),c("true.index", "source"))
# 
# # we remove A and B duplicates to ensure returning 1-2-1 true matches:
# df <- rbind.fill(A[!duplicated(A[names(identifiers)], keep=FALSE),], B[!duplicated(B[names(identifiers)], keep=FALSE),])
# #df <- merge(, B[!duplicated(B[names(identifiers)], keep=FALSE),], all=TRUE, join='inner')
# duplicata <- df[duplicated(df[names(identifiers)]),]
# df[duplicated(df[names(identifiers)]),]
# df[duplicated(df[names(identifiers)], fromlast=TRUE),]
# 
# anyDuplicated(df[names(identifiers)])
# to_vec(df[names(identifiers)])
# dup.pairs(df[names(identifiers)])
