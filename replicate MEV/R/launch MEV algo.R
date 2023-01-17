
### DO SIMULATIONS AND REPLICATE MEV ALGORITHM ###

# install.packages("comparator")
# used for: elementwise
library(comparator) 

# install.packages("comprehenr")
# used for: to_vec
library(comprehenr) 

# install.packages("testit")
# used for: assert
library(testit)

# install.packages("dplyr")
# used for: %>%
library(dplyr)

# install.packages("RoughSets")
# used for: D.discretize.quantiles.RST
library(RoughSets)

# install.packages("ggplot2")
# used for: ggplot
library(ggplot2)

# install.packages("Hmisc")
# used for: mean_ci... confidence intervals
library(Hmisc)

### ALGORITHM ###

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
  # Compare one record in A with all record in B.
  # Return the binary comparison of the identifiers for one record in A with all records in B.
  
  # A_record:     series of one row,
  # B:            dataframe,
  # identifiers:  list dict: k = column name,
  #                          v = method in {'large', 'strict', 'levenshtein', 'jaro-winkler'}
  methods <- list('jaro-winkler' = jaro_winkler_similarity, 'levenshtein' = levenshtein_similarity, 'strict' = strict_equality, 'large' = large_equality)
  comparisons <- list()
  for (linking_var in names(identifiers)) {
    method <- methods[[identifiers[[linking_var]]]]
    comparisons[[linking_var]] <- as.array(apply(B, 1, function(row) method(A_record[[linking_var]], row[[linking_var]]) ))
  }
  return(do.call(cbind, comparisons))
}

compute_max_linking_score <- function(A_record, B, identifiers, match, unmatch) {
  # Compute the argmax and max linking score of records in B with A_record.
  
  # A_record:     series of one row,
  # B:            dataframe,
  # identifiers:  list dict: k = column name,
  #                          v = method in {'large', 'strict', 'levenshtein', 'jaro-winkler'},
  # match:        array of proba of having same linking variables when being a match,
  # unmatch:      array of proba of having same linking variables (at all, among the nA x nB pairs of record)
  similarities <- comparison_vector(A_record, B, identifiers)
  if (sum(apply(similarities, 1, all)) == 1) { # there is one unique match: probably a true match
    return(c(rownames(B[apply(similarities, 1, all),]), NA))
  } else {
    linking_score <- apply( t(apply(similarities, 1, function(row) row * log2(match/unmatch))) + t(apply(1-similarities, 1, function(row) row * log2((1-match)/(1-unmatch)))), 1, sum)
    return(c(which.max(linking_score), max(linking_score))) # assign the first element of the argmax set
  }
}

linking_score <- function(A, B, identifiers, match, unmatch) {
  # Compare records in A with records in B, computing all linking scores for records in A with records in B.
  # Return the indices of records in A with the best match index for record in B.
  
  # A:            dataframe,
  # B:            dataframe,
  # identifiers:  list dict: k = column name,
  #                          v = method in {'large', 'strict', 'levenshtein', 'jaro-winkler'},
  # match:        array of proba of having same linking variables when being a match,
  # unmatch:      array of proba of having same linking variables (at all, among the nA x nB pairs of record)
  links <- apply(A, 1, function(row) compute_max_linking_score(row[names(identifiers)], B, identifiers, match, unmatch))
  idx_in_A <- as.array(1:nrow(A))
  idx_in_B <- as.array(as.numeric(links[1,]))
  matching_scores <- as.array(as.numeric(links[2,]))
  max_score <- max(matching_scores, na.rm = TRUE)
  matching_scores[is.na(matching_scores)] <- max_score + 1
  return(list("A"=idx_in_A, "B"=idx_in_B, "scores"=matching_scores))
}

stratified_ATE <- function(DF_group, pop_size) {
  # Compute the average treatment effect on the specific stratum represented in DF_group.
  
  # DF_group:  dataframe,
  # pop_size:  size of the entire population from which the group comes from.
  n_treated <- nrow(DF_group[DF_group$treatment == 1,])
  n_untreated <- nrow(DF_group[DF_group$treatment == 0,])
  assert("One group among treated/untreated is empty with this stratification", {
    (n_treated!=0)
    (n_untreated!=0)
  })      
  avg_outcome_treated <- sum(DF_group$treatment * DF_group$Y) /  n_treated
  avg_outcome_untreated <- sum((1-DF_group$treatment) * DF_group$Y) /  n_untreated
  var_treated <- sum((DF_group$treatment * DF_group$Y - avg_outcome_treated)^2) / (n_treated-1)
  var_untreated <- sum(((1-DF_group$treatment) * DF_group$Y - avg_outcome_untreated)^2) / (n_untreated-1)
  size_group <- nrow(DF_group)
  ATE_group <- (size_group/pop_size) * (avg_outcome_treated - avg_outcome_untreated)
  variance_group <- (size_group/pop_size)^2 * (var_treated/n_treated + var_untreated/n_untreated)
  return(c(ATE_group, variance_group))
}

logit <- function(p) {
  return(log(p/(1-p)))
}

minmaxscaler <- function(X){
  return((X - min(X)) / (max(X) - min(X)))
}

propensity_score <- function(DF, covariates, scaler, convert_to_logit){
  # Compute propensity score estimates: the probability (logistic regression) that an observation is treated or not conditioned on some covariates.
  
  # DF:                dataframe,
  # covariates:        list of strings for the names of covariates variables in DF,
  # scaler:            scaler function,
  # convert_to_logit:  boolean for converting probabilities to logit when building the propensity score estimates based on a logistic regression.
  if (is.null(scaler)) {
    model <- glm(treatment ~ ., family="binomial", data=DF[append(covariates,'treatment')])
    #model <- glm(treatment ~ X1, family="binomial", data=DF) # binomial(link='logit')
  } else {
    DF[covariates] <- scaler(DF[covariates])
    model <- glm(treatment ~ ., family="binomial", data=DF[append(covariates,'treatment')])
  }
  fitted.results <- predict(model, newdata=DF[append(covariates,'treatment')],type='response')
  if (convert_to_logit) {
    fitted.results <- logit(fitted.results)
    return(fitted.results)
  }
  return(fitted.results)
}

ATE <- function(DF, strata){
  # Compute the average treatment effect in DF according to the stratification method:
  # no stratification when strata is NULL,
  # stratified dataframe build based on the variables given (should be modified manually in the 3rd part of the function),
  # or propensity score stratification.
  # propensity score estimates and propensity score quantiles assignment should be already built for this last method.
  
  # DF:      dataframe,
  # strata:  value in {NULL, 'propensity score', [...]}
  #          [...] variables to compute the strata on.
  pop_size = nrow(DF)
  if (is.null(strata)) {
    return(stratified_ATE(DF, pop_size))
  } else if ((length(strata) == 1) && (strata == 'propensity stratification')) {
    assert("For propensity score stratification you need first to add the propensity_score and prop_Score_quantile columns into the dataframe.", {
      ("propensity_score" %in% names(DF))
      ("prop_score_quantile" %in% names(DF))
    })
    ATE <- 0
    variance <- 0
    for (q in unique(DF$prop_score_quantile)){
      stratum_data <- DF[DF$prop_score_quantile == q,]
      stratified_ATE_variance <- stratified_ATE(stratum_data, pop_size)
      ATE_stratum <- stratified_ATE_variance[1]
      variance_stratum <- stratified_ATE_variance[2]
      ATE = ATE + ATE_stratum
      variance = variance + variance_stratum} 
    return(c(ATE, variance))} else {
      ATE <- 0
      variance <- 0
      for (stratum_data in split(DF, interaction(strata))){
        stratified_ATE_variance <- stratified_ATE(stratum_data, pop_size)
        ATE_stratum <- stratified_ATE_variance[1]
        variance_stratum <- stratified_ATE_variance[2]
        ATE = ATE + ATE_stratum
        variance = variance + variance_stratum}
      return(c(ATE, variance))}
}

Minimum_Estimated_Variance <- function(A, B, identifiers, match, unmatch){
  # Compute the Average Treatment Effect (and the variance) of each subset defined by linking scores levels.
  # Algorithm based on 
  # Returns the indices of records in A with the best match index for record in B.
  
  # DF:      dataframe,
  # strata:  value in {NULL, 'propensity score', [...]}
  #          [...] variables to compute the strata on.
  ATE_list <- c()
  variance_list <- c()
  matchings <- linking_score(A, B, identifiers, match, unmatch)
  for (score in rev(sort(unique(matchings$scores)))){
    best_matches <- matchings$scores>=score
    from_A <- A[best_matches,]
    from_B <- B[matchings$B[best_matches],]
    Y <- from_B$Y
    linked_records <- cbind(from_A, Y)
    ATE_variance <- ATE(linked_records, 'propensity stratification')
    ATE_links <- ATE_variance[1]
    variance_links <- ATE_variance[2]
    ATE_list <- append(ATE_list, ATE_links)
    variance_list <- append(variance_list, variance_links)}
  return(c(ATE_list, variance_list))
}

### GENERATE ASSOCIATION AND LAUNCH ALGORITHM, BUILD CONFIDENCE INTERVALS ###

nbr_iter = 5 # more than 4 minutes per iteration for the following

name_DF <- 'DF_N=4401_2023-01-16.csv'

DF <- read.csv(file.path('..','..','datasets', name_DF), header=TRUE)
DF <- DF[!duplicated(DF),] # remove duplicates
# remove NaN values or empty elements:
DF <- DF[DF$name!="",]
DF <- DF[DF$family_name!="",]
DF <- DF[DF$was_born_female!="",]
DF <- DF[DF$country!="",]
DF <- DF[DF$birth_year!="",]
DF <- na.omit(DF)

identifiers <- list("family_name" = "jaro-winkler",
                    "was_born_female" = "strict",
                    "country" = "strict",
                    "birth_year" = "large")

covariates <- c('X1', 'X2', 'X3', 'X4', 'X5')

columns <- c("set","iter","variance","ate") 
results <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(results) <- columns

for(i in 1:nbr_iter) {
    
  ### GENERATE ASSOCIATION ### 
  
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
  
  ### fixed treatement effect:
  DF$Y <- - 10 + a*DF$treatment + b*exp(DF$X4) + c*DF$X3*DF$X1 + d*DF$X5 
  ### non-fixed treatement effect:
  # DF$Y <- - 10 + a*DF$treatment*DF$X2 + b*exp(DF$X4) + c*DF$X3*DF$X1 + d*DF$X5
  
  overlap <- DF[sample(nrow(DF),200),]
  
  A <- rbind(DF[sample(nrow(DF),350),], overlap)
  A <- A[,!names(A) %in% c("Y")]
  rownames(A) <- NULL
  
  B <- rbind(DF[sample(nrow(DF),500),], overlap)[,append(names(identifiers),"Y")]
  rownames(B) <- NULL
  
  ### MEV ALGORITHM ###
  
  # Parameters
  # match:    probability of having same linking variable when being true matches
  # unmatch:  probability of having same linking variable (at all)
  match <- rep(c(0.95), times=length(identifiers))
  unmatch <- rowSums( as.array( apply( A, 1, function(row) colSums(comparison_vector(row, B, identifiers)) ) ) ) / (nrow(A)*nrow(B))
  
  A$propensity_score <- propensity_score(A, covariates, minmaxscaler, TRUE)
  q <- 5
  A$prop_score_quantile <- with(A, cut(propensity_score, breaks = qu <- quantile(propensity_score, probs = seq(0,1, by=1/q)), labels = names(qu)[-1], include.lowest=TRUE))
  
  MEV <- Minimum_Estimated_Variance(A, B, identifiers, match, unmatch)
  ub_ate <- length(MEV)/2
  lb_var <- length(MEV)/2 + 1
  ATE_list <- MEV[1:ub_ate]
  variance_list <- MEV[lb_var:length(MEV)]
  
  len <- length(ATE_list)
  new_results <- data.frame(set=as.array(1:len), iter=rep(c(i), times=len), variance=variance_list, ate=ATE_list)
  results <- rbind(results, new_results)
}

darkorange <- rgb(0.804, 0.4, 0)
airforceblue <- rgb(0.36, 0.54, 0.66)

ggplot(data=results, aes(x=set, y=ate)) + stat_summary(geom = "line", fun = mean, color = airforceblue) + stat_summary(geom = "ribbon", fun.data = mean_cl_normal, alpha = 0.3) + geom_hline(yintercept=a, color = darkorange)

ggplot(data=results, aes(x=set, y=variance)) + stat_summary(geom = "line", fun = mean, color = airforceblue) + stat_summary(geom = "ribbon", fun.data = mean_cl_normal, alpha = 0.3) 
