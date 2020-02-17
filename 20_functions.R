library(InformationValue)
library(ROSE)
library(dplyr)
library(rpart)
library(BBmisc)
library(e1071)
#library(mltools)

#' @title remove_missing_levels
#' @description Accounts for missing factor levels present only in test data
#' but not in train data by setting values to NA
#'
#' @import magrittr
#' @importFrom gdata unmatrix
#' @importFrom stringr str_split
#'
#' @param fit fitted model on training data
#'
#' @param test_data data to make predictions for
#'
#' @return data.frame with matching factor levels to fitted model
#'
#' @keywords internal
#'
#' @export
remove_missing_levels <- function(fit, test_data) {
  
  # https://stackoverflow.com/a/39495480/4185785
  
  # drop empty factor levels in test data
  test_data %>%
    droplevels() %>%
    as.data.frame() -> test_data
  
  # 'fit' object structure of 'lm' and 'glmmPQL' is different so we need to
  # account for it
  if (any(class(fit) == "glmmPQL")) {
    # Obtain factor predictors in the model and their levels
    factors <- (gsub("[-^0-9]|as.factor|\\(|\\)", "",
                     names(unlist(fit$contrasts))))
    # do nothing if no factors are present
    if (length(factors) == 0) {
      return(test_data)
    }
    map(fit$contrasts, function(x) names(unmatrix(x))) %>%
      unlist() -> factor_levels
    factor_levels %>% str_split(":", simplify = TRUE) %>%
      extract(, 1) -> factor_levels
    model_factors <- as.data.frame(cbind(factors, factor_levels))
  } else {
    # Obtain factor predictors in the model and their levels
    factors <- (gsub("[-^0-9]|as.factor|\\(|\\)", "",
                     names(unlist(fit$xlevels))))
    # do nothing if no factors are present
    if (length(factors) == 0) {
      return(test_data)
    }
    
    factor_levels <- unname(unlist(fit$xlevels))
    model_factors <- as.data.frame(cbind(factors, factor_levels))
  }
  
  # Select column names in test data that are factor predictors in
  # trained model
  
  predictors <- names(test_data[names(test_data) %in% factors])
  
  # For each factor predictor in your data, if the level is not in the model,
  # set the value to NA
  
  for (i in 1:length(predictors)) {
    found <- test_data[, predictors[i]] %in% model_factors[
      model_factors$factors == predictors[i], ]$factor_levels
    if (any(!found)) {
      # track which variable
      var <- predictors[i]
      # set to NA
      test_data[!found, predictors[i]] <- NA
      # drop empty factor levels in test data
      test_data %>%
        droplevels() -> test_data
      # issue warning to console
      #message(sprintf(paste0("Setting missing levels in '%s', only present",
      #                      " in test data but missing in train data,",
      #                     " to 'NA'."),
      #             var))
    }
  }
  return(test_data)
}


mylogit_orign_mcc = function(mydata){
  mydata = as.data.frame(mydata)
  mydata$defect = factor(mydata$defect)
  train = cbind.data.frame(normalize(mydata[,-dim(mydata)[2]]),mydata$defect)
  names(train) = names(mydata)
  names(train)[dim(mydata)[2]] = "defect"
  n = dim(mydata)[2]
  ratio = table(mydata$defect)[1]/n
  #resample all data
  if(ratio>0.3|ratio<0.7){
    train.rose = train
  }
  else{train.rose <- ROSE(defect ~ ., data = train, seed = 123)$data}
  test = mydata
  mylogit <- glm(defect ~ ., data = train.rose, family = "binomial")
  newtest = remove_missing_levels (fit = mylogit, test_data=test)
  newtest = na.omit(newtest)
  glm.probs <- predict(mylogit,newdata=newtest, type = "response")
  #optCutOff <- optimalCutoff(newtest$defect, glm.probs, optimiseFor = "Zeros", returnDiagnostics = F)[1]
  optCutOff = 0.5
  glm.pred <- ifelse(glm.probs > optCutOff, 1, 0)
  n00 = as.numeric(sum(glm.pred == 0 &  newtest$defect == 0))
  n10 = as.numeric(sum(glm.pred == 1 &  newtest$defect == 0))
  n01 = as.numeric(sum(glm.pred == 0 &  newtest$defect == 1))
  n11 = as.numeric(sum(glm.pred == 1 &  newtest$defect == 1))
  num = (n00 * n11 - n10 * n01)
  denom = sqrt((n00 + n10)*(n00 + n01)*(n11 + n10)*(n11 + n01))
  #denom = sqrt(n00 + n10)/sqrt(n00 + n01)/sqrt(n11 + n10)/sqrt(n11 + n01)
  if(denom == 0){
    denom = 1
  }
  mcc = num / denom
  return(mcc)
}


mcc.boot = function(mydata,boot.times){
  k = sample(1:boot.times*100000, 1)
  set.seed(k)
  mydata$defect = factor(mydata$defect)
  mydf_boot = cbind.data.frame(normalize(mydata[,-dim(mydata)[2]]),mydata$defect)
  names(mydf_boot) = names(mydata)
  names(mydf_boot)[dim(mydata)[2]] = "defect"
  #mydf_boot = mydata
  n = dim(mydf_boot)[1]
  ratio = table(mydata$defect)[1]/n
  # = sum(as.numeric(as.character(mydata$defect)))/n
  Stage_permute1 = data.frame()
  Stage_permute1 = tbl_df(mydf_boot[,-dim(mydf_boot)[2]])
  mcc_boot = numeric(0)
  for (i in 1:boot.times){
    Stage_permute1  = sample_n(Stage_permute1, dim(Stage_permute1)[1])
    #Stage_permute1$defect_permute = rbinom(n, size=1, prob=ratio)
    #Stage_permute1 = sample(Stage_permute1)
    Stage_permute1$defect_permute = mydata$defect
    train <-Stage_permute1
    names(train)[-dim(train)[2]] = names(mydf_boot)[-dim(mydf_boot)[2]]
    if(ratio>0.3|ratio<0.7){
      Stage_permute = train
    }
    else{Stage_permute <- ROSE(defect_permute ~ ., data = train)$data}
    test <- mydf_boot
    mylogit_boot <- glm(defect_permute ~ ., data = Stage_permute, family = "binomial") #model is dynamic, mcc not uniform
    newtest = remove_missing_levels (fit = mylogit_boot, test_data=test)
    newtest = na.omit(newtest)
    glm.probs_boot <- predict(mylogit_boot,newdata = newtest)
    #optCutOff <- optimalCutoff(newtest$defect, glm.probs_boot)[1]
    optCutOff = 0.5
    glm.pred_boot <- ifelse(glm.probs_boot > optCutOff, 1, 0)
    n00 = sum(glm.pred_boot == 0 &  newtest$defect == 0)
    n10 = sum(glm.pred_boot == 1 &  newtest$defect == 0)
    n01 = sum(glm.pred_boot == 0 &  newtest$defect == 1)
    n11 = sum(glm.pred_boot == 1 &  newtest$defect == 1)
    num = (n00 * n11 - n10 * n01)
    denom = sqrt((n00 + n10)*(n00 + n01)*(n11 + n10)*(n11 + n01))
    #denom = sqrt(n00 + n10)/sqrt(n00 + n01)/sqrt(n11 + n10)/sqrt(n11 + n01)
    if(denom == 0){
      denom = 1
    }
    mcc_boot[i] = num / denom
    #mcc_boot[i] = (n00 * n11 - n10 * n01) / sqrt(n00 + n10)/sqrt(n00 + n01)/sqrt(n11 + n10)/sqrt(n11 + n01)
    #mcc_boot[i] = abs(mcc_boot[i])
  }
  mcc_boot_backup = mcc_boot
  return(mcc_boot_backup)
}

depidx = function(mcc, mcc.bootstrap){
  idx = (sum(mcc.bootstrap >= mcc, na.rm = T)+1 )/ (sum(!is.nan(mcc.bootstrap), na.rm = T)+1)
  if(is.na(mcc)==T){
    idx = NaN
  }
  return(idx)
}


##CART
mytree_orign_mcc = function(mydata){
  mydata = as.data.frame(mydata)
  mydata$defect = factor(mydata$defect)
  train = cbind.data.frame(normalize(mydata[,-dim(mydata)[2]]),mydata$defect)
  names(train) = names(mydata)
  names(train)[dim(mydata)[2]] = "defect"
  n = dim(mydata)[1]
  ratio = table(mydata$defect)[1]/n
  #train <- mydata
  #resample all data
  if(ratio>0.3|ratio<0.7){
    train.rose = train
  }
  else{train.rose <- ROSE(defect ~ ., data = train, seed = 123)$data}
  test <- mydata
  #grow the tree
  mytree <- rpart(defect ~ ., data = train.rose)
  # prune the tree
  pfit_tree <- prune(mytree, cp=mytree$cptable[which.min(mytree$cptable[,"xerror"]),"CP"]) #minimizes the cross-validated error
  newtest = remove_missing_levels (fit = pfit_tree, test_data=test)
  newtest = na.omit(newtest)
  tree.probs <- predict(pfit_tree,newdata=newtest,type = "vector")
  optCutOff <- optimalCutoff(newtest$defect, tree.probs)[1]
  tree.pred <- ifelse(tree.probs > optCutOff, "1", "0")
  n00 = sum(tree.pred == 0 &  newtest$defect == 0)
  n10 = sum(tree.pred == 1 &  newtest$defect == 0)
  n01 = sum(tree.pred == 0 &  newtest$defect == 1)
  n11 = sum(tree.pred == 1 &  newtest$defect == 1)
  #mcc = (n00 * n11 - n10 * n01) / sqrt(n00 + n10) / sqrt(n00 + n01) / sqrt(n11 + n10) / sqrt(n11 + n01)
  num = (n00 * n11 - n10 * n01)
  denom = sqrt((n00 + n10)*(n00 + n01)*(n11 + n10)*(n11 + n01))
  #denom = sqrt(n00 + n10)/sqrt(n00 + n01)/sqrt(n11 + n10)/sqrt(n11 + n01)
  if(denom == 0){
    denom = 1
  }
  mcc = num / denom
  return(mcc)
}

# bootstrap
mcc.boot.tree = function(mydata,boot.times){
  k = sample(1:boot.times*100000, 1)
  set.seed(k)
  mydata$defect = factor(mydata$defect)
  #mydf_boot = mydata
  mydf_boot = cbind.data.frame(normalize(mydata[,-dim(mydata)[2]]),mydata$defect)
  names(mydf_boot) = names(mydata)
  names(mydf_boot)[dim(mydata)[2]] = "defect"
  n = dim(mydf_boot)[1]
  ratio = sum(as.numeric(as.character(mydata$defect)))/n
  Stage_permute1 = data.frame()
  Stage_permute1 = tbl_df(mydf_boot[,-dim(mydf_boot)[2]])
  mcc_boot = numeric(0)
  for (i in 1:boot.times){
    Stage_permute1  = sample_n(Stage_permute1, dim(Stage_permute1)[1])
    Stage_permute1$defect_permute = mydata$defect
    train <-Stage_permute1
    #Stage_permute1$defect_permute = rbinom(n, size=1, prob=ratio)
    names(train)[-dim(train)[2]] = names(mydf_boot)[-dim(mydf_boot)[2]]
    if(ratio>0.3|ratio<0.7){
      Stage_permute = train
    }
    else{Stage_permute <- ROSE(defect_permute ~ ., data = train)$data}
    test <- mydf_boot
    mytree_boot <- rpart(defect_permute ~ ., data = Stage_permute) #model is dynamic, mcc not uniform
    pfit_tree_boot <- prune(mytree_boot, cp=mytree_boot$cptable[which.min(mytree_boot$cptable[,"xerror"]),"CP"]) #minimizes the cross-validated error
    newtest = remove_missing_levels (fit = pfit_tree_boot, test_data=test)
    newtest = na.omit(newtest)
    tree.probs_boot <- predict(pfit_tree_boot,newdata=newtest,type = "vector")
    optCutOff <- optimalCutoff(newtest$defect, tree.probs_boot)[1]
    tree.pred_boot <- ifelse(tree.probs_boot > optCutOff, 1, 0)
    n00 = sum(tree.pred_boot == 0 &  newtest$defect == 0)
    n10 = sum(tree.pred_boot == 1 &  newtest$defect == 0)
    n01 = sum(tree.pred_boot == 0 &  newtest$defect == 1)
    n11 = sum(tree.pred_boot == 1 &  newtest$defect == 1)
    num = (n00 * n11 - n10 * n01)
    #denom = sqrt(n00 + n10)/sqrt(n00 + n01)/sqrt(n11 + n10)/sqrt(n11 + n01)
    denom = sqrt((n00 + n10)*(n00 + n01)*(n11 + n10)*(n11 + n01))
    if(denom == 0){
      denom = 1
    }
    mcc_boot[i] = num / denom
    #mcc_boot[i] = (n00 * n11 - n10 * n01) / sqrt(n00 + n10) / sqrt(n00 + n01) / sqrt(n11 + n10) / sqrt(n11 + n01)
  }
  mcc_boot_backup = mcc_boot
  return(mcc_boot_backup)
}

##SVM
library(e1071)
mysvm_orign_mcc = function(mydata){
  mydata = as.data.frame(mydata)
  mydata$defect = factor(mydata$defect)
  train <- mydata
  #resample all data
  n = dim(mydata)[2]
  ratio = table(mydata$defect)[1]/n
  if(ratio>0.3|ratio<0.7){
    train.rose = train
  }
  else{train.rose <- ROSE(defect ~ ., data = train, seed = 123)$data}
  test <- mydata
  #fit svm
  mysvm <- svm(defect ~ ., data = train.rose, type = "C-classification", kernel = "radial")
  newtest = remove_missing_levels (fit = mysvm, test_data=test)
  newtest = na.omit(newtest)
  svm.pred <- predict(mysvm, newdata = newtest)
  n00 = sum(svm.pred == 0 &  newtest$defect == 0)
  n10 = sum(svm.pred == 1 &  newtest$defect == 0)
  n01 = sum(svm.pred == 0 &  newtest$defect == 1)
  n11 = sum(svm.pred == 1 &  newtest$defect == 1)
  #mcc = (n00 * n11 - n10 * n01) / sqrt(n00 + n10) / sqrt(n00 + n01) / sqrt(n11 + n10) / sqrt(n11 + n01)
  num = (n00 * n11 - n10 * n01)
  denom = sqrt((n00 + n10)*(n00 + n01)*(n11 + n10)*(n11 + n01))
  if(denom == 0){
    denom = 1
  }
  mcc = num / denom
  return(mcc)
}

# bootstrap
mcc.boot.svm = function(mydata,boot.times){
  #k = sample(1:boot.times*1000000, 1)
  #set.seed(k)
  set.seed(Sys.time())
  #mydata$defect = factor(mydata$defect)
  #mydf_boot = mydata
  mydf_boot = cbind.data.frame(normalize(mydata[,-dim(mydata)[2]]),mydata$defect)
  names(mydf_boot) = names(mydata)
  names(mydf_boot)[dim(mydata)[2]] = "defect"
  n = dim(mydf_boot)[1]
  Stage_permute1 = data.frame()
  Stage_permute1 = tbl_df(mydf_boot[,-dim(mydf_boot)[2]])
  mcc_boot = numeric(0)
  ratio = table(mydata$defect)[1]/n
  mcc_boot = numeric(0)
  for (i in 1:boot.times){
    Stage_permute1  = sample_n(Stage_permute1, dim(Stage_permute1)[1])
    Stage_permute1$defect_permute = mydata$defect
    #Stage_permute1$defect_permute = as.factor(rbinom(n, size=1, prob=ratio)) 
    #print(Stage_permute1$defect_permute)
    train <-Stage_permute1
    names(train)[-dim(train)[2]] = names(mydf_boot)[-dim(mydf_boot)[2]]
    if(ratio>0.3|ratio<0.7){
      Stage_permute = train
    }
    else{Stage_permute <- ROSE(defect_permute ~ ., data = train)$data}
    test <- mydf_boot
    mysvm_boot <- svm(defect_permute ~ ., data = Stage_permute, type = "C-classification", kernel = "radial") #model is dynamic, mcc not uniform
    newtest = remove_missing_levels (fit = mysvm_boot, test_data=test)
    newtest = na.omit(newtest)
    svm.pred_boot <- predict(mysvm_boot,newdata=newtest)
    n00 = sum(svm.pred_boot == 0 &  newtest$defect == 0)
    n10 = sum(svm.pred_boot == 1 &  newtest$defect == 0)
    n01 = sum(svm.pred_boot == 0 &  newtest$defect == 1)
    n11 = sum(svm.pred_boot == 1 &  newtest$defect == 1)
    num = (n00 * n11 - n10 * n01)
    denom = sqrt((n00 + n10)*(n00 + n01)*(n11 + n10)*(n11 + n01))
    if(denom == 0){
      denom = 1
    }
    mcc_boot[i] = num / denom
    #mcc_boot[i] = (n00 * n11 - n10 * n01) / sqrt(n00 + n10) / sqrt(n00 + n01) / sqrt(n11 + n10) / sqrt(n11 + n01)
  }
  mcc_boot_backup = mcc_boot
  return(mcc_boot_backup)
}

