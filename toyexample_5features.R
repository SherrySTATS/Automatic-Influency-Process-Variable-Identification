library(MASS)
library(parallel)
library(nsprcomp)
#generate data with collinearity 
n = 100
m = 5 # num of features
S_mean = c(0,0,0,0)
sig = matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4)
mydata = mvrnorm(n = 100, mu=S_mean, Sigma=sig, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
col_f = 2*mydata[,1]-1
mydf = cbind.data.frame(mydata,col_f)
names(mydf) = c("x1","x2","x3","x4","x5")
mystage = list()
mystage[[1]] = mydf$x1
mystage[[2]] = mydf$x2
mystage[[3]] = mydf$x3
mystage[[4]] = mydf$x4
mystage[[5]] = mydf$x5 #colinear with x1

#Case 1
#colinearity 
z = mydf$x1 + mydf$x5       # linear combination with a bias
pr = 1/(1+exp(-z))         # pass through an inv-logit function
y = rbinom(n,1,pr)      # bernoulli response variable
rp = y
#apend response to each stage
for (i in 1:5) {
  mystage[[i]] = cbind.data.frame(mystage[[i]],rp)
  names(mystage[[i]])[1]<-paste("x",i,sep = "")
  names(mystage[[i]])[length(names(mystage[[i]]))]<-"defect"
}

cl <- makeCluster(detectCores())
clusterEvalQ(cl, list(library(lme4),library(ROSE),library(dplyr),library(InformationValue),library(BBmisc)))
clusterExport(cl, "remove_missing_levels")
mcc.logit = parLapply(cl,mystage, mylogit_orign_mcc)
mcc.boot.backup.logit = parLapply(cl,mystage,boot.times = 10000, mcc.boot)
stopCluster(cl)
mcc.logit.abs = lapply(mcc.logit, abs)
mcc.boot.backup.logit.abs =  lapply(mcc.boot.backup.logit, abs)
depIDX.logit = mapply(depidx, mcc.logit.abs, mcc.bootstrap = mcc.boot.backup.logit.abs)


cl <- makeCluster(detectCores())
clusterEvalQ(cl, list(library(lme4),library(ROSE),library(dplyr),library(InformationValue),library(BBmisc),library(rpart)))
clusterExport(cl, "remove_missing_levels")
mcc.tree = parLapply(cl,mystage, mytree_orign_mcc)
mcc.boot.backup.tree = parLapply(cl,mystage,boot.times = 10000, mcc.boot.tree)
stopCluster(cl)

mcc.tree.abs = lapply(mcc.tree, abs)
mcc.boot.backup.tree.abs =  lapply(mcc.boot.backup.tree, abs)
depIDX.tree = mapply(depidx, mcc.tree.abs, mcc.bootstrap = mcc.boot.backup.tree.abs)

cl <- makeCluster(detectCores())
clusterEvalQ(cl, list(library(lme4),library(ROSE),library(dplyr),library(InformationValue),library(BBmisc),library(e1071)))
clusterExport(cl, "remove_missing_levels")
mcc.svm = parLapply(cl,mystage, mysvm_orign_mcc)
mcc.boot.backup.svm = parLapply(cl,mystage,boot.times = 10000, mcc.boot.svm)
stopCluster(cl)
mcc.svm.abs = lapply(mcc.svm, abs)
mcc.boot.backup.svm.abs =  lapply(mcc.boot.backup.svm, abs)
depIDX.svm = mapply(depidx, mcc.svm.abs, mcc.bootstrap = mcc.boot.backup.svm.abs)

dep.all = cbind(unlist(depIDX.logit),unlist(depIDX.tree),unlist(depIDX.svm))
mcc.all = cbind(unlist(mcc.logit.abs),unlist(mcc.tree.abs),unlist(mcc.svm.abs))
dep.all[which(is.na(dep.all))]=1
dep.all.adj = matrix(nrow = dim(dep.all)[1], ncol = 3)
for (i in 1:dim(dep.all)[2]) {
  dep.all.adj[,i] = p.adjust(dep.all[,i],method = "holm")
}

dep.all.adj[which(is.na(dep.all.adj))]=1
mcc.all[which(is.na(mcc.all),arr.ind = T)]=0
case1.score = apply(matrix(nsprcomp(cbind(c(1-dep.all.adj),c(mcc.all)),ncomp = 1,nneg = T)$x[,1],nrow = dim(dep.all)[1],ncol = dim(dep.all)[2]),1,max)
case1.score 
barplot(t(case1.score), col = "blue",names = c("x1","x2","x3","x4","x5"))

#Case 2
#interaction
z = mydf$x1 * mydf$x2       # interaction between x1 and x2
pr = 1/(1+exp(-z))         # pass through an inv-logit function
y = rbinom(n,1,pr)      # bernoulli response variable
rp = y
#apend response to each stage
for (i in 1:5) {
  mystage[[i]]$defect = rp
  names(mystage[[i]])[length(names(mystage[[i]]))]<-"defect"
}

cl <- makeCluster(detectCores())
clusterEvalQ(cl, list(library(lme4),library(ROSE),library(dplyr),library(InformationValue),library(BBmisc)))
clusterExport(cl, "remove_missing_levels")
mcc.logit = parLapply(cl,mystage, mylogit_orign_mcc)
stopCluster(cl)

cl <- makeCluster(detectCores())
clusterEvalQ(cl, list(library(lme4),library(ROSE),library(dplyr),library(InformationValue),library(BBmisc),library(rpart)))
clusterExport(cl, "remove_missing_levels")
mcc.boot.backup.logit = parLapply(cl,mystage,boot.times = 10000, mcc.boot)
stopCluster(cl)
mcc.logot.abs = lapply(mcc.logit, abs)
mcc.boot.backup.logit.abs =  lapply(mcc.boot.backup.logit, abs)
depIDX.logit = mapply(depidx, mcc.logit.abs, mcc.bootstrap = mcc.boot.backup.logit.abs)


cl <- makeCluster(detectCores())
clusterEvalQ(cl, list(library(lme4),library(ROSE),library(dplyr),library(InformationValue),library(BBmisc),library(rpart)))
clusterExport(cl, "remove_missing_levels")
mcc.tree = parLapply(cl,mystage, mytree_orign_mcc)
stopCluster(cl)

cl <- makeCluster(detectCores())
clusterEvalQ(cl, list(library(lme4),library(ROSE),library(dplyr),library(InformationValue),library(BBmisc),library(rpart)))
clusterExport(cl, "remove_missing_levels")
mcc.boot.backup.tree = parLapply(cl,mystage,boot.times = 10000, mcc.boot.tree)
stopCluster(cl)

mcc.tree.abs = lapply(mcc.tree, abs)
mcc.boot.backup.tree.abs =  lapply(mcc.boot.backup.tree, abs)
depIDX.tree = mapply(depidx, mcc.tree.abs, mcc.bootstrap = mcc.boot.backup.tree.abs)

cl <- makeCluster(detectCores())
clusterEvalQ(cl, list(library(lme4),library(ROSE),library(dplyr),library(InformationValue),library(BBmisc),library(e1071)))
clusterExport(cl, "remove_missing_levels")
mcc.svm = parLapply(cl,mystage, mysvm_orign_mcc)
stopCluster(cl)

cl <- makeCluster(detectCores())
clusterEvalQ(cl, list(library(lme4),library(ROSE),library(dplyr),library(InformationValue),library(BBmisc),library(e1071)))
clusterExport(cl, "remove_missing_levels")
mcc.boot.backup.svm = parLapply(cl,mystage,boot.times = 10000, mcc.boot.svm)
stopCluster(cl)
mcc.svm.abs = lapply(mcc.svm, abs)
mcc.boot.backup.svm.abs =  lapply(mcc.boot.backup.svm, abs)
depIDX.svm = mapply(depidx, mcc.svm.abs, mcc.bootstrap = mcc.boot.backup.svm.abs)

dep.all = cbind(unlist(depIDX.logit),unlist(depIDX.tree),unlist(depIDX.svm))
mcc.all = cbind(unlist(mcc.logit.abs),unlist(mcc.tree.abs),unlist(mcc.svm.abs))
dep.all[which(is.na(dep.all))]=1
dep.all.adj = matrix(nrow = dim(dep.all)[1], ncol = 3)
for (i in 1:dim(dep.all)[2]) {
  dep.all.adj[,i] = p.adjust(dep.all[,i],method = "holm")
}

dep.all.adj[which(is.na(dep.all.adj))]=1
mcc.all[which(is.na(mcc.all),arr.ind = T)]=0
case2.score = apply(matrix(nsprcomp(cbind(c(1-dep.all.adj),c(mcc.all)),ncomp = 1,nneg = T)$x[,1],nrow = dim(dep.all)[1],ncol = dim(dep.all)[2]),1,max)
case2.score 
barplot(t(case2.score), col = "blue",names = c("x1","x2","x3","x4","x5"))


#Case 3
#simple case
z = mydf$x2       #only x2 is informative
pr = 1/(1+exp(-z))         # pass through an inv-logit function
y = rbinom(n,1,pr)      # bernoulli response variable
rp = y
#apend response to each stage
for (i in 1:5) {
  mystage[[i]]$defect = rp
  names(mystage[[i]])[length(names(mystage[[i]]))]<-"defect"
}

cl <- makeCluster(detectCores())
clusterEvalQ(cl, list(library(lme4),library(ROSE),library(dplyr),library(InformationValue),library(BBmisc)))
clusterExport(cl, "remove_missing_levels")
mcc.logit = parLapply(cl,mystage, mylogit_orign_mcc)
mcc.boot.backup.logit = parLapply(cl,mystage,boot.times = 10000, mcc.boot)
stopCluster(cl)
mcc.logit.abs = lapply(mcc.logit, abs)
mcc.boot.backup.logit.abs =  lapply(mcc.boot.backup.logit, abs)
depIDX.logit = mapply(depidx, mcc.logit.abs, mcc.bootstrap = mcc.boot.backup.logit.abs)


cl <- makeCluster(detectCores())
clusterEvalQ(cl, list(library(lme4),library(ROSE),library(dplyr),library(InformationValue),library(BBmisc),library(rpart)))
clusterExport(cl, "remove_missing_levels")
mcc.tree = parLapply(cl,mystage, mytree_orign_mcc)
mcc.boot.backup.tree = parLapply(cl,mystage,boot.times = 10000, mcc.boot.tree)
stopCluster(cl)

mcc.tree.abs = lapply(mcc.tree, abs)
mcc.boot.backup.tree.abs =  lapply(mcc.boot.backup.tree, abs)
depIDX.tree = mapply(depidx, mcc.tree.abs, mcc.bootstrap = mcc.boot.backup.tree.abs)

cl <- makeCluster(detectCores())
clusterEvalQ(cl, list(library(lme4),library(ROSE),library(dplyr),library(InformationValue),library(BBmisc),library(e1071)))
clusterExport(cl, "remove_missing_levels")
mcc.svm = parLapply(cl,mystage, mysvm_orign_mcc)
mcc.boot.backup.svm = parLapply(cl,mystage,boot.times = 10000, mcc.boot.svm)
stopCluster(cl)
mcc.svm.abs = lapply(mcc.svm, abs)
mcc.boot.backup.svm.abs =  lapply(mcc.boot.backup.svm, abs)
depIDX.svm = mapply(depidx, mcc.svm.abs, mcc.bootstrap = mcc.boot.backup.svm.abs)

dep.all = cbind(unlist(depIDX.logit),unlist(depIDX.tree),unlist(depIDX.svm))
mcc.all = cbind(unlist(mcc.logit.abs),unlist(mcc.tree.abs),unlist(mcc.svm.abs))
dep.all[which(is.na(dep.all))]=1
dep.all.adj = matrix(nrow = dim(dep.all)[1], ncol = 3)
for (i in 1:dim(dep.all)[2]) {
  dep.all.adj[,i] = p.adjust(dep.all[,i],method = "holm")
}

dep.all.adj[which(is.na(dep.all.adj))]=1
mcc.all[which(is.na(mcc.all),arr.ind = T)]=0
case3.score = apply(matrix(nsprcomp(cbind(c(1-dep.all.adj),c(mcc.all)),ncomp = 1,nneg = T)$x[,1],nrow = dim(dep.all)[1],ncol = dim(dep.all)[2]),1,max)
case3.score 
barplot(t(case3.score), col = "blue",names = c("x1","x2","x3","x4","x5"))

#Case 4
#more complicated case
z = sin(mydf$x2)       #only x2 is informative
pr = 1/(1+exp(-z))         # pass through an inv-logit function
y = rbinom(n,1,pr)      # bernoulli response variable
rp = y
#apend response to each stage
for (i in 1:5) {
  mystage[[i]]$defect = rp
  names(mystage[[i]])[length(names(mystage[[i]]))]<-"defect"
}

cl <- makeCluster(detectCores())
clusterEvalQ(cl, list(library(lme4),library(ROSE),library(dplyr),library(InformationValue),library(BBmisc)))
clusterExport(cl, "remove_missing_levels")
mcc.logit = parLapply(cl,mystage, mylogit_orign_mcc)
mcc.boot.backup.logit = parLapply(cl,mystage,boot.times = 10000, mcc.boot)
stopCluster(cl)
mcc.logit.abs = lapply(mcc.logit, abs)
mcc.boot.backup.logit.abs =  lapply(mcc.boot.backup.logit, abs)
depIDX.logit = mapply(depidx, mcc.logit.abs, mcc.bootstrap = mcc.boot.backup.logit.abs)


cl <- makeCluster(detectCores())
clusterEvalQ(cl, list(library(lme4),library(ROSE),library(dplyr),library(InformationValue),library(BBmisc),library(rpart)))
clusterExport(cl, "remove_missing_levels")
mcc.tree = parLapply(cl,mystage, mytree_orign_mcc)
mcc.boot.backup.tree = parLapply(cl,mystage,boot.times = 10000, mcc.boot.tree)
stopCluster(cl)

mcc.tree.abs = lapply(mcc.tree, abs)
mcc.boot.backup.tree.abs =  lapply(mcc.boot.backup.tree, abs)
depIDX.tree = mapply(depidx, mcc.tree.abs, mcc.bootstrap = mcc.boot.backup.tree.abs)

cl <- makeCluster(detectCores())
clusterEvalQ(cl, list(library(lme4),library(ROSE),library(dplyr),library(InformationValue),library(BBmisc),library(e1071)))
clusterExport(cl, "remove_missing_levels")
mcc.svm = parLapply(cl,mystage, mysvm_orign_mcc)
mcc.boot.backup.svm = parLapply(cl,mystage,boot.times = 10000, mcc.boot.svm)
stopCluster(cl)
mcc.svm.abs = lapply(mcc.svm, abs)
mcc.boot.backup.svm.abs =  lapply(mcc.boot.backup.svm, abs)
depIDX.svm = mapply(depidx, mcc.svm.abs, mcc.bootstrap = mcc.boot.backup.svm.abs)

dep.all = cbind(unlist(depIDX.logit),unlist(depIDX.tree),unlist(depIDX.svm))
mcc.all = cbind(unlist(mcc.logit.abs),unlist(mcc.tree.abs),unlist(mcc.svm.abs))
dep.all[which(is.na(dep.all))]=1
dep.all.adj = matrix(nrow = dim(dep.all)[1], ncol = 3)
for (i in 1:dim(dep.all)[2]) {
  dep.all.adj[,i] = p.adjust(dep.all[,i],method = "holm")
}

dep.all.adj[which(is.na(dep.all.adj))]=1
mcc.all[which(is.na(mcc.all),arr.ind = T)]=0
case4.score = apply(matrix(nsprcomp(cbind(c(1-dep.all.adj),c(mcc.all)),ncomp = 1,nneg = T)$x[,1],nrow = dim(dep.all)[1],ncol = dim(dep.all)[2]),1,max)
case4.score 
barplot(t(case4.score), col = "blue",names = c("x1","x2","x3","x4","x5"))
#cor_mat = as.matrix(cbind(mydf,z))
#corrplot::corrplot(cor(cor_mat))


#Case 5
#more complicated case
z = mydf$x3*sin(mydf$x2)       #only x2 is informative
pr = 1/(1+exp(-z))         # pass through an inv-logit function
y = rbinom(n,1,pr)      # bernoulli response variable
rp = y
#apend response to each stage
for (i in 1:5) {
  mystage[[i]]$defect = rp
  names(mystage[[i]])[length(names(mystage[[i]]))]<-"defect"
}
mystage[[6]] = list()
mystage[[6]] = cbind.data.frame(z,rp)
names(mystage[[6]]) = c("z","defect")
cl <- makeCluster(detectCores())
clusterEvalQ(cl, list(library(lme4),library(ROSE),library(dplyr),library(InformationValue),library(BBmisc)))
clusterExport(cl, "remove_missing_levels")
mcc.logit = parLapply(cl,mystage, mylogit_orign_mcc)
mcc.boot.backup.logit = parLapply(cl,mystage,boot.times = 50000, mcc.boot)
stopCluster(cl)
mcc.logit.abs = lapply(mcc.logit, abs)
mcc.boot.backup.logit.abs =  lapply(mcc.boot.backup.logit, abs)
depIDX.logit = mapply(depidx, mcc.logit.abs, mcc.bootstrap = mcc.boot.backup.logit.abs)


cl <- makeCluster(detectCores())
clusterEvalQ(cl, list(library(lme4),library(ROSE),library(dplyr),library(InformationValue),library(BBmisc),library(rpart)))
clusterExport(cl, "remove_missing_levels")
mcc.tree = parLapply(cl,mystage, mytree_orign_mcc)
mcc.boot.backup.tree = parLapply(cl,mystage,boot.times = 50000, mcc.boot.tree)
stopCluster(cl)

mcc.tree.abs = lapply(mcc.tree, abs)
mcc.boot.backup.tree.abs =  lapply(mcc.boot.backup.tree, abs)
depIDX.tree = mapply(depidx, mcc.tree.abs, mcc.bootstrap = mcc.boot.backup.tree.abs)

cl <- makeCluster(detectCores())
clusterEvalQ(cl, list(library(lme4),library(ROSE),library(dplyr),library(InformationValue),library(BBmisc),library(e1071)))
clusterExport(cl, "remove_missing_levels")
mcc.svm = parLapply(cl,mystage, mysvm_orign_mcc)
mcc.boot.backup.svm = parLapply(cl,mystage,boot.times = 50000, mcc.boot.svm)
stopCluster(cl)
mcc.svm.abs = lapply(mcc.svm, abs)
mcc.boot.backup.svm.abs =  lapply(mcc.boot.backup.svm, abs)
depIDX.svm = mapply(depidx, mcc.svm.abs, mcc.bootstrap = mcc.boot.backup.svm.abs)

dep.all = cbind(unlist(depIDX.logit),unlist(depIDX.tree),unlist(depIDX.svm))
#mcc.all = cbind(unlist(mcc.logit),unlist(mcc.tree),unlist(mcc.svm))
mcc.all = cbind(unlist(mcc.logit.abs),unlist(mcc.tree.abs),unlist(mcc.svm.abs))
dep.all[which(is.na(dep.all))]=1
dep.all.adj = matrix(nrow = dim(dep.all)[1], ncol = 3)
for (i in 1:dim(dep.all)[2]) {
  dep.all.adj[,i] = p.adjust(dep.all[,i],method = "holm")
}

dep.all.adj[which(is.na(dep.all.adj))]=1
mcc.all[which(is.na(mcc.all),arr.ind = T)]=0
case5.score = apply(matrix(nsprcomp(cbind(c(1-dep.all.adj),c(mcc.all)),ncomp = 1,nneg = T)$x[,1],nrow = dim(dep.all)[1],ncol = dim(dep.all)[2]),1,max)
case5.score 
barplot(t(case5.score), col = "blue",names = c("x1","x2","x3","x4","x5","sin(x2)x3"))

