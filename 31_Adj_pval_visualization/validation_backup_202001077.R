#validation
mypath =  "/Users/xinranshi/Dropbox/Share\ folders/Andi/ShuhData/Code"
setwd(mypath)
setwd("../")
setwd("./Paper/plots/plot")
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(magrittr)
predName = c('NA','NA','stand5sidetemp','stand5bottomtemp','stand10flow',
             'stand16flow','stand21temp','stand21speed','pntm1waterboxflow','pntm2waterboxflow',
             'ntmentryspeed','ntmentrytemp','ntmexittemp','stand26speed','s26waterbox1',
             's26waterbox2','s26waterbox3','s26waterbox4')

#load data
setwd(mypath)
setwd('./10_Data-processing/10_Generate_Features/data')
#load feature name
library(stringr)
setwd(mypath)
setwd("./10_Data-processing/10_Generate_Features/data")
FeatureName = NA
FeatureName = readLines("FeatureName.txt",n=18)
FeatureName = c(NA,NA,FeatureName)
FeatureName = FeatureName[-3] #Feature name of each sensor
sub_feature_name = list()
for (i in 3:18) {
  sub_feature_name[[i]] = str_split(FeatureName[i],",")
}
#read data
mystage = list()
for (i in 3:18) {
  fileID = paste('Stage',i,'.txt',sep = "")
  mystage[[i]] = read.table(fileID, sep = ",")
  names(mystage[[i]]) = sub_feature_name[[i]][[1]]
  mymatch = match(sub_feature_name[[i]][[1]],"Group_label")
  mymatch[is.na(mymatch)] = 0
  if(mymatch == 1)
  {
    mystage[[i]]$Group_label = as.factor(mystage[[i]]$Group_label)
  }
}
rp = read.table("response.txt", sep = ",")
rp = as.factor(rp$V1)
#apend response to each stage
for (i in 3:18) {
  mystage[[i]] = cbind(mystage[[i]],rp)
  names( mystage[[i]])[length(names(mystage[[i]]))]<-"defect"
}

setwd(mypath)
setwd("./30_P_val/result")
inf_stage = read.csv("inf_stage_20200107.csv",header = F)
inf_stage = data.matrix(inf_stage)
inf_stage = c(5,8,15)
#k-s test
#sensor 5
myks.p = NA
myks.d = NA
test = mystage[[5]] %>% group_split(defect)
myks.p[1] = ks.test(test[[1]]$PCA_score1,test[[2]]$PCA_score1)$p.value
myks.p[2] = ks.test(test[[1]]$PCA_score2,test[[2]]$PCA_score2)$p.value
myks.p[3] = ks.test(test[[1]]$Mean,test[[2]]$Mean)$p.value
myks.p[4] = ks.test(test[[1]]$Total_variation,test[[2]]$Total_variation)$p.value
myks.p[5] = ks.test(test[[1]]$Curvature,test[[2]]$Curvature)$p.value

myks.d[1] = ks.test(test[[1]]$PCA_score1,test[[2]]$PCA_score1)$statistic
myks.d[2] = ks.test(test[[1]]$PCA_score2,test[[2]]$PCA_score2)$statistic
myks.d[3] = ks.test(test[[1]]$Mean,test[[2]]$Mean)$statistic
myks.d[4] = ks.test(test[[1]]$Total_variation,test[[2]]$Total_variation)$statistic
myks.d[5] = ks.test(test[[1]]$Curvature,test[[2]]$Curvature)$statistic
myks = cbind(myks.d,myks.p)

mymean = NA
temp1 = c(mean(test[[1]]$PCA_score1),mean(test[[2]]$PCA_score1))
temp2 = c(mean(test[[1]]$PCA_score2),mean(test[[2]]$PCA_score2))
temp3 = c(mean(test[[1]]$Mean),mean(test[[2]]$Mean))
temp4 = c(mean(test[[1]]$Total_variation),mean(test[[2]]$Total_variation))
temp5 = c(mean(test[[1]]$Curvature),mean(test[[2]]$Curvature))
mymean = rbind(t(temp1),t(temp2),t(temp3),t(temp4),t(temp5))

mysd = NA
temp1 = c(sd(test[[1]]$PCA_score1),sd(test[[2]]$PCA_score1))
temp2 = c(sd(test[[1]]$PCA_score2),sd(test[[2]]$PCA_score2))
temp3 = c(sd(test[[1]]$Mean),sd(test[[2]]$Mean))
temp4 = c(sd(test[[1]]$Total_variation),sd(test[[2]]$Total_variation))
temp5 = c(sd(test[[1]]$Curvature),sd(test[[2]]$Curvature))
mysd = rbind(t(temp1),t(temp2),t(temp3),t(temp4),t(temp5))


F_5 = cbind(mymean,mysd,myks)
colnames(F_5)  = c("mean_non_def","mean_def","std_non_def","std_def","ks.d","ks.p")


#sensor 8
myks_8 = NA
myks_8.p = NA
myks_8.d = NA
test = mystage[[8]] %>% group_split(defect)
myks_8.p[1] = ks.test(test[[1]]$Mode,test[[2]]$Mode)$p.value
myks_8.p[2] = ks.test(test[[1]]$Total_variation,test[[2]]$Total_variation)$p.value
myks_8.p[3] = ks.test(test[[1]]$Range,test[[2]]$Range)$p.value
myks_8.d[1] = ks.test(test[[1]]$Mode,test[[2]]$Mode)$statistic
myks_8.d[2] = ks.test(test[[1]]$Total_variation,test[[2]]$Total_variation)$statistic
myks_8.d[3] = ks.test(test[[1]]$Range,test[[2]]$Range)$statistic
myks_8 = cbind(myks_8.d,myks_8.p)

mymean_8 = NA
temp1 = c(mean(test[[1]]$Mode),mean(test[[2]]$Mode))
temp2 = c(mean(test[[1]]$Total_variation),mean(test[[2]]$Total_variation))
temp3 = c(mean(test[[1]]$Range),mean(test[[2]]$Range))
mymean_8 = rbind(t(temp1),t(temp2),t(temp3))

mysd_8 = NA
temp1 = c(sd(test[[1]]$Mode),sd(test[[2]]$Mode))
temp2 = c(sd(test[[1]]$Total_variation),sd(test[[2]]$Total_variation))
temp3 = c(sd(test[[1]]$Range),sd(test[[2]]$Range))
mysd_8 = rbind(t(temp1),t(temp2),t(temp3))

F_8 = cbind(mymean_8,mysd_8,myks_8)
colnames(F_8)  = c("mean_non_def","mean_def","std_non_def","std_def","ks.d","ks.p")

#sensor 15
myks_15 = NA
myks_15.p = NA
myks_15.d = NA
test = mystage[[15]] %>% group_split(defect)
myks_15.p[1] = ks.test(test[[1]]$Total_variation,test[[2]]$Total_variation)$p.value
myks_15.p[2] = ks.test(test[[1]]$Range,test[[2]]$Range)$p.value
myks_15.p[3] = ks.test(test[[1]]$Curvature,test[[2]]$Curvature)$p.value
myks_15.p[4] = ks.test(test[[1]]$Median,test[[2]]$Median)$p.value
myks_15.p[5] = ks.test(test[[1]]$Slope,test[[2]]$Slope)$p.value
myks_15.p[6] = ks.test(test[[1]]$Fitting_error_of_Huber_regression_model,test[[2]]$Fitting_error_of_Huber_regression_model)$p.value
myks_15.d[1] = ks.test(test[[1]]$Total_variation,test[[2]]$Total_variation)$statistic
myks_15.d[2] = ks.test(test[[1]]$Range,test[[2]]$Range)$statistic
myks_15.d[3] = ks.test(test[[1]]$Curvature,test[[2]]$Curvature)$statistic
myks_15.d[4] = ks.test(test[[1]]$Median,test[[2]]$Median)$statistic
myks_15.d[5] = ks.test(test[[1]]$Slope,test[[2]]$Slope)$statistic
myks_15.d[6] = ks.test(test[[1]]$Fitting_error_of_Huber_regression_model,test[[2]]$Fitting_error_of_Huber_regression_model)$statistic
myks_15 = cbind(myks_15.d,myks_15.p)

mymean_15 = NA
temp1 = c(mean(test[[1]]$Total_variation),mean(test[[2]]$Total_variation))
temp2 = c(mean(test[[1]]$Range),mean(test[[2]]$Range))
temp3 = c(mean(test[[1]]$Curvature),mean(test[[2]]$Curvature))
temp4 = c(mean(test[[1]]$Median),mean(test[[2]]$Median))
temp5 = c(mean(test[[1]]$Slope),mean(test[[2]]$Slope))
temp6 = c(mean(test[[1]]$Fitting_error_of_Huber_regression_model),mean(test[[2]]$Fitting_error_of_Huber_regression_model))
mymean_15 = rbind(t(temp1),t(temp2),t(temp3),t(temp4),t(temp5),t(temp6))

mysd_15 = NA
temp1 = c(sd(test[[1]]$Total_variation),sd(test[[2]]$Total_variation))
temp2 = c(sd(test[[1]]$Range),sd(test[[2]]$Range))
temp3 = c(sd(test[[1]]$Curvature),sd(test[[2]]$Curvature))
temp4 = c(sd(test[[1]]$Median),sd(test[[2]]$Median))
temp5 = c(sd(test[[1]]$Slope),sd(test[[2]]$Slope))
temp6 = c(sd(test[[1]]$Fitting_error_of_Huber_regression_model),sd(test[[2]]$Fitting_error_of_Huber_regression_model))
mysd_15 = rbind(t(temp1),t(temp2),t(temp3),t(temp4),t(temp5),t(temp6))

F_15 = cbind(mymean_15,mysd_15,myks_15)
colnames(F_15)  = c("mean_non_def","mean_def","std_non_def","std_def","ks.d","ks.p")

#visalization
for (i in 1:length(inf_stage)) {
  temp = as.data.frame(mystage[[inf_stage[i]]])
  #temp = temp[,-1]
  temp$defect = as.factor(temp$defect)
  temp1 = inf_stage[i]
  inf_feature = dim(temp)[2]-1
  temp$ID = 1:dim(temp)[1]
  plot_name = paste("Checking_cdf_",predName[inf_stage[i]],"_",sub_feature_name[[inf_stage[i]]][[1]][1:inf_feature],sep = "")
    for (j in 1:inf_feature) {
      mymatch = match(sub_feature_name[[inf_stage[i]]][[1]][j],"Group_label")
      mymatch[is.na(mymatch)] = 0
      if(mymatch == 0)
      { 
        myplot_name = plot_name[j]
        mtemp = melt(temp, id = c("ID","defect"))
        a = (j-1)*dim(temp)[1]+1
        b = j*dim(temp)[1]
        mytemp = mtemp[a:b,]
        mytemp$defect = as.factor(mytemp$defect)
        mytemp$value = as.numeric(mytemp$value)
        compare_mean = mytemp %>% group_by(defect) %>% summarise(Mean = mean(value))
        p = ggplot(mytemp, aes(x = value, colour = defect, linetype = defect)) +
          stat_ecdf(geom = "step",size = 2) + xlab(names(temp)[j]) 
        #p = p + geom_vline(data = compare_mean, aes(xintercept = Mean, color = defect),
        #             linetype = "dashed", size = 1) + geom_density(alpha = .5)
        p + theme(text=element_text(size=24,  family="Times New Roman"))
        ggsave(paste(myplot_name,".png",sep = ""))
    }
  else
  {
    myplot_name = plot_name[j]
    temp$Group_label = as.factor(temp$Group_label)
    p = ggplot(temp, aes(x = defect, fill = Group_label)) +
      stat_ecdf(geom = "step",size = 2)
    p + theme(text=element_text(size=24,  family="Times New Roman"))
    ggsave(paste(myplot_name,".png",sep = ""))
  }
}
}

#sensor 16
myks_16 = NA
myks_16.p = NA
myks_16.d = NA
test = mystage[[16]] %>% group_split(defect)
myks_16.p[1] = ks.test(test[[1]]$Mean,test[[2]]$Mean)$p.value
myks_16.p[2] = ks.test(test[[1]]$Total_variation,test[[2]]$Total_variation)$p.value
myks_16.p[3] = ks.test(test[[1]]$Curvature,test[[2]]$Curvature)$p.value
myks_16.d[1] = ks.test(test[[1]]$Mean,test[[2]]$Mean)$statistic
myks_16.d[2] = ks.test(test[[1]]$Total_variation,test[[2]]$Total_variation)$statistic
myks_16.d[3] = ks.test(test[[1]]$Curvature,test[[2]]$Curvature)$statistic
myks_16 = cbind(myks_16.d,myks_16.p)




