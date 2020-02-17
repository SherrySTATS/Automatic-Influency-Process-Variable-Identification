##------------Author: Xinran Shi-------------------------------------##
##------------Date:2019-08-27----------------------------------------##
#setwd("../..") #as we all have different root folder in our own pc, this method won't fit all of us
#mypath = getwd()
mypath = "/Users/xinranshi/Dropbox/Share\ folders/Andi/ShuhData/Code"
# load CCA result
setwd(mypath)
setwd('./20_Models/21_CCA/Result')
#pval_cca = read.csv("pval_cca.csv",header = F, na.strings = "NA")
pval_cca = read.csv("pval_cca_seams.csv",header = F, na.strings = "NA") #seams
colnames(pval_cca) = "pval_cca"
pval_cca = pval_cca[-c(1,2),]
#load CART result
setwd(mypath)
setwd('./20_Models/22_CART/Result')
dir()
pval_tree = read.csv("pval_cart_20190917.csv",header = F, na.strings = "NA")
#pval_tree = read.csv("pval_cart_20191210_seams.csv",header = F, na.strings = "NA") #seams
colnames(pval_tree) = "pval_tree"
mcc_tree = read.csv("mcc_CART_20190917.csv",header = F,sep = " ")
#pval_tree = pval_tree[-c(1,2),]
#load Logistic regression result
setwd(mypath)
setwd('./20_Models/23_logistic_reg/Result')
dir()
pval_logit_reg = read.csv("pval_logistic_reg_20190910.csv",header = F, na.strings = "NA")
#pval_logit_reg = read.csv("pval_logistic_reg_20191210_seams.csv",header = F, na.strings = "NA") #seams
colnames(pval_logit_reg) = "pval_logit_reg"
#pval_logit_reg = pval_logit_reg[-c(1,2),]
mcc_logit = read.csv("mcc_logit_20190910.csv",header = F,sep = " ")

#load name
predName = c('stand5sidetemp','stand5bottomtemp','stand10flow',
             'stand16flow','stand21temp','stand21speed','pntm1waterboxflow','pntm2waterboxflow',
             'ntmentryspeed','ntmentrytemp','ntmexittemp','stand26speed','s26waterbox1',
             's26waterbox2','s26waterbox3','s26waterbox4')

ID = c(3:18)

#compute adj_pval
#P_matrix = rbind(pval_cca, t(pval_tree),t(pval_logit_reg))
#colnames(P_matrix) = predName
#adjpval = p.adjust(c(pval_cca, t(pval_tree), t(pval_logit_reg)), method = "holm")
#adjpval = p.adjust(P_matrix, method = "holm")
adj_cca = p.adjust(pval_cca, method = "holm")
adj_tree = p.adjust(t(pval_tree), method = "holm")
adj_logit = p.adjust(t(pval_logit_reg), method = "holm")
adjpval = rbind(adj_cca,adj_tree,adj_logit)
tadjpval = t(adjpval)

mcc = cbind.data.frame(t(mcc_tree),t(mcc_tree),t(mcc_logit))
names(mcc) = c("CCA","tree","logit")


library(nsprcomp)
temp = matrix(nrow = dim(mcc)[1], ncol = 1)
for (i in 1:dim(mcc)[1]) {
  test = cbind.data.frame((1-tadjpval[i,]),t(mcc[i,]))
  test2 = nsprcomp(test,nneg = T)
  temp[i] = max(test2$x[,1])
}
aggPIN = temp
rownames(aggPIN) = predName
plot(aggPIN)



#rank matrix
pval_rank = rbind(rank(adj_cca),rank(adj_tree),rank(adj_logit))
pval_rank_sum = colSums(pval_rank)

#adj_cca = adjpval[1:16]
#adj_tree = adjpval[17:32]
#adj_logit = adjpval[33:48]
#adj_data = cbind.data.frame(ID, predName, adj_logit, adj_cca, adj_tree)
adj_data = cbind.data.frame(ID, predName, aggPIN)


#visualization
m_adj_data = NA
library(reshape)
library(ggplot2)
library(ggrepel)

m_adj_data = melt(adj_data, id = c("ID","predName"))

p = ggplot(data=m_adj_data,
           aes(x=ID, y=value, colour=variable)) + geom_line() + geom_point(aes(shape = variable)) 
#p1 = p + geom_text_repel(data = subset(m_adj_data, value < 0.1), aes(label = predName,
#                                                                                fill = factor(variable)), color = 'black',
#                         size = 6,family = "Times New Roman",
#                         box.padding = 0.4) 

#p + geom_text_repel(data = subset(m_adj_data[17:32,], adj_cca < 0.05), aes(label = predName)) +  theme(text=element_text(size=20,  family="Times New Roman"))
p + geom_text_repel(data = subset(m_adj_data, adj_cca > 0.6), aes(label = predName)) +  theme(text=element_text(size=20,  family="Times New Roman"))

setwd(mypath)
setwd("../")
setwd("./Paper/plots/plot")
ggsave("aggPIN.png")
ggsave("20191106_adj_pval.png")

data = cbind.data.frame(ID, predName, pval_logit_reg, pval_cca, pval_tree)
mdata <- melt(data, id=c("ID","predName"))


p1 = ggplot(data, aes(x=ID, y=pval_logit_reg)) +
  geom_point() + geom_text_repel(data = subset(data, pval_logit_reg > 0), aes(label = predName))

p1

p2  = ggplot(data, aes(x=ID, y=pval_cca)) +
  geom_point() + geom_text_repel(data = subset(data, pval_cca > 0), aes(label = predName))
p2

p3 = ggplot(data, aes(x=ID, y=pval_tree)) +
  geom_point() + geom_text_repel(data = subset(data, pval_tree > 0), aes(label = predName))
p3


## result: 
## CCA < 0.05 : Stage  7 10 16
## Tree < 0.05: Stage 16
## Logit < 0.05: NA


#######result explore
inf_stage = NA
#inf_stage = c(which(adj_data$adj_logit<0.05), which(adj_data$adj_cca<0.05),which(adj_data$adj_tree<0.001)) + 2
#inf_stage = unique(inf_stage)
inf_stage = c(5,8,15)
setwd(mypath)
setwd("./30_P_val")
write.table(inf_stage,'./result/inf_stage_20200107.csv', row.names = F, col.names = F)
#predName[inf_stage - 2]
#"stand21temp"       "pntm2waterboxflow" "s26waterbox2" 
