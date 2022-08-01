#Phylogenetic regression 

library(phylolm)
library(phytools)
library(treedata.table)

#For first tree
tree_hap1_pr<-read.nexus("Datafile S1.nexus")
data_hap1_pr<-read.csv("Datafile S3.csv",header=TRUE,row.names = 1)
data_hap1_pr$tr1<-row.names(data_hap1_pr)
dataset_hap1<-as.treedata.table(tree_hap1_pr,data_hap1_pr)

#Running phylogenetic logistic regression using IG10 method
hap_glm_1a<-phyloglm(Eusociality~Haplodiploidy,data =data_hap1_pr,phy= tree_hap1_pr,method = c("logistic_MPLE"),boot=1000)
summary(hap_glm_1a)
hap_glm_1b<-phyloglm(Haplodiploidy~Eusociality,data=data_hap1_pr,phy=tree_hap1_pr,method = c("logistic_MPLE"),boot=1000)
summary(hap_glm_1b)

#For the second tree (where 11 species were dropped)
data_hap2_pr<-read.csv("Datafile S4.csv",header=TRUE,row.names = 1)
data_hap2_pr<-data_hap2_pr[!is.na(data_hap2_pr$Haplodiploidy),]
data_hap2_pr$tr2<-row.names(data_hap2_pr)
dataset_hap2<-as.treedata.table(tree_hap1_pr,data_hap2_pr)

#Running phylogenetic logistic regression using IG10 method

hap_glm_2a<-phyloglm(Eusociality~Haplodiploidy,data =data_hap2_pr,phy=dataset_hap2$phy,method = c("logistic_MPLE"))
summary(hap_glm_2a)
hap_glm_2b<-phyloglm(Haplodiploidy~Eusociality,data =data_hap2_pr,phy=dataset_hap2$phy,method = c("logistic_IG10"))
summary(hap_glm_2b)


#For the third tree
data_hap3_pr<-read.csv("Datafile S5.csv",header=TRUE,row.names = 1)
data_hap3_pr$tr3<-row.names(data_hap3_pr)
dataset_hap3<-as.treedata.table(tree_hap1_pr,data_hap3_pr)

#Running phylogenetic logistic regression using IG10 method
hap_glm_3a<-phyloglm(Eusociality~Haplodiploidy,data =data_hap3_pr,phy=dataset_hap3$phy,method = c("logistic_IG10"))
summary(hap_glm_3a)
hap_glm_3b<-phyloglm(Haplodiploidy~Eusociality,data =data_hap3_pr,phy=dataset_hap3$phy,method = c("logistic_IG10"))
summary(hap_glm_3b)

