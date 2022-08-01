library(phylolm)
library(phytools)
library(treedata.table)


tree_hap1_pr<-read.nexus("Datafile S1.nexus")
data_hap1_pr<-read.csv("Datafile S3_pc.csv",header=TRUE,row.names = 1)
 

#Running phylogenetic logistic regression using IG10 method
hap_glm_1a<-phyloglm(Eusociality~Haplodiploidy,data =data_hap1_pr,phy=tree_hap1_pr,method = c("logistic_IG10"),boot=1000)
summary(hap_glm_1a)


#Running phylogenetic logistic regression using MPLE method
hap_glm_1b<-phyloglm(Eusociality~Haplodiploidy,data =data_hap1_pr,phy=tree_hap1_pr,method = c("logistic_MPLE"),boot=1000)
summary(hap_glm_1b)



