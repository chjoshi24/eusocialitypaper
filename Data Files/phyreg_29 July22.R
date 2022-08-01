library(phylolm)
library(phytools)
library(treedata.table)

#FOR THE ORIGINAL DATASET 

#Doing phylogenetic regression without any manipulations in the data
tree_hap1<-read.nexus("Datafile S1.nexus")
data_hap1<-read.csv("Datafile S3.csv",header=TRUE,row.names = 1)

#Phylogenetic regression with MPLE model
#Eusociality as a function of haplodiploidy
hap_glm_1<-phyloglm(Eusociality~Haplodiploidy,data =data_hap1,phy= tree_hap1,method = c("logistic_MPLE"),boot=1000)
summary(hap_glm_1)
#Haplodiploidy as a function of eusocuality
hap_glm_2<-phyloglm(Haplodiploidy~Eusociality,data =data_hap1,phy= tree_hap1,method = c("logistic_MPLE"),boot=1000)
summary(hap_glm_2)

#For first analysis: Removing eusociality presence for aphids
tree_hap1_pr<-read.nexus("Datafile S1.nexus")
data_hap1_pr<-read.csv("Datafile S3_aphid.csv",header=TRUE,row.names = 1)


#Phylogenetic regression with MPLE model
#Eusociality as a function of haplodiploidy
hap_glm_1a<-phyloglm(Eusociality~Haplodiploidy,data =data_hap1_pr,phy= tree_hap1_pr,method = c("logistic_MPLE"),boot=1000)
summary(hap_glm_1a)
#Haplodiploidy as a function of eusocuality
hap_glm_1b<-phyloglm(Haplodiploidy~Eusociality,data =data_hap1_pr,phy= tree_hap1_pr,method = c("logistic_MPLE"),boot=1000)
summary(hap_glm_1b)


#For second analysis: Removing eusociality presence for aphids and termites
tree_hap2_pr<-read.nexus("Datafile S1.nexus")
data_hap2_pr<-read.csv("Datafile S3_aphid_termite.csv",header=TRUE,row.names = 1)

#Phylogenetic regression with MPLE model
#Eusociality as a function of haplodiploidy
hap_glm_2a<-phyloglm(Eusociality~Haplodiploidy,data =data_hap2_pr,phy= tree_hap2_pr,method = c("logistic_MPLE"),boot=1000)
summary(hap_glm_2a)
#Haplodiploidy as a function of eusocuality
hap_glm_2b<-phyloglm(Haplodiploidy~Eusociality,data =data_hap2_pr,phy= tree_hap2_pr,method = c("logistic_MPLE"),boot=1000)
summary(hap_glm_2b)


#For third analysis: Removing haplodiploidy presence for some of the families
tree_hap3_pr<-read.nexus("Datafile S1.nexus")
data_hap3_pr<-read.csv("Datafile S3_hap_removal.csv",header=TRUE,row.names = 1)

#Phylogenetic regression with MPLE model
#Eusociality as a function of haplodiploidy
hap_glm_3a<-phyloglm(Eusociality~Haplodiploidy,data =data_hap3_pr,phy= tree_hap3_pr,method = c("logistic_MPLE"),boot=1000)
summary(hap_glm_3a)
#Haplodiploidy as a function of eusocuality
hap_glm_3b<-phyloglm(Haplodiploidy~Eusociality,data =data_hap3_pr,phy= tree_hap3_pr,method = c("logistic_MPLE"),boot=1000)
summary(hap_glm_3b)

#FOR THE DATASET WHERE THE SPECIES 

