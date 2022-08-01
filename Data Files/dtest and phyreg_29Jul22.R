 
library(phytools)
library(treedata.table)
library(phylolm)

#FOR THE ORIGINAL DATASET 

#No manipulations in the data
tree_phyreg1<-read.nexus("Datafile S1.nexus")
data_phyreg1<-read.csv("Datafile S3.csv",header=TRUE,row.names = 1)

#Phylogenetic regression with MPLE model
#Eusociality as a function of haplodiploidy
phyreg_glm1a<-phyloglm(Eusociality~Haplodiploidy,data =data_phyreg1,phy=tree_phyreg1,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm1a)
#Haplodiploidy as a function of eusocuality
phyreg_glm1b<-phyloglm(Haplodiploidy~Eusociality,data =data_phyreg1,phy= tree_phyreg1,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm1b)

#Removing eusociality presence for aphids
tree_phyreg2<-read.nexus("Datafile S1.nexus")
data_phyreg2<-read.csv("Datafile S3_aphid.csv",header=TRUE,row.names = 1)


#Phylogenetic regression with MPLE model
#Eusociality as a function of haplodiploidy
phyreg_glm2a<-phyloglm(Eusociality~Haplodiploidy,data =data_phyreg2,phy= tree_phyreg2,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm2a)
#Haplodiploidy as a function of eusocuality
phyreg_glm2b<-phyloglm(Haplodiploidy~Eusociality,data =data_phyreg2,phy= tree_phyreg2,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm2b)


#Removing eusociality presence for aphids and termites
tree_phyreg3<-read.nexus("Datafile S1.nexus")
data_phyreg3<-read.csv("Datafile S3_aphid_termite.csv",header=TRUE,row.names = 1)

#Phylogenetic regression with MPLE model
#Eusociality as a function of haplodiploidy
phyreg_glm3a<-phyloglm(Eusociality~Haplodiploidy,data =data_phyreg3,phy= tree_phyreg3,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm3a)
#Haplodiploidy as a function of eusocuality
phyreg_glm3b<-phyloglm(Haplodiploidy~Eusociality,data =data_phyreg3,phy= tree_phyreg3,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm3b)


#Removing haplodiploidy presence for some of the families
tree_phyreg4<-read.nexus("Datafile S1.nexus")
data_phyreg4<-read.csv("Datafile S3_hap_removal.csv",header=TRUE,row.names = 1)

#Phylogenetic regression with MPLE model
#Eusociality as a function of haplodiploidy
phyreg_glm4a<-phyloglm(Eusociality~Haplodiploidy,data =data_phyreg4,phy= tree_phyreg4,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm4a)
#Haplodiploidy as a function of eusocuality
phyreg_glm4b<-phyloglm(Haplodiploidy~Eusociality,data =data_phyreg4,phy= tree_phyreg4,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm4b)


#FOR THE DATASET THAT HAS 11 FAMILIES REMOVED 

#No manipulations in the data
tree_phyreg5<-read.tree("Datafile S7.tre")
data_phyreg5<-read.csv("Datafile S4.csv",header=TRUE,row.names = 1)

 
#Phylogenetic regression with MPLE model
#Eusociality as a function of haplodiploidy
phyreg_glm5a<-phyloglm(Eusociality~Haplodiploidy,data =data_phyreg5,phy=tree_phyreg5,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm5a)
#Haplodiploidy as a function of eusocuality
phyreg_glm5b<-phyloglm(Haplodiploidy~Eusociality,data =data_phyreg5,phy=tree_phyreg5,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm5b)


#Removing eusociality presence for aphids
tree_phyreg6<-read.tree("Datafile S7.tre")
data_phyreg6<-read.csv("Datafile S4_aphid.csv",header=TRUE,row.names = 1)

#Phylogenetic regression with MPLE model
#Eusociality as a function of haplodiploidy
phyreg_glm6a<-phyloglm(Eusociality~Haplodiploidy,data =data_phyreg6,phy=tree_phyreg6,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm6a)
#Haplodiploidy as a function of eusocuality
phyreg_glm6b<-phyloglm(Haplodiploidy~Eusociality,data =data_phyreg6,phy=tree_phyreg6,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm6b)

#Removing eusociality presence for aphids and termites
tree_phyreg7<-read.tree("Datafile S7.tre")
data_phyreg7<-read.csv("Datafile S4_aphid_termite.csv",header=TRUE,row.names = 1)

#Phylogenetic regression with MPLE model
#Eusociality as a function of haplodiploidy
phyreg_glm7a<-phyloglm(Eusociality~Haplodiploidy,data =data_phyreg7,phy=tree_phyreg7,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm7a)
#Haplodiploidy as a function of eusocuality
phyreg_glm7b<-phyloglm(Haplodiploidy~Eusociality,data =data_phyreg7,phy=tree_phyreg7,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm7b)


#Removing haplodiploidy presence for some of the families

tree_phyreg8<-read.tree("Datafile S7.tre")
data_phyreg8<-read.csv("Datafile S4_hap_removal.csv",header=TRUE,row.names = 1)


#Phylogenetic regression with MPLE model
#Eusociality as a function of haplodiploidy
phyreg_glm8a<-phyloglm(Eusociality~Haplodiploidy,data =data_phyreg8,phy= tree_phyreg8,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm8a)
#Haplodiploidy as a function of eusocuality
phyreg_glm8b<-phyloglm(Haplodiploidy~Eusociality,data =data_phyreg8,phy= tree_phyreg8,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm8b)



#FOR THE DATASET WHERE THE BEETLE FAMILY THAT HAS EUSOCIALITY IS CONSIDERED AS DIPLODIPLOID


#No manipulations in the data
tree_phyreg9<-read.nexus("Datafile S1.nexus")
data_phyreg9<-read.csv("Datafile S5.csv",header=TRUE,row.names = 1)

#Phylogenetic regression with MPLE model
#Eusociality as a function of haplodiploidy
phyreg_glm9a<-phyloglm(Eusociality~Haplodiploidy,data =data_phyreg9,phy=tree_phyreg9,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm9a)
#Haplodiploidy as a function of eusocuality
phyreg_glm9b<-phyloglm(Haplodiploidy~Eusociality,data =data_phyreg9,phy=tree_phyreg9,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm9b)

#Removing eusociality presence for aphids
tree_phyreg10<-read.nexus("Datafile S1.nexus")
data_phyreg10<-read.csv("Datafile S5_aphid.csv",header=TRUE,row.names = 1)

#Phylogenetic regression with MPLE model
#Eusociality as a function of haplodiploidy
phyreg_glm10a<-phyloglm(Eusociality~Haplodiploidy,data =data_phyreg10,phy=tree_phyreg10,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm10a)
#Haplodiploidy as a function of eusocuality
phyreg_glm10b<-phyloglm(Haplodiploidy~Eusociality,data =data_phyreg10,phy=tree_phyreg10,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm10b)

#Removing eusociality presence for aphids and termites
tree_phyreg11<-read.nexus("Datafile S1.nexus")
data_phyreg11<-read.csv("Datafile S5_aphid_termite.csv",header=TRUE,row.names = 1)

#Phylogenetic regression with MPLE model
#Eusociality as a function of haplodiploidy
phyreg_glm11a<-phyloglm(Eusociality~Haplodiploidy,data =data_phyreg11,phy=tree_phyreg11,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm11a)
#Haplodiploidy as a function of eusocuality
phyreg_glm11b<-phyloglm(Haplodiploidy~Eusociality,data =data_phyreg11,phy=tree_phyreg11,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm11b)

#Removing haplodiploidy presence for some of the families

tree_phyreg12<-read.nexus("Datafile S1.nexus")
data_phyreg12<-read.csv("Datafile S5_hap_removal.csv",header=TRUE,row.names = 1)

#Phylogenetic regression with MPLE model
#Eusociality as a function of haplodiploidy
phyreg_glm12a<-phyloglm(Eusociality~Haplodiploidy,data =data_phyreg12,phy=tree_phyreg12,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm12a)
#Haplodiploidy as a function of eusocuality
phyreg_glm12b<-phyloglm(Haplodiploidy~Eusociality,data =data_phyreg12,phy=tree_phyreg12,method = c("logistic_MPLE"),boot=1000)
summary(phyreg_glm12b)

#D-TEST 

#FOR THE ORIGINAL DATASET 

tree_hap1_pr<-read.nexus("Datafile S1.nexus")
data_hap1_pr<-read.csv("Datafile S3.csv",header=TRUE,row.names = 1)


#Extracting the hapldoiploidy and eusociality with tips
data_hap1_pr_hap<-setNames(data_hap1_pr[,1],rownames(data_hap1_pr))
data_hap1_pr_eus<-setNames(data_hap1_pr[,2],rownames(data_hap1_pr))

#Making simmaps first

#simmaps for haplodiploidy
map_hap_ER1<-make.simmap(tree_hap1_pr,data_hap1_pr_hap,model="ER",nsim=100)
map_hap_ARD1<-make.simmap(tree_hap1_pr,data_hap1_pr_hap,model="ARD",nsim=100)

#simmaps for eusociality
map_eus_ER1<-make.simmap(tree_hap1_pr,data_hap1_pr_eus,model="ER",nsim=100)
map_eus_ARD1<-make.simmap(tree_hap1_pr,data_hap1_pr_eus,model="ARD",nsim=100)

#Running D test
dtest_ER1<-Dtest(map_hap_ER1,map_eus_ER1,nsim=100)
dtest_ARD1<-Dtest(map_hap_ARD1,map_eus_ARD1,nsim=100)

dtest_ER_ARD_1<-Dtest(map_hap_ARD1,map_eus_ER1,nsim=100)
dtest_ER_ARD_2<-Dtest(map_hap_ER1,map_eus_ARD1,nsim=100)


#FOR THE DATASET THAT HAS 11 FAMILIES REMOVED 

tree_hap2_pr<-read.tree("Datafile S7.tre")
data_hap2_pr<-read.csv("Datafile S4.csv",header=TRUE,row.names = 1)

#Extracting the hapldoiploidy and eusociality with tips
data_hap2_pr_hap<-setNames(data_hap2_pr[,1],rownames(data_hap2_pr))
data_hap2_pr_eus<-setNames(data_hap2_pr[,2],rownames(data_hap2_pr))

#Making simmaps first

#simmaps for haplodiploidy
map_hap_ER2<-make.simmap(tree_hap2_pr,data_hap2_pr_hap,model="ER",nsim=100)
map_hap_ARD2<-make.simmap(tree_hap2_pr,data_hap2_pr_hap,model="ARD",nsim=100)

#simmaps for eusociality
map_eus_ARD2<-make.simmap(tree_hap2_pr,data_hap2_pr_eus,model="ARD",nsim=100)

#Running D test
dtest_ARD2<-Dtest(map_hap_ARD2,map_eus_ARD2,nsim=100)
dtest_ER_ARD_3<-Dtest(map_hap_ER2,map_eus_ARD2,nsim=100)

#FOR THE DATASET WHERE THE BEETLE FAMILY THAT HAS EUSOCIALITY IS CONSIDERED AS DIPLODIPLOID

d<-2
