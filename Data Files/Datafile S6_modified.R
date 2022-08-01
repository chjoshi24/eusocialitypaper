library(phytools)
library(treedata.table)

###### Haplodiploidy-eusociality analysis ######

#Loading the data and the tree. 

tree_hap1<-read.nexus("Datafile S1-family-level tree for hexapoda.nexus")
data_hap1<-read.csv("Datafile S3.csv",header=TRUE,row.names = 1)
data_hap1$t1<-row.names(data_hap1)
td1<-as.treedata.table(tree_hap1,data_hap1)

#Checking which model is better for each of the traits
x1<-setNames(data_hap1$Haplodiploidy,rownames(data_hap1))
y1<-setNames(data_hap1$Eusociality,rownames(data_hap1))
#Eusociality
eus_ER1<-fitMk(tree_hap1,y,model="ER")
eus_ER1
eus_ARD1<-fitMk(tree_hap1,y,model="ARD")
eus_ARD1
aic_eus1<-setNames(sapply(list(eus_ER1,eus_ARD1),AIC),c("ER","ARD"))
aic.w(aic_eus1)

#Haplodiploidy 
hap_ER1<-fitMk(tree_hap1,x,model="ER")
hap_ER1
hap_ARD1<-fitMk(tree_hap1,x,model="ARD")
hap_ARD1
aic_hap1<-setNames(sapply(list(hap_ER1,hap_ARD1),AIC),c("ER","ARD"))
aic.w(aic_hap1)

#Running the pagel's test as an ARD model(by default)

#Eusociality as a function of haplodiploidy
model_hap_ARD1<-tdt(td1, phytools::fitPagel(phy, x=extractVector(td1,"Haplodiploidy"),y=extractVector(td1, "Eusociality"), 
                                       dep.var="y"))
model_hap_ARD1

#Haploidiploidy as a function of eusociality
model_hap_ARD2<-tdt(td1, phytools::fitPagel(phy, x=extractVector(td1,"Haplodiploidy"),y=extractVector(td1, "Eusociality"), 
                                       dep.var="x"))
model_hap_ARD2

#Both traits depending on each other
model_hap_ARD3<-tdt(td1, phytools::fitPagel(phy, x=extractVector(td1,"Haplodiploidy"),y=extractVector(td1, "Eusociality"), 
                                       dep.var="xy"))
model_hap_ARD3

#AIC weight calculation to find which model has the highest AIC weight
aic_hap_eus1<-c(model_hap_ARD1$dependent.AIC,model_hap_ARD2$dependent.AIC,model_hap_ARD3$dependent.AIC,model_hap_ARD3$independent.AIC)
aic.w(aic_hap_eus1)

#Running the pagel's test as an ER model 

#Eusociality as a function of haplodiploidy
model_hap_ER1<-tdt(td1, phytools::fitPagel(phy, x=extractVector(td1,"Haplodiploidy"),y=extractVector(td1, "Eusociality"), 
                                      dep.var="y", model="ER"))
model_hap_ER1

#Haploidploidy as a function of eusociality
model_hap_ER2<-tdt(td1, phytools::fitPagel(phy, x=extractVector(td1,"Haplodiploidy"),y=extractVector(td1, "Eusociality"), 
                                      dep.var="x", model="ER"))
model_hap_ER2

#Both traits depending on each other
model_hap_ER3<-tdt(td1, phytools::fitPagel(phy, x=extractVector(td1,"Haplodiploidy"),y=extractVector(td1, "Eusociality"), 
                                      dep.var="xy",model="ER"))
model_hap_ER3

#AIC weight calculation to find which model has the highest AIC weight
aic_hap_eus2<-c(model_hap_ER1$dependent.AIC,model_hap_ER2$dependent.AIC,model_hap_ER3$dependent.AIC,model_hap_ER3$independent.AIC)
aic.w(aic_hap_eus2)

#Sequential bonferroni correction to correct for p-values across a table.
#ARD model
p_1<-c(0.000150128 ,0.001148938 ,0.001394496)
adjusted_1<-p.adjust(p_1,method="holm",length(p_1)) 
#ER model
p_2<-c(2.756891*10^-6,0.1925039 ,1.57633*10^-5)
adjusted_2<-p.adjust(p_2,method="holm",length(p_2)) 

##########################################################################

###### Haplodiploidy-eusociality analysis ######

#Here the seven families in Phthiraptera and four families in Protura were
#removed from the analysis due to lack of sex determination system data.

#Loading the data and removing the missing taxa for the analysis
data_hap2<-read.csv("Datafile S4.csv",header=TRUE,row.names = 1)
data_hap2$t2<-row.names(data_hap2)
data_hap2<-data_hap2[!is.na(data_hap2$Haplodiploidy),]
td2<-as.treedata.table(tree_hap1,data_hap2)
tree_hap2<-td2$phy

#Checking which model is better for each of the traits
x2<-setNames(data_hap2$Haplodiploidy,rownames(data_hap2))
y2<-setNames(data_hap2$Eusociality,rownames(data_hap2))

#Eusociality
eus_ER2<-fitMk(tree_hap2,y2,model="ER")
eus_ARD2<-fitMk(tree_hap2,y2,model="ARD")
aic_eus2<-setNames(sapply(list(eus_ER2,eus_ARD2),AIC),c("ER","ARD"))
aic.w(aic_eus2)

#Haplodiploidy
hap_ER2<-fitMk(tree_hap2,x2,model="ER")
hap_ARD2<-fitMk(tree_hap2,x2,model="ARD")
aic_hap2<-setNames(sapply(list(hap_ER2,hap_ARD2),AIC),c("ER","ARD"))
aic.w(aic_hap2)

#Running the pagel's test as an ARD model(by default)

#Eusociality as a function of haplodiploidy
model_hap_ARD4<-tdt(td2, phytools::fitPagel(phy, x=extractVector(td2,"Haplodiploidy"),y=extractVector(td2, "Eusociality"), 
                                       dep.var="y"))
model_hap_ARD4

#Haplodiploidy as a function of eusociality
model_hap_ARD5<-tdt(td2, phytools::fitPagel(phy, x=extractVector(td2,"Haplodiploidy"),y=extractVector(td2, "Eusociality"), 
                                            dep.var="x"))
model_hap_ARD5

#Both traits depending on each other
model_hap_ARD6<-tdt(td2, phytools::fitPagel(phy, x=extractVector(td2,"Haplodiploidy"),y=extractVector(td2, "Eusociality"), 
                                            dep.var="xy"))
model_hap_ARD6


#AIC weight calculation to find which model has the highest AIC weight
aic_hap_eus3<-c(model_hap_ARD4$dependent.AIC,model_hap_ARD5$dependent.AIC,model_hap_ARD6$dependent.AIC,model_hap_ARD6$independent.AIC)
aic.w(aic_hap_eus3)

#Running the pagel's test as an ER model

#Eusociality as a function of haplodiploidy
model_hap_ER4<-tdt(td2, phytools::fitPagel(phy, x=extractVector(td2,"Haplodiploidy"),y=extractVector(td2, "Eusociality"), 
                                            dep.var="y",model="ER"))
model_hap_ER4

#Haplodiploidy as a function of eusociality
model_hap_ER5<-tdt(td2, phytools::fitPagel(phy, x=extractVector(td2,"Haplodiploidy"),y=extractVector(td2, "Eusociality"), 
                                           dep.var="x",model="ER"))
model_hap_ER5

#Both traits depending on each other
model_hap_ER6<-tdt(td2, phytools::fitPagel(phy, x=extractVector(td2,"Haplodiploidy"),y=extractVector(td2, "Eusociality"), 
                                           dep.var="xy",model="ER"))
model_hap_ER6


#AIC weight calculation to find which model has the highest AIC weight
aic_hap_eus4<-c(model_hap_ER4$dependent.AIC,model_hap_ER5$dependent.AIC,model_hap_ER6$dependent.AIC,model_hap_ER6$independent.AIC)
aic.w(aic_hap_eus4)

#Sequential bonferroni correction to correct for p-values across a table.

#ARD model
p_3<-c(0.000792046,0.000160087,0.00149033)
adjusted_3<-p.adjust(p_3,method="holm",length(p_3)) 

#ER model
p_4<-c(3.26388*10^-6,0.196287 ,1.85176*10^-5)
adjusted_4<-p.adjust(p_4,method="holm",length(p_4)) 

###########################################################################

###### Haplodiploidy-eusociality analysis ######

#Here the haplodiploidy presence was removed for the Curculonidae family 
#because the eusocial beetle is diplodiploid.


#Loading the data and using the same tree used for the main haplodiploidy analysis
data_hap3<-read.csv("Datafile S5.csv",header=TRUE,row.names = 1)
data_hap3$t3<-row.names(data_hap3)
td3<-as.treedata.table(tree_hap1,data_hap3)

#Checking which model is better for each of the traits

x3<-setNames(data_hap3$Haplodiploidy,rownames(data_hap3))
y3<-setNames(data_hap3$Eusociality,rownames(data_hap3))

#Eusociality
eus_ER3<-fitMk(tree_hap1,y3,model="ER")
eus_ARD3<-fitMk(tree_hap1,y3,model="ARD")
aic_eus3<-setNames(sapply(list(eus_ER3,eus_ARD3),AIC),c("ER","ARD"))
aic.w(aic_eus3)

#Haplodiploidy
hap_ER3<-fitMk(tree_hap1,x3,model="ER")
hap_ARD3<-fitMk(tree_hap1,x3,model="ARD")
aic_hap3<-setNames(sapply(list(hap_ER3,hap_ARD3),AIC),c("ER","ARD"))
aic.w(aic_hap3)

#Running the pagel's test as an ARD model(by default)

#Eusociality as a function of haplodiploidy

model_hap_ARD7<-tdt(td3, phytools::fitPagel(phy, x=extractVector(td3, "Haplodiploidy"),
                                       y=extractVector(td3, "Eusociality"), 
                                       dep.var="y"))
model_hap_ARD7


#Haploidploidy as a function of eusociality
model_hap_ARD8<-tdt(td3, phytools::fitPagel(phy, x=extractVector(td3, "Haplodiploidy"),
                                            y=extractVector(td3, "Eusociality"), 
                                            dep.var="x"))
model_hap_ARD8


#Both traits depending on each other
model_hap_ARD9<-tdt(td3, phytools::fitPagel(phy, x=extractVector(td3, "Haplodiploidy"),
                                            y=extractVector(td3, "Eusociality"), 
                                            dep.var="xy"))
model_hap_ARD9

#AIC weight calculation to find which model has the highest AIC weight

aic_hap_eus5<-c(model_hap_ARD7$dependent.AIC,model_hap_ARD8$dependent.AIC,model_hap_ARD9$dependent.AIC,model_hap_ARD9$independent.AIC)
aic.w(aic_hap_eus5)

#Running the pagel's test as an ER model 

#Eusociality as a function of haplodiploidy

model_hap_ER7<-tdt(td3, phytools::fitPagel(phy, x=extractVector(td3, "Haplodiploidy"),
                                            y=extractVector(td3, "Eusociality"), 
                                            dep.var="y",model="ER"))
model_hap_ER7


#Haploidploidy as a function of eusociality
model_hap_ER8<-tdt(td3, phytools::fitPagel(phy, x=extractVector(td3, "Haplodiploidy"),
                                            y=extractVector(td3, "Eusociality"), 
                                            dep.var="x",model="ER"))
model_hap_ER8


#Both traits depending on each other
model_hap_ER9<-tdt(td3, phytools::fitPagel(phy, x=extractVector(td3, "Haplodiploidy"),
                                            y=extractVector(td3, "Eusociality"), 
                                            dep.var="xy",model="ER"))
model_hap_ER9

#AIC weight calculation to find which model has the highest AIC weight

aic_hap_eus6<-c(model_hap_ER7$dependent.AIC,model_hap_ER8$dependent.AIC,model_hap_ER9$dependent.AIC,model_hap_ER9$independent.AIC)
aic.w(aic_hap_eus6)

#Sequential bonferroni correction to correct for p-values across a table.
#ARD model
p_5<-c(0.00289994 ,1,0.0189305)
adjusted_p_5<-p.adjust(p_5,method="holm",length(p_5))

#ER model
p_6<-c(3.74279*10^-5,0.668376,0.000185882)
adjusted_p_6<-p.adjust(p_6,method="holm",length(p_6))

########################################################################

###### Monogamy-eusociality analysis ######

#Here the tree has families and superfamilies which are constrained
#to be monophyletic.

#Loading the tree and data files
data_mon1<-read.csv("Datafile S11.csv",header=TRUE,row.names = 1)
tree_mon1<-read.tree("Datafile S8-species-level tree for hymenoptera-1.tre")
data_mon1$t4<-row.names(data_mon1)
td4<-as.treedata.table(tree_mon1,data_mon1)

#Checking which model is better for each of the traits
tree_mon2<-td4$phy
x4<-setNames(data_mon1$Eusociality,rownames(data_mon1))
y4<-setNames(data_mon1$Monogamy,rownames(data_mon1))

#Eusociality
eus_ER4<-fitMk(tree_mon2,x4,model="ER")
eus_ARD4<-fitMk(tree_mon2,x4,model="ARD")
aic_eus4<-setNames(sapply(list(eus_ER4,eus_ARD4),AIC),c("ER","ARD"))
aic.w(aic_eus4)

#Monogamy
mon_ER1<-fitMk(tree_mon2,y4,model="ER")
mon_ARD1<-fitMk(tree_mon2,y4,model="ARD")
aic_mon1<-setNames(sapply(list(mon_ER1,mon_ARD1),AIC),c("ER","ARD"))
aic.w(aic_mon1)

#Running the models as ARD (by default)

#Eusociality as a function of monogamy
model_mon_ARD1<-tdt(td4, phytools::fitPagel(phy, x=extractVector(td4, "Monogamy"),
                                           y=extractVector(td4, "Eusociality"), 
                                           dep.var="y"))
model_mon_ARD1

#Monogamy as a function of eusociality
model_mon_ARD2<-tdt(td4, phytools::fitPagel(phy, x=extractVector(td4, "Monogamy"),
                                            y=extractVector(td4, "Eusociality"), 
                                            dep.var="x"))
model_mon_ARD2

#Both traits depending on each other
model_mon_ARD3<-tdt(td4, phytools::fitPagel(phy, x=extractVector(td4, "Monogamy"),
                                            y=extractVector(td4, "Eusociality"), 
                                            dep.var="xy"))
model_mon_ARD3

#AIC weight calculation to find which model has the highest AIC weight
aic_mon_eus1<-c(model_mon_ARD1$dependent.AIC,model_mon_ARD2$dependent.AIC,model_mon_ARD3$dependent.AIC,model_mon_ARD3$independent.AIC)
aic.w(aic_mon_eus1)


#Running the models as ER

#Eusociality as a function of monogamy
model_mon_ER1<-tdt(td4, phytools::fitPagel(phy, x=extractVector(td4, "Monogamy"),
                                            y=extractVector(td4, "Eusociality"), 
                                            dep.var="y",model="ER"))
model_mon_ER1

#Monogamy as a function of eusociality
model_mon_ER2<-tdt(td4, phytools::fitPagel(phy, x=extractVector(td4, "Monogamy"),
                                            y=extractVector(td4, "Eusociality"), 
                                            dep.var="x",model="ER"))
model_mon_ER2

#Both traits depending on each other
model_mon_ER3<-tdt(td4, phytools::fitPagel(phy, x=extractVector(td4, "Monogamy"),
                                            y=extractVector(td4, "Eusociality"), 
                                            dep.var="xy",model="ER"))
model_mon_ER3

#AIC weight calculation to find which model has the highest AIC weight
aic_mon_eus2<-c(model_mon_ER1$dependent.AIC,model_mon_ER2$dependent.AIC,model_mon_ER3$dependent.AIC,model_mon_ER3$independent.AIC)
aic.w(aic_mon_eus2)

#Sequential bonferroni correction to correct for p-values across a table.

#ARD model
p_7<-c(0.894726,7.32278*10^-6,3.42805*10^-5)
adjusted_7<-p.adjust(p_7,method="holm",length(p_7)) 

#ER model
p_8<-c(0.502152,1.90461*10^-6,3.94085*10^-6)
adjusted_8<-p.adjust(p_8,method="holm",length(p_8)) 


#########################################################################

###### Monogamy-eusociality analysis ######

#Running similar analysis for another tree where only the superfamilies were 
#constrained to be monophyletic.

#Loading the tree and data files
#Using the same datafile as used for previous tree

tree_mon3<-read.tree("Datafile S9- species-level tree for hymenoptera-2.tre")
td5<-as.treedata.table(tree_mon3,data_mon1)

#Checking which model is better for each of the traits
tree_mon4<-td5$phy
x5<-setNames(data_mon1$Eusociality,rownames(data_mon1))
y5<-setNames(data_mon1$Monogamy,rownames(data_mon1))

 
#Eusociality
eus_ER5<-fitMk(tree_mon4,x5,model="ER")
eus_ARD5<-fitMk(tree_mon4,x5,model="ARD")
aic_eus5<-setNames(sapply(list(eus_ER5,eus_ARD5),AIC),c("ER","ARD"))
aic.w(aic_eus5)

#Monogamy
mon_ER2<-fitMk(tree_mon4,y5,model="ER")
mon_ARD2<-fitMk(tree_mon4,y5,model="ARD")
aic_mon2<-setNames(sapply(list(mon_ER2,mon_ARD2),AIC),c("ER","ARD"))
aic.w(aic_mon2)

#Running the models as ARD (by default)

#Eusociality as a function of monogamy

model_mon_ARD4<-tdt(td5, phytools::fitPagel(phy, x=extractVector(td5, "Monogamy"),
                                            y=extractVector(td5, "Eusociality"), 
                                            dep.var="y"))
model_mon_ARD4

#Monogamy as a function of eusociality

model_mon_ARD5<-tdt(td5, phytools::fitPagel(phy, x=extractVector(td5, "Monogamy"),
                                            y=extractVector(td5, "Eusociality"), 
                                            dep.var="x"))
model_mon_ARD5

#Both traits depending on each other

model_mon_ARD6<-tdt(td5, phytools::fitPagel(phy, x=extractVector(td5, "Monogamy"),
                                            y=extractVector(td5, "Eusociality"), 
                                            dep.var="xy"))
model_mon_ARD6

#AIC weight calculation to find which model has the highest AIC weight
aic_mon_eus3<-c(model_mon_ARD4$dependent.AIC,model_mon_ARD5$dependent.AIC,model_mon_ARD6$dependent.AIC,model_mon_ARD6$independent.AIC)
aic.w(aic_mon_eus3)

#Sequential bonferroni correction to correct for p-values across a table.

#ARD model
p_9<-c(0.74399,3.32243*10^-8,2.67928*10^-8)
adjusted_9<-p.adjust(p_9,method="holm",length(p_9)) 

#########################################################################

###### Monogamy-eusociality analysis ######

#Rerunning the monogamy-eusociality analysis for both trees.
#Some of the species belonging to same genus were represented by single label
# Two of genera (Apanteles and Nomia) had variation in monogamy presence.
# So the monogamy presence is reassigned in the datasheets and the analysis is run again.

#Here the tree has families and superfamilies which are constrained
#to be monophyletic.

#Loading the tree and data files
data_mon2<-read.csv("Datafile S12.csv",header=TRUE,row.names = 1)
tree_mon5<-read.tree("Datafile S8-species-level tree for hymenoptera-1.tre")
data_mon2$t5<-row.names(data_mon2)
td6<-as.treedata.table(tree_mon5,data_mon2)

#Checking which model is better for each of the traits
tree_mon6<-td6$phy
x6<-setNames(data_mon2$Eusociality,rownames(data_mon2))
y6<-setNames(data_mon2$Monogamy,rownames(data_mon2))


#Eusociality
eus_ER6<-fitMk(tree_mon6,x6,model="ER")
eus_ARD6<-fitMk(tree_mon6,x6,model="ARD")
aic_eus6<-setNames(sapply(list(eus_ER6,eus_ARD6),AIC),c("ER","ARD"))
aic.w(aic_eus6)

#Monogamy
mon_ER3<-fitMk(tree_mon6,y6,model="ER")
mon_ARD3<-fitMk(tree_mon6,y6,model="ARD")
aic_mon3<-setNames(sapply(list(mon_ER3,mon_ARD3),AIC),c("ER","ARD"))
aic.w(aic_mon3)

#Running the models as ARD (by default)

#Eusociality as a function of monogamy

model_mon_ARD7<-tdt(td6, phytools::fitPagel(phy, x=extractVector(td6, "Monogamy"),
                                            y=extractVector(td6, "Eusociality"), 
                                            dep.var="y"))
model_mon_ARD7

#Monogamy as a function of eusociality

model_mon_ARD8<-tdt(td6, phytools::fitPagel(phy, x=extractVector(td6, "Monogamy"),
                                            y=extractVector(td6, "Eusociality"), 
                                            dep.var="x"))
model_mon_ARD8


#Both traits depending on each other

model_mon_ARD9<-tdt(td6, phytools::fitPagel(phy, x=extractVector(td6, "Monogamy"),
                                            y=extractVector(td6, "Eusociality"), 
                                            dep.var="xy"))
model_mon_ARD9


#AIC weight calculation to find which model has the highest AIC weight
aic_mon_eus4<-c(model_mon_ARD7$dependent.AIC,model_mon_ARD8$dependent.AIC,model_mon_ARD9$dependent.AIC,model_mon_ARD9$independent.AIC)
aic.w(aic_mon_eus4)

#Running the models as ER

model_mon_ER4<-tdt(td6, phytools::fitPagel(phy, x=extractVector(td6, "Monogamy"),
                                            y=extractVector(td6, "Eusociality"), 
                                            dep.var="y",model="ER"))
model_mon_ER4

#Monogamy as a function of eusociality

model_mon_ER5<-tdt(td6, phytools::fitPagel(phy, x=extractVector(td6, "Monogamy"),
                                            y=extractVector(td6, "Eusociality"), 
                                            dep.var="x",model="ER"))
model_mon_ER5


#Both traits depending on each other

model_mon_ER6<-tdt(td6, phytools::fitPagel(phy, x=extractVector(td6, "Monogamy"),
                                            y=extractVector(td6, "Eusociality"), 
                                            dep.var="xy",model="ER"))
model_mon_ER6


#AIC weight calculation to find which model has the highest AIC weight
aic_mon_eus4<-c(model_mon_ER4$dependent.AIC,model_mon_ER5$dependent.AIC,model_mon_ER6$dependent.AIC,model_mon_ER6$independent.AIC)
aic.w(aic_mon_eus4)

#Sequential bonferroni correction to correct for p-values across a table.

#ARD model
p_10<-c(0.838275,1.26259*10^-5,5.06072*10^-5)
adjusted_10<-p.adjust(p_10,method="holm",length(p_10)) 

#ER model
p_11<-c(0.508471,3.56815*10^-6,6.37051*10^-6)
adjusted_11<-p.adjust(p_11,method="holm",length(p_11)) 

######################################################################

###### Monogamy-eusociality analysis ######


#Running similar analysis for another tree where only the superfamilies were 
#constrained to be monophyletic.

#Loading the tree and data files
#Using the same datafile as used for previous tree
tree_mon7<-read.tree("Datafile S9- species-level tree for hymenoptera-2.tre")
td7<-as.treedata.table(tree_mon7,data_mon2)

#Checking which model is better for each of the traits
tree_mon8<-td7$phy
x7<-setNames(data_mon2$Eusociality,rownames(data_mon2))
y7<-setNames(data_mon2$Monogamy,rownames(data_mon2))


#Eusociality
eus_ER7<-fitMk(tree_mon8,x7,model="ER")
eus_ARD7<-fitMk(tree_mon8,x7,model="ARD")
aic_eus7<-setNames(sapply(list(eus_ER7,eus_ARD7),AIC),c("ER","ARD"))
aic.w(aic_eus7)

#Monogamy
mon_ER3<-fitMk(tree_mon8,y7,model="ER")
mon_ARD3<-fitMk(tree_mon8,y7,model="ARD")
aic_mon3<-setNames(sapply(list(mon_ER3,mon_ARD3),AIC),c("ER","ARD"))
aic.w(aic_mon3)

#Running the models as ARD (by default)

#Eusociality as a function of monogamy

model_mon_ARD10<-tdt(td7, phytools::fitPagel(phy, x=extractVector(td7, "Monogamy"),
                                            y=extractVector(td7, "Eusociality"), 
                                            dep.var="y"))
model_mon_ARD10

#Monogamy as a function of eusociality

model_mon_ARD11<-tdt(td7, phytools::fitPagel(phy, x=extractVector(td7, "Monogamy"),
                                             y=extractVector(td7, "Eusociality"), 
                                             dep.var="x"))
model_mon_ARD11

#Both traits depending on each other
model_mon_ARD12<-tdt(td7, phytools::fitPagel(phy, x=extractVector(td7, "Monogamy"),
                                             y=extractVector(td7, "Eusociality"), 
                                             dep.var="xy"))
model_mon_ARD12


#AIC weight calculation to find which model has the highest AIC weight
aic_mon_eus5<-c(model_mon_ARD10$dependent.AIC,model_mon_ARD11$dependent.AIC,model_mon_ARD12$dependent.AIC,model_mon_ARD12$independent.AIC)
aic.w(aic_mon_eus5)

#Sequential bonferroni correction to correct for p-values across a table.

#ARD model
p_12<-c(1,6.00335*10^-8,3.50135*10^-8)
adjusted_12<-p.adjust(p_12,method="holm",length(p_12)) 

########################################################################
#Phylogenetic regression 
tree_hap1_pr<-read.nexus("Datafile S1.nexus")
data_hap1_pr<-read.csv("Datafile S3.csv",header=TRUE,row.names = 1)
library(phylolm)
hap_glm_1<-phyloglm(Eusociality~Haplodiploidy,data =data_hap1_pr,phy=tree_hap1_pr,method = c("logistic_IG10"))
    
