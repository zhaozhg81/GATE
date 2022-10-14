
## #################
## #################
## #################
## #################
## #################
## #################
## #################
## #################
## #################
## This part of the code extract the data from AYP 2013
AYP_temp <- read.csv("Independent_AYP_2013_SD.csv", sep=",",header=T)
AYP <- AYP_temp[AYP_temp[,3]=="E",c(3,7,8,122,123,124,194,195,196)]

for( i in c(4:9) )
  AYP[,i] <- as.numeric( AYP[,i] )

AYP.SAN <- AYP[,4]-AYP[,7]
AYP.SAP <- 100*(AYP[,5]-AYP[,8])/AYP.SAN

AYP1 <- cbind(AYP[,-c(5,8)],AYP.SAN,AYP.SAP)

AYP2 <- AYP1[AYP1[,6]>=2 & AYP1[,8] >=2 & !is.na(AYP1[,9]),]
AYP2 <- AYP1[AYP1[,6]>=20 & AYP1[,8] >=20 & !is.na(AYP1[,9]) & ((AYP1[,7] > 0 & AYP1[,9] >=0) | AYP1[,7]>=0 & AYP1[,9]>0),]
length(AYP2[,1])
DIFF <- median(AYP2[,9]) - median(AYP2[,7])
Z <- ((AYP2[,9] - AYP2[,7] - DIFF)/100)/sqrt((AYP2[,7]/100*(1-AYP2[,7]/100))/AYP2[,6]+(AYP2[,9]/100*(1-AYP2[,9]/100))/AYP2[,8])

analysis <- cbind(AYP2,Z)

analysis <- analysis[abs(analysis$Z) < 10,]
analysis <- analysis[order(analysis$DNAME),]

AYP13 <- analysis[, c(2,3,10) ]

## #################
## #################
## #################
## #################
## #################
## #################
## #################
## #################
## #################
## This part of the code extract the data from AYP 2015
AYP_temp <- read.csv("apr15a.txt", sep="\t",header=T)

## 6 School name
## 7 School district
## 17 number of students who tested math
## 113 number of students who are proficient or above in math
## 114 Percentage of students who are proficient or above in math
## 89 number of students with socioeconomically disadvantaged who tested in math
## 185 number of students with SED who are proficient in math
## 186 percentage of students with SED who are proficient in math
AYP1 <- AYP_temp[ , c(3,6,7,17,113,114, 89,185, 186)]
for( i in c(4:9) )
  AYP1[,i] <- as.numeric( as.character(AYP1[,i] ) )


AYP.SAN <- AYP1[,4]-AYP1[,7]
AYP.SAP <- 100*(AYP1[,5]-AYP1[,8])/AYP.SAN

AYP1 <- cbind(AYP1[,-c(5,8)],AYP.SAN,AYP.SAP)

AYP2 <- AYP1[AYP1[,6]>=2 & AYP1[,8] >=2 & !is.na(AYP1[,9]),]
AYP2 <- AYP1[AYP1[,6]>=20 & AYP1[,8] >=20 & !is.na(AYP1[,9]) & ((AYP1[,7] > 0 & AYP1[,9] >=0) | AYP1[,7]>=0 & AYP1[,9]>0),]
length(AYP2[,1])
DIFF <- median(AYP2[,9]) - median(AYP2[,7])
Z <- ((AYP2[,9] - AYP2[,7] - DIFF)/100)/sqrt((AYP2[,7]/100*(1-AYP2[,7]/100))/AYP2[,6]+(AYP2[,9]/100*(1-AYP2[,9]/100))/AYP2[,8])
AYP2 <- cbind( AYP2, Z )

## Keep only the elementary school
keep.ind <- which( AYP2[,1]=="E" )
AYP2 <- AYP2[keep.ind,]
AYP2 <- AYP2[order(AYP2$dname),]
AYP15 <- AYP2[, c(2,3,10)]

write.csv("AYP13.csv", AYP13)

write.table(AYP13, "AYP13_processed.csv", sep=",")
write.table(AYP15, "AYP15_processed.csv", sep=",")
