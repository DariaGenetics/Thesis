#Extremely Condensed Script

library("data.table")
library(corrplot)
library("ggpubr")
library(dplyr)
library(tidyr)
library(ggfortify)
DGRPData <-fread("/Users/dgg/Desktop/CompliedDariaP12_013122.txt") # reads the doc, expands library

#DGRPDataNoNAs -- DGRPData with NAs replaced with respective phenotypes mean 
DGRPDataNoNAs <- DGRPData %>% 
  mutate_if(is.numeric, ~replace_na(.,mean(., na.rm = TRUE)))

#DGRPCor -- Correlation coefficients for all variables from DGRPDataNoNAs
DGRPCor <- cor(DGRPDataNoNAs)
DGRPCorFrame <- data.frame(DGRPCor) #DGRPCorFrame -- transfers DGRPCor from atomic vectors into a data frame 

##Normalicy Check ~~ Shapiro Test, p-value<0.05 is normal
DGRPNormalicyST <- apply(DGRPDataNoNAs,2,shapiro.test)
DGRPNormalicySTTable <- sapply(DGRPNormalicyST, unlist)
DGRPNormalicySTTableDF <- data.table(DGRPNormalicyTable)
write_xlsx(DGRPNormalicySTTableDF, "C:\\Users\\dgg\\Desktop\\Thesis\\DGRPNormalicy.xlsx")
   #with data containing NAs
   DGRPNormalicySTWNAs <- apply(DGRPData,2,shapiro.test)
   DGRPNormalicySTWNAsTable <- sapply(DGRPNormalicySTWNAs, unlist)
   DGRPNormalicySTWNAsTableDF <- data.table(DGRPNormalicySTWNAsTable)
   write_xlsx(DGRPNormalicySTWNAsTableDF, "C:\\Users\\dgg\\Desktop\\Thesis\\DGRPNormalicySTWNAsTableDF.xlsx")

##Normalicy Check ~~ Anderson-Darling test, p-value > 0.05 is normal
DGRPNormalicyAD <- sapply(DGRPDataNoNAs, ad.test)
write.csv(DGRPNormalicyAD, "C:\\Users\\dgg\\Desktop\\Thesis\\DGRPNormalityADNoNAs.csv")
   #with data containing NAs
   DGRPNormalicyADWNAs <- sapply(DGRPData, ad.test)
   write.csv(DGRPNormalicyADWNAs, "C:\\Users\\dgg\\Desktop\\Thesis\\DGRPNormalityADWNAs.csv")
   
## Correlation dendrogram to see what phenotypes are related 
sim.by.hclustDGRPCor <- hclust(dist(DGRPCor))
plot(sim.by.hclustDGRPCor, cex=0.40, main = "Phenotype All Lines Correlation Dendrogram")

#PCA
DGRPpca <- princomp(DGRPCor, cor = FALSE, scores = TRUE )
autoplot(DGRPpca, loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)


~~~~~~~~~~TTrying to get full Dataset~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FullDGRPData <-readRDS("/Users/dgg/Desktop/wideform.phenotypedata.RDS") # reads the doc, expands library
#FDGRPDataNoNAs -- DGRPData with NAs replaced with respective phenotypes mean 
FullDGRPDataNoNAs <- FullDGRPData %>% 
   mutate_if(is.numeric, ~replace_na(.,mean(., na.rm = TRUE)))
data_numeric <- as.numeric(FullDGRPDataNoNAs)


#DGRPCor -- Correlation coefficients for all variables from DGRPDataNoNAs ----- ## can't figure out how to coerce FullDGRPDataNoNas into a numeric 
FullDGRPCor <- cor(FullDGRPFrame)
FullDGRPCorFrame <- data.frame(FullDGRPCor) #DGRPCorFrame -- transfers DGRPCor from atomic vectors into a data frame 

~~~~~~~~~~~Specific~Lines INV or not~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ##Finding specific lines with the inversion
   DGRPData <-fread("/Users/dgg/Desktop/CompliedDariaP12_013122.txt") # reads the doc, expands library
InversionData <- fread("/Users/dgg/Desktop/inversion2.csv")
InversionDatatable <- as.data.table(InversionData)
#ST = staight 
#INV = inversion
#Plan: mesh the inversion table to the DGRP data table, then split by what they possess in the ln(2L)t table 
DGRPDataWithInversions <- merge(DGRPData,InversionDatatable,by="ral_id")
write_xlsx(DGRPDataWithInversions,"/Users/dgg/Desktop/DGRPDataWithInversions.xlsx")

##reupload separated dataset 
DGRPDataSTRAIGHTLINES <-fread("/Users/dgg/Desktop/DGRPDataSTRAIGHT.csv")
DGRPDataCURVYLINES <-fread("/Users/dgg/Desktop/DGRPDataCURVY.csv")


#DGRPDataNoNAs -- DGRPData with NAs replaced with respective phenotypes mean 
DGRPDataSTRAIGHTLINESNoNAs <- DGRPDataSTRAIGHTLINES %>% 
   mutate_if(is.numeric, ~replace_na(.,mean(., na.rm = TRUE)))
DGRPDataSTRAIGHTLINESNoNAs = subset (DGRPDataSTRAIGHTLINESNoNAs, select = - (V105)) #<-- deleting this extra column takes out NAs, allows calculations
   is.na(DGRPDataSTRAIGHTLINESNoNAs) %>% table() #Resolved deleting col(V105) trying to find why hclust() for cor doesn't work - 158 contained NAS

DGRPDataCURVYLINESNoNAs <- DGRPDataCURVYLINES %>% 
   mutate_if(is.numeric, ~replace_na(.,mean(., na.rm = TRUE)))
DGRPDataCURVYLINESNoNAs = subset (DGRPDataCURVYLINESNoNAs, select = - (V105))
   is.na(DGRPDataCURVYLINESNoNAs) %>% table()

#DGRPCor -- Correlation coefficients for all variables from DGRPDataNoNAs
DGRPDataSTRAIGHTLINESNoNAsCor <- cor(DGRPDataSTRAIGHTLINESNoNAs) #<- note that 'ral_id' treated as a phenotype here 
DGRPDataCURVYLINESNoNAsCor <- cor(DGRPDataCURVYLINESNoNAs) # <- note that 'ral_id' treated as a phenotype here


## Correlation dendrogram to see what phenotypes are related 
sim.by.hclust.DGRPDataSTRAIGHTLINESNoNAsCor <- hclust(dist(DGRPDataSTRAIGHTLINESNoNAsCor))
plot(sim.by.hclust.DGRPDataSTRAIGHTLINESNoNAsCor, cex=0.40, main = "Phenotype ST Lines Correlation Dendrogram")#<- non-square dataset - snippping a line?

sim.by.hclust.DGRPDataCURVYLINESNoNAs <- hclust(dist(DGRPDataCURVYLINESNoNAs))
plot(sim.by.hclust.DGRPDataCURVYLINESNoNAs, cex=0.40, main = "Phenotype INV Lines Correlation Dendrogram") #<- non-square datset - snipping the data?

~~~~~~~~~~~~PCA with INV/ST lines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#PCA
DGRPDataSTRAIGHTLINESNoNAsPCA <- princomp(DGRPDataSTRAIGHTLINESNoNAsCor, cor = FALSE, scores = TRUE )
autoplot(DGRPDataSTRAIGHTLINESNoNAsPCA, loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)

DGRPDataCURVYLINESNoNAsCorPCA <- princomp(DGRPDataCURVYLINESNoNAsCor, cor = FALSE, scores = TRUE )
autoplot(DGRPDataCURVYLINESNoNAsCorPCA, loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)
