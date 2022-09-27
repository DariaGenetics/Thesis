library(corrplot)
library("data.table")
library(dplyr)
library("ggpubr")
library(ggfortify)
library(tidyr)
library("writexl")
library("nortest")
library(mycor)
library("writexl")
DGRPData <-fread("/Users/dgg/Desktop/CompliedDariaP12_013122.txt") # reads the doc, expands library
DGRPDataNoNAs #DGRPData with NAs replaced with respective phenotypes mean 
DGRPCor #Correlation coefficients for all variables from DGRPDataNoNAs

## Creating a corrplot for 4 variables with NAs removed  -----------------------------------------

library(tidyverse)
DGRPDataTest <- (DGRPData[, c(1:5)])
names(DGRPDataTest)[names(DGRPDataTest) == "1Week-AcclimationSurvival_NegativeSixDegreesC_F"] <- "FA"
names(DGRPDataTest)[names(DGRPDataTest) == "1Week-AcclimationSurvival_NegativeSixDegreesC_M"] <- "MA"
names(DGRPDataTest)[names(DGRPDataTest) == "1Week-Non-acclimationSurvival_NegativeSixDegrees_F"] <- "FNA"
names(DGRPDataTest)[names(DGRPDataTest) == "1Week-Non-acclimationSurvival_NegativeSixDegrees_M"] <- "MNA"
na.omit(DGRPDataTest)

VarFour <- var(DGRPDataTest[,2:5], na.rm = TRUE)
Correlation <- cor(DGRPDataTest [,2:5], use="pairwise.complete.obs")
corrplot(Correlation, method = "ellipse")

## Creating a corrplot for 4 variables without subdividing table -----------------------------------------

VarFour <- var(DGRPData[,2:5], na.rm = TRUE)
names(DGRPData)[names(DGRPData) == "1Week-AcclimationSurvival_NegativeSixDegreesC_F"] <- "FA"
names(DGRPData)[names(DGRPData) == "1Week-AcclimationSurvival_NegativeSixDegreesC_M"] <- "MA"
names(DGRPData)[names(DGRPData) == "1Week-Non-acclimationSurvival_NegativeSixDegrees_F"] <- "FNA"
names(DGRPData)[names(DGRPData) == "1Week-Non-acclimationSurvival_NegativeSixDegrees_M"] <- "MNA"
Correlation <- cor(DGRPData [,2:5], use="pairwise.complete.obs")
corrplot(Correlation, method = "ellipse")

## Creating a corrplot for 4 variables without subdividing table and with NA's Present -----------------------------------------

names(DGRPData)[names(DGRPData) == "1Week-AcclimationSurvival_NegativeSixDegreesC_F"] <- "OneWeekFA"
names(DGRPData)[names(DGRPData) == "1Week-AcclimationSurvival_NegativeSixDegreesC_M"] <- "OneWeekMA"
names(DGRPData)[names(DGRPData) == "1Week-Non-acclimationSurvival_NegativeSixDegrees_F"] <- "OneWeekFNA"
names(DGRPData)[names(DGRPData) == "1Week-Non-acclimationSurvival_NegativeSixDegrees_M"] <- "OneWeekMNA"
names(DGRPData)[names(DGRPData) == "1Week-StarvationResistance_StarvationMedium_F"] <- "StarvationF"
names(DGRPData)[names(DGRPData) == "1Week-StarvationResistance_StarvationMedium_M"] <- "StarvationM"
    colnames(DGRPData)
Correlation <- cor(DGRPData [,2:7], use="pairwise.complete.obs")
corrplot(Correlation, method = "ellipse")

Correlation <- cor(DGRPData [,c(2,3,4,5)], use="pairwise.complete.obs")
corrplot(Correlation, method = "ellipse")

## Covariance matrixes -----------------------------------------------------

DGRPData <-fread("/Users/dgg/Desktop/CompliedDariaP12_013122.txt") # reads the doc, expands library
names(DGRPData)[names(DGRPData) == "1Week-AcclimationSurvival_NegativeSixDegreesC_F"] <- "OneWeekFA"
names(DGRPData)[names(DGRPData) == "1Week-AcclimationSurvival_NegativeSixDegreesC_M"] <- "OneWeekMA"
names(DGRPData)[names(DGRPData) == "1Week-Non-acclimationSurvival_NegativeSixDegrees_F"] <- "OneWeekFNA"
names(DGRPData)[names(DGRPData) == "1Week-Non-acclimationSurvival_NegativeSixDegrees_M"] <- "OneWeekMNA"
names(DGRPData)[names(DGRPData) == "1Week-StarvationResistance_StarvationMedium_F"] <- "StarvationF"
names(DGRPData)[names(DGRPData) == "1Week-StarvationResistance_StarvationMedium_M"] <- "StarvationM"

CovarianceTestPartialDataset <- cov(DGRPData[2:5], use = 'pairwise.complete.obs')
    #caution using this strategy: http://bwlewis.github.io/covar/missing.html
    #note: the above calculated the covariance test of columns 2:5 to *all other columns*
CovarianceTestWholeDataset <- cov(DGRPData, use = 'pairwise.complete.obs')

## Isolating desired covariance matrixes -----------------------------------------------------

  #finding maximum and minimum correlation values 
   DGRPData <-fread("/Users/dgg/Desktop/CompliedDariaP12_013122.txt") # reads the doc, expands library
   DGRPData$ral_id <- NULL
   CovarianceTestWholeDataset <- cov(DGRPData, use = 'pairwise.complete.obs')
   summary(CovarianceTestWholeDataset)
   
   DTCovarianceTestWholeDataset <- data.frame(CovarianceTestWholeDataset)
   DTCovarianceTestWholeDataset[DTCovarianceTestWholeDataset >= "0.5"]
   summary(DTCovarianceTestWholeDataset)

   DTCovarianceTestWholeDataset[rowSums(DTCovarianceTestWholeDataset[1:103] > 1500) > 0, ]
     #Egg production 
     # not sure what the '1500' value is actually selecting for 
   
   
   Correlation <- cor(DGRPData [,c(2,3,4,5)], use="pairwise.complete.obs")
   corrplot(Correlation, method = "ellipse")
   
   CorrelationEGG <- cor(DGRPData [,c(31, 72, 82, 96)], use="pairwise.complete.obs", )
   corrplot(CorrelationEGG, method = "ellipse")

colnames(DGRPData)
   
   EggProduction_Standard_F #31
   Mortality_Oneμg.vialForMalathionTenμg.vialForPermethrin_F #72
   ResistanceToDDT_15HourExposureTo200μlAcetone.DDTSolutionConcentration05μg.ml_F #82
   StarvationResistance_NoFood_Larvae.mixed #94
   SurvivalLength_Ma549Infection_F #96
   
   names(DGRPData)[names(DGRPData) == "EggProduction_Standard_F"] <- "A"
   names(DGRPData)[names(DGRPData) == "Mortality_Oneμg-vialForMalathionTenμg-vialForPermethrin_F"] <- "B"
   names(DGRPData)[names(DGRPData) == "ResistanceToDDT_15HourExposureTo200μlAcetone-DDTSolutionConcentration05μg-ml_F"] <- "C"
   names(DGRPData)[names(DGRPData) == "StarvationResistance_NoFood_Larvae-mixed"] <- "D"
   names(DGRPData)[names(DGRPData) == "SurvivalLength_Ma549Infection_F"] <- "E"

## Attempting correlation coefficients -----------------------------------------------------
   
   #Fa and Ma
   ggscatter(DGRPData, x = "FA", y = "MA", 
             add = "reg.line", conf.int = TRUE, 
             cor.coef = TRUE, cor.method = "pearson",
             xlab = "Survival female", ylab = "survival male")
   cor.test(DGRPData$FA, DGRPData$FA, 
            method = "pearson", use = "complete.obs")
   
   #EggProduction_Standard_F and Mortality_Oneμg-vialForMalathionTenμg-vialForPermethrin_F
   
   ggscatter(DGRPData, x = "EggProduction_Standard_F", y = "Mortality_Oneμg-vialForMalathionTenμg-vialForPermethrin_F", 
             add = "reg.line", conf.int = TRUE, 
             cor.coef = TRUE, cor.method = "pearson", use = "complete.obs",
             xlab = "X", ylab = "Y")
   cor.test(DGRPData$A, DGRPData$B, 
            method = "pearson", use = "complete.obs")
       #note - swapped to 'A' and 'B' since B was having trouble in the extended form
   
  #correlation for all variables (1: remove NAs in one phenotype, 2: remove NAs throughout datatable, 3: run round(cor(...))) again)
       round(cor(DGRPData, use = "everything"))
          #'no complete pairs' 
          #'pairwise.complete.pairs' returns a lot of zeros, also sketchy way of dealing with NAs 
          #'everything' returns NAs for basically everything 

      #replacing NAs in SurvivalLength_Standard_F (easier phenotype to see when looking at the readout)
         DGRPData %>% 
            mutate(SurvivalLength_Standard_F = replace_na(SurvivalLength_Standard_F,mean(SurvivalLength_Standard_F, na.rm = TRUE)))
         
      #replacing NAs through DGRPData 
         DGRPDataNoNAs <- DGRPData %>% 
            mutate_if(is.numeric, ~replace_na(.,mean(., na.rm = TRUE)))
      
      #finding correlation for the whole table
         #with centered means as NAs
           DGRPCor <- cor(DGRPDataNoNAs)
         #with NAs omitted   ------------######
           DGRPCarOmittedNas <- cor.test(DGRPData[,2], DGRPData[,3], na.action=na.omit)
           DGRPCarOmittedNas <- mycor(DGRPData, method = "pearson", na.action=na.omit)
        
         #counting values above 0.7, below -0.7 (considering >0.7, <-0.7 to be significant) 
         sum(DGRPCor > 0.7 & DGRPCor < 1) #168 (168 results of positive correlation)
         sum(DGRPCor < -0.7 & DGRPCor > -1) #2 (two numbers of strong negative correlation)
           min(DGRPCor) #-0.8244451
            
         #Print as xlsx file
            write_xlsx(DGRPCorFrame, 'Test.xlsx')

##-------------covariance UNFINISHED
            
         DGRPCov <- cov(DGRPDataNoNAs)   
         DGRPCovPrint <- data.frame(DGRPCov)
         write_xlsx(DGRPCovPrint, 'DGRP Covariance.xlsx')

# finding high coveriances 
colnames(DGRPData)
 #groups to chekc:
   #1) col 59-62 ([59]"Lifespan_TwentyEightDegreesC_F"                                                
              "Lifespan_TwentyEightDegreesC_M"                                                
              "Lifespan_TwentyFiveDegreesC_F"                                                 
              "Lifespan_TwentyFiveDegreesC_M"  
       #col 1-11 (B-N)
       # col 63-64 "HighThermalToleranceExtreme_VaryingWithTemperature_F"                          
               "HighThermalToleranceExtreme_VaryingWithTemperature_M"           
   
       #col 101-102 "SurvivalLength_Standard_F"                                                     
                    "SurvivalLength_Standard_M"

               
## Correlation dendrogram to see what phenotypes are related 
 sim.by.hclust <- hclust(dist(DGRPCor))
plot(sim.by.hclust, cex=0.50, main = "Phenotype Correlation Dendrogram")


##Covariance dendrogram ##cannot nullify "EggProduction_Standard_F" or "ral_id" DOES NOT WORK
DGRPcovprep <- DGRPDataNoNAs
DGRPcovprep$"EggProduction_Standard_F" <- NULL
DGRPcovprep$"ral_id" <- NULL
DGRPcov <- cov(DGRPcovprep)
as.data.frame(DGRPcov)
sim.by.hclust <- hclust(dist(DGRPCov))
plot(sim.by.hclust, cex=0.50, main = "Phenotype Covariance Dendrogram") 

##Normalicy Check ~~ Shapiro Test, p-value > 0.05 is normal
DGRPNormalicyST <- apply(DGRPDataNoNAs,2,shapiro.test)
capture.output(summary(DGRPNormalicy), file = "DGRPNormalicy.txt") ##doesn't save the results, only that the test was conducted
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
  
#PCA
  library(ggfortify)
  DGRPpca <- princomp(DGRPCor, cor = FALSE, scores = TRUE )
  autoplot(DGRPpca)

##loading the full dataset:
  FullDGRPData <-readRDS("/Users/dgg/Desktop/wideform.phenotypedata.RDS") # reads the doc, expands library

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
#ST
   DGRPDataSTRAIGHTLINESNoNAsPCA <- princomp(DGRPDataSTRAIGHTLINESNoNAsCor, cor = FALSE, scores = TRUE )
autoplot(DGRPDataSTRAIGHTLINESNoNAsPCA, loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3) #<-- generates PCA plot
  #minimizing principle components
    DGRPDataSTRAIGHTLINESNoNAsPrComp <- prcomp(DGRPDataSTRAIGHTLINESNoNAsCor, center = TRUE, scale = FALSE)
    plot(cumsum(DGRPDataSTRAIGHTLINESNoNAsPrComp$sdev^2/sum(DGRPDataSTRAIGHTLINESNoNAsPrComp$sdev^2)))
        #takeaways: the first 6 principle components explain ~70% of the variation 
    pc.ST.use <- 6
    trunc <- DGRPDataSTRAIGHTLINESNoNAsPrComp[,1:pc.ST.use] #### Ongoing problems 

   
#INV
DGRPDataCURVYLINESNoNAsCorPCA <- princomp(DGRPDataCURVYLINESNoNAsCor, cor = FALSE, scores = TRUE )
autoplot(DGRPDataCURVYLINESNoNAsCorPCA, loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3) #<-- generates PCA plot
  #minimizing principle components
DGRPDataCURVYLINESNoNAsCorPCA <- as.numeric(unlist(DGRPDataCURVYLINESNoNAsCorPCA))
    DGRPDataCURVYLINESNoNAsPrComp <- prcomp(DGRPDataCURVYLINESNoNAsCorPCA, center = TRUE, scale = FALSE)
    plot(cumsum(DGRPDataSTRAIGHTLINESNoNAsPrComp$sdev^2/sum(DGRPDataSTRAIGHTLINESNoNAsPrComp$sdev^2)))
        #takeaways: the first 6 principle compoents explain ~70% of the variation 
  

~~~~~~~~~~~~cPCA with INV/ST lines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ALL BELOW IS CODE COPIED NONADAPTED
       
         scpca_df <- scpca_sim$x %>%
   as_tibble() %>%
   mutate(label = toy_df[, 31] %>% as.character)
colnames(scpca_df) <- c("scPC1", "scPC2", "label")

# plot the results
p_scpca <- ggplot(scpca_df, aes(x = scPC1, y = scPC2, colour = label)) +
   geom_point(alpha = 0.5) +
   ggtitle("scPCA of Simulated Data") +
   theme_minimal()
p_scpca
   
   
   
   
   

write_xlsx(DGRPDataSTRAIGHTLINESNoNAsCor,"/Users/dgg/Desktop/DGRPDataSTRAIGHTLINESNoNAsCor.xlsx")
  
  
  
  
  
  
  
  
  
  


                    