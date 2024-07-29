# Header -----------------------------------------------------------------------

# AUTHOR: John Eric Hamilton
# PURPOSE: Turn data files into a more informative format to include in publication
# NOTES: 
#   *Inputs:   - Harmonised Data Files  
#              - MR Result Files
#              - GTEx Summary Statistics Files
#              - Steiger Directionality Files
#
#   *Outputs:  - CSV files summaries for data 

# Load Libraries ---------------------------------------------------------------
library(openxlsx)

# Load in relevant files -------------------------------------------------------

# Raw MR Results ---------------------------------------------------------------
#allMR <- readRDS("Final/Pancreas/Data_MR1/MR_Result/allResult_allNonSignifAndSignif.rds")
#write.xlsx(allMR, 'Final/Pancreas/Data_MR1/allMR.xlsx')

# MR Filtered Results ----------------------------------------------------------
# Pancreas MR filtered results 
allFinalMR_Pancreas <- data.frame()
for (i in c(1, 2, 4)){
  finalMR <- readRDS(sprintf("Final/Pancreas/Data_MR%s/MR_Result/allResults_filteredWithNames.rds", i))
  if (i == 1){
    finalMR$MR <- "Chiou"
  } else if (i == 2){
    finalMR$MR <- "Michalek Eur"
  } else if (i == 4){
    finalMR$MR <- "Michalek Transethnic"
  }
  
  allFinalMR_Pancreas <- rbind(allFinalMR_Pancreas, finalMR)
}
write.xlsx(allFinalMR_Pancreas, 'Final/allMR_filtered_Pancreas.xlsx')

# WB MR filtered results
allFinalMR_WB <- data.frame()
for (i in c(5, 6, 7)){
  finalMR <- readRDS(sprintf("Final/WB/Data_MR%s/MR_Result/allResults_filteredWithNames.rds", i))
  if (i == 5){
    finalMR$MR <- "Chiou"
  } else if (i == 6){
    finalMR$MR <- "Michalek Eur"
  } else if (i == 7){
    finalMR$MR <- "Michalek Transethnic"
  }
  
  allFinalMR_WB <- rbind(allFinalMR_WB, finalMR)
}

# Spleen MR filtered results
write.xlsx(allFinalMR_WB, 'Final/allMR_filtered_WB.xlsx')

allFinalMR_Spleen <- data.frame()
for (i in c(8, 9, 10)){
  finalMR <- readRDS(sprintf("Final/Spleen/Data_MR%s/MR_Result/allResults_filteredWithNames.rds", i))
  if (i == 8){
    finalMR$MR <- "Chiou"
  } else if (i == 9){
    finalMR$MR <- "Michalek Eur"
  } else if (i == 10){
    finalMR$MR <- "Michalek Transethnic"
  }
  
  allFinalMR_Spleen <- rbind(allFinalMR_Spleen, finalMR)
}

write.xlsx(allFinalMR_Spleen, 'Final/allMR_filtered_Spleen.xlsx')


# Harmonised Data Results ------------------------------------------------------
# Harmonised Pancreas Data (for only SNPs whose gene_id are in the finalMR data)
allHarm_Pancreas <- data.frame()
allHarm_Pancreas_MR2 <- data.frame()

for (i in c(1)){
  for (j in 1:22){
    harm <- readRDS(sprintf("Final/Pancreas/Data_MR%s/Harmonised_Data/Har_Chr%s.rds", i, j))
    harm$MR <- "Chiou"
    allHarm_Pancreas <- rbind(allHarm_Pancreas, harm)
  }
}

for (i in c(2, 4)){
  for (j in 1:22){
    harm <- readRDS(sprintf("Final/Pancreas/Data_MR%s/Harmonised_Data/Har_Chr%s.rds", i, j))
    if (i == 2){
      harm$MR <- "Michalek Eur"
    } else if (i == 4){
      harm$MR <- "Michalek Transethnic"
    }
    allHarm_Pancreas_MR2 <- rbind(allHarm_Pancreas_MR2, harm)
  }
}

colnames(allHarm_Pancreas_MR2)[c(14)] <- c("variant_id")
allHarm_Pancreas_MR2 <- select(allHarm_Pancreas_MR2, -chromosome, -base_pair_location)

allHarm_Pancreas <- rbind(allHarm_Pancreas, allHarm_Pancreas_MR2)
allHarm_Pancreas_filtered <- allHarm_Pancreas[allHarm_Pancreas$id.exposure%in%allFinalMR_Pancreas$id.exposure,]

write.xlsx(allHarm_Pancreas_filtered, 'Final/allHarm_Pancreas.xlsx')



# Harmonised WB Data (for only SNPs whose gene_id are in the finalMR data)
allHarm_WB <- data.frame()
allHarm_WB_MR2 <- data.frame()

for (i in c(5)){
  for (j in 1:22){
    harm <- readRDS(sprintf("Final/WB/Data_MR%s/Harmonised_Data/Har_Chr%s.rds", i, j))
    harm$MR <- "Chiou"
    allHarm_WB <- rbind(allHarm_WB, harm)
  }
}

for (i in c(6, 7)){
  for (j in 1:22){
    harm <- readRDS(sprintf("Final/WB/Data_MR%s/Harmonised_Data/Har_Chr%s.rds", i, j))
    if (i == 6){
      harm$MR <- "Michalek Eur"
    } else if (i == 7){
      harm$MR <- "Michalek Transethnic"
    }
    allHarm_WB_MR2 <- rbind(allHarm_WB_MR2, harm)
  }
}

colnames(allHarm_WB_MR2)[c(14)] <- c("variant_id")
allHarm_WB_MR2 <- select(allHarm_WB_MR2, -chromosome, -base_pair_location)

allHarm_WB <- rbind(allHarm_WB, allHarm_WB_MR2)
allHarm_WB_filtered <- allHarm_WB[allHarm_WB$id.exposure%in%allFinalMR_WB$id.exposure,]

write.xlsx(allHarm_WB_filtered, 'Final/allHarm_WB.xlsx')


# Harmonised Spleen Data (for only SNPs whose gene_id are in the finalMR data)
allHarm_Spleen <- data.frame()
allHarm_Spleen_MR2 <- data.frame()

for (i in c(8)){
  for (j in 1:22){
    harm <- readRDS(sprintf("Final/Spleen/Data_MR%s/Harmonised_Data/Har_Chr%s.rds", i, j))
    harm$MR <- "Chiou"
    allHarm_Spleen <- rbind(allHarm_Spleen, harm)
  }
}

for (i in c(9, 10)){
  for (j in 1:22){
    harm <- readRDS(sprintf("Final/Spleen/Data_MR%s/Harmonised_Data/Har_Chr%s.rds", i, j))
    if (i == 9){
      harm$MR <- "Michalek Eur"
    } else if (i == 10){
      harm$MR <- "Michalek Transethnic"
    }
    allHarm_Spleen_MR2 <- rbind(allHarm_Spleen_MR2, harm)
  }
}

colnames(allHarm_Spleen_MR2)[c(14)] <- c("variant_id")
allHarm_Spleen_MR2 <- select(allHarm_Spleen_MR2, -chromosome, -base_pair_location)

allHarm_Spleen <- rbind(allHarm_Spleen, allHarm_Spleen_MR2)
allHarm_Spleen_filtered <- allHarm_Spleen[allHarm_Spleen$id.exposure%in%allFinalMR_Spleen$id.exposure,]

write.xlsx(allHarm_Spleen_filtered, 'Final/allHarm_Spleen.xlsx')

# GTEx Summary Statistics Data -------------------------------------------------
# Pancreas GTEx summary statistics
GTEX_Pancreas <- df[df$gene_id%in%allFinalMR_Pancreas$id.exposure,]
write.xlsx(GTEX_Pancreas, 'Final/GTEX_Pancreas.xlsx')

# WB GTEx summary statistics
GTEX_WB <- df_WB[df_WB$gene_id%in%allFinalMR_WB$id.exposure,]
write.xlsx(GTEX_WB, 'Final/GTEX_WB.xlsx')

# Spleen GTEx summary statistics
GTEX_Spleen <- df_Spleen[df_Spleen$gene_id%in%allFinalMR_Spleen$id.exposure,]
write.xlsx(GTEX_Spleen, 'Final/GTEX_Spleen.xlsx')

# Steiger ----------------------------------------------------------------------

# Pancreas
steiger_pancreas <- data.frame()
steiger_pancreas_MR2 <- data.frame()

for (i in c(1)){
  steiger <- readRDS(sprintf('Final/Steiger/steiger_MR%s.rds', i))
  steiger_pancreas <- rbind(steiger_pancreas, steiger)
  steiger_pancreas$MR <- "Chiou"
}

for (i in c(2, 4)){
  steiger <- readRDS(sprintf('Final/Steiger/steiger_MR%s.rds', i))
  if (i == 2){
    steiger$MR <- "Michalek Eur"
  } else if (i == 4){
    steiger$MR <- "Michalek Transethnic"
  }
  steiger_pancreas_MR2 <- rbind(steiger_pancreas_MR2, steiger)
}

colnames(steiger_pancreas_MR2)[c(14)] <- c("variant_id")
steiger_pancreas_MR2 <- select(steiger_pancreas_MR2, -chromosome, -base_pair_location)

steiger_pancreas <- rbind(steiger_pancreas, steiger_pancreas_MR2)

write.xlsx(steiger_pancreas, 'Final/Steiger_Pancreas.xlsx')


# WB
steiger_WB <- data.frame()
steiger_WB_MR2 <- data.frame()

for (i in c(5)){
  steiger <- readRDS(sprintf('Final/Steiger/steiger_MR%s.rds', i))
  steiger_WB <- rbind(steiger_WB, steiger)
  steiger_WB$MR <- "Chiou"
}

for (i in c(6, 7)){
  steiger <- readRDS(sprintf('Final/Steiger/steiger_MR%s.rds', i))
  if (i == 6){
    steiger$MR <- "Michalek Eur"
  } else if (i == 7){
    steiger$MR <- "Michalek Transethnic"
  }
  steiger_WB_MR2 <- rbind(steiger_WB_MR2, steiger)
}

colnames(steiger_WB_MR2)[c(14)] <- c("variant_id")
steiger_WB_MR2 <- select(steiger_WB_MR2, -chromosome, -base_pair_location)

steiger_WB <- rbind(steiger_WB, steiger_WB_MR2)

write.xlsx(steiger_WB, 'Final/Steiger_WB.xlsx')



# Spleen
steiger_Spleen <- data.frame()
steiger_Spleen_MR2 <- data.frame()

for (i in c(8)){
  steiger <- readRDS(sprintf('Final/Steiger/steiger_MR%s.rds', i))
  steiger_Spleen <- rbind(steiger_Spleen, steiger)
  steiger_Spleen$MR <- "Chiou"
}

for (i in c(9, 10)){
  steiger <- readRDS(sprintf('Final/Steiger/steiger_MR%s.rds', i))
  if (i == 9){
    steiger$MR <- "Michalek Eur"
  } else if (i == 10){
    steiger$MR <- "Michalek Transethnic"
  }
  steiger_Spleen_MR2 <- rbind(steiger_Spleen_MR2, steiger)
}

colnames(steiger_Spleen_MR2)[c(14)] <- c("variant_id")
steiger_Spleen_MR2 <- select(steiger_Spleen_MR2, -chromosome, -base_pair_location)

steiger_Spleen <- rbind(steiger_Spleen, steiger_Spleen_MR2)

write.xlsx(steiger_Spleen, 'Final/Steiger_Spleen.xlsx')
