# Header ----------------------------------------------------------------------

# AUTHOR: John Eric Hamilton
# PURPOSE: Plotting results from the Data collection (03)
# NOTES: 
#   *Inputs:   - MR Results filtered and non-filtered 
#
#   *Outputs:  - Forest Plot of significant genes based on OR 
#              - Gene Name Conversion Table (from Gene Ensemble ID to common gene name)

# Library and Data Loading ----------------------------------------------------
library(tidyverse)

# Load the filtered (for significant p-value) MR results
result_signif_total <- readRDS('New/Pancreas/Data_MR1/MR_Result/allResultSignif.rds')
filtered <- readRDS('New/Pancreas/Data_MR1/MR_Result/allResults_filtered.rds')
filtered_MHC <- readRDS('New/Pancreas/Data_MR1/MR_Result/allResults_filtered_MHC_region.rds')

# Making a Conversion Table --------------------------------------------------
gtf <- read.table('Homo_sapiens.GRCh38.112.gene.gtf', sep='\t')
gtf_reduced <- data.frame()

for (i in 1:dim(gtf)[1]){
  gtf_temp <- strsplit(gtf[[9]][i], split = '; ')
  
  # Length 4 means that there is no gene name in the gtf
  if (length(gtf_temp[[1]]) == 4){
    id <- strsplit(gtf_temp[[1]], split = ' ')[[1]][2]
    biotype <- strsplit(gtf_temp[[1]], split = ' ')[[4]][2]
    name <- paste(id, biotype, sep='_')
    
    gtf_temp_row <- data.frame(id, name, biotype)
    gtf_reduced <- rbind(gtf_reduced, gtf_temp_row)
    
  }
  
  # Length 5 means that there is no gene name in the gtf
  if (length(gtf_temp[[1]]) == 5){
    id <- strsplit(gtf_temp[[1]], split = ' ')[[1]][2]
    name <- strsplit(gtf_temp[[1]], split = ' ')[[3]][2]
    biotype <- strsplit(gtf_temp[[1]], split = ' ')[[5]][2]
    
    gtf_temp_row <- data.frame(id, name, biotype)
    gtf_reduced <- rbind(gtf_reduced, gtf_temp_row)
  }
  
}

saveRDS(gtf_reduced, 'Data/GTF/GeneCode_ConversionTable.rds')
gtf_reduced <- readRDS('Pancreas/Data/GTF/GeneCode_ConversionTable.rds')


# Creating a column in the result_signif_total table to make the conversion ----
result_signif_total_WithNames <- data.frame()

for (j in 1:dim(result_signif_total)[1]){
  #entry <- gtf_reduced %>% filter(gtf_reduced$id == strsplit(result_signif_total$id.exposure[j], split='\\.')[[1]][1])
  entry <- gtf_reduced[gtf_reduced$id == strsplit(result_signif_total$id.exposure[j], split='\\.')[[1]][1],]
  name <- entry$name
  
  # To fix bug when the ensembleId is not in the GTF Conversion table
  if (dim(entry)[1] == 0){
    print("Not found")
    name <- result_signif_total$id.exposure[j]
  }
  
  temp_name <- cbind(result_signif_total[j,], name)
  result_signif_total_WithNames <- rbind(result_signif_total_WithNames, temp_name) 
  print("found")
}

# Change column Name
names(result_signif_total_WithNames)[names(result_signif_total_WithNames) == "name"] <- "geneName"
result_signif_total <- result_signif_total_WithNames

saveRDS(result_signif_total, 'New/Spleen/Data_MR10/MR_Result/allResultSignifWithNames.rds')

# Forest Plot statistics -------------------------------------------------------

# Calculating the Odd Ratio
result_signif_total$OR <- exp(result_signif_total$b)

# Calculating the confidence interval
result_signif_total$CI_lower <- result_signif_total$b -  result_signif_total$se*qnorm(0.975)
result_signif_total$CI_upper <- result_signif_total$b +  result_signif_total$se*qnorm(0.975)

# Order properly and Create factors for ordering the entries for ggplot
result_signif_total <- result_signif_total[order(result_signif_total$b),]
result_signif_total$geneName <- factor(result_signif_total$geneName,levels=rev(unique(result_signif_total$geneName)))

saveRDS(result_signif_total, 'New/Spleen/Data_MR10/MR_Result/allResults_OR.rds')


### Filtering result_signif_total to get only non-MHC region
filtered_WithNames <- result_signif_total[(result_signif_total$id.exposure %in% filtered$gene_id),]

filtered_WithNames <- filtered_WithNames[order(filtered_WithNames$b),]
filtered_WithNames$geneName <- factor(filtered_WithNames$geneName,levels=rev(unique(filtered_WithNames$geneName)))

saveRDS(filtered_WithNames, 'New/Spleen/Data_MR10/MR_Result/allResults_filteredWithNames.rds')

### Filtering result_signif_total to get only MHC region
filtered_MHC_WithNames <- result_signif_total[(result_signif_total$id.exposure %in% filtered_MHC$gene_id),]

filtered_MHC_WithNames <- filtered_MHC_WithNames[order(filtered_MHC_WithNames$b),]
filtered_MHC_WithNames$geneName <- factor(filtered_MHC_WithNames$geneName,levels=rev(unique(filtered_MHC_WithNames$geneName)))

saveRDS(filtered_MHC_WithNames, 'New/Spleen/Data_MR10/MR_Result/allResults_filteredMHCWithNames.rds')


# Remove bad columns (duplicated geneName.1 for some reason) - (Might not always need)
#result_signif_total <- subset(result_signif_total, select=-c(15))
#filtered_WithNames <- subset(filtered_WithNames, select=-c(15))
#filtered_MHC_WithNames <- subset(filtered_MHC_WithNames, select=-c(15))



# Making a forest plot ---------------------------------------------------------
# 3 plots for total_signif, MHC, and non-MHC plots

data <- result_signif_total
###ALL
p_all <- ggplot(data, aes(x = b, y = geneName, xmin = CI_lower, xmax = CI_upper)) +
  geom_vline(xintercept = 0) +
  geom_point() +
  geom_errorbarh() + 
  xlab("Effect") + 
  ylab("Gene Name") +
  ggtitle("All Gene Significatif")

# Save Plot
pdf('New/Spleen/Data_MR10/Plots/AllSignifForestPlot.pdf', height = 12, width = 10)

#Show Plot
p_all
dev.off()

###MHC
p_MHC <- ggplot(filtered_MHC_WithNames, aes(x = b, y = geneName, xmin = CI_lower, xmax = CI_upper)) +
  geom_vline(xintercept = 0) +
  geom_point() +
  geom_errorbarh() + 
  xlab("Effect") + 
  ylab("Gene Name") +
  ggtitle("All Gene MHC Region Significatif")

# Save Plot
pdf('New/Spleen/Data_MR10/Plots/AllSignif_MHCForestPlot.pdf', height = 12, width = 10)

#Show Plot
p_MHC
dev.off()

### NON-MHC
p_nonMHC <- ggplot(filtered_WithNames, aes(x = b, y = geneName, xmin = CI_lower, xmax = CI_upper)) +
  geom_vline(xintercept = 0) +
  geom_point() +
  geom_errorbarh() + 
  xlab("Effect") + 
  ylab("Gene Name") +
  ggtitle("All Gene Non MHC Significatif") 

# Save Plot
pdf('New/Spleen/Data_MR10/Plots/AllSignif_nonMHCForestPlot.pdf', height = 12, width = 10)

#Show Plot
p_nonMHC
dev.off()

