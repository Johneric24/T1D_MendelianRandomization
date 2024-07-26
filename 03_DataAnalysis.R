# Header ---------------------------------------------------------------------

# AUTHOR: John Eric Hamilton
# PURPOSE: Interpretting the results of MR
# NOTES: 
#   *Inputs:   - MR Results 
#              - Seuil (Bonferroni calculated by dividing 0.05 by the number of entries in the exposure)
#
#   *Outputs:  - All significant resulting SNPs from the MR (called result_signif_total)
#              - HLA/MHC region filtered significant SNPs (called filtered_MHC)
#              - non-HLA region filtered significant SNPs (called filtered) 
#   *Resulting files contain merged chromosome results


# Generating Results ----------------------------------------------------------

# Recall the Bonferroni threshold (eQTL Pancreas)
seuil <- 1

# Importing the results for chromosome 1 to 22 and merging all significant genes in one table
result_signif_total <- data.frame()
result_all <- data.frame()
for (i in 1:22){
  result_temp <- readRDS(sprintf('New/Pancreas/Data_MR1/MR_Result/Raw_MR_Chr%s.rds', i))
  
  if (dim(result_temp %>% filter(result_temp$pval < seuil))[1] > 0){
    result_signif <- result_temp %>% filter(result_temp$pval < seuil)
    result_signif$chromosome <- sprintf("%s", i)
    
    result_signif_total <- rbind(result_signif_total, result_signif)
  }
  
  result_temp$chromosome <- sprintf("%s", i)
  result_all <- rbind(result_all, result_temp)
}

saveRDS(result_all,'New/Spleen/Data_MR10/MR_Result/allResult_allNonSignifAndSignif.rds')


## Filtering genes that are in the MHC region (chromosome 6)
filtered <- data.frame()
filtered_MHC <- data.frame()

for (j in 1:length(result_signif_total$id.exposure)){
  temp <- df %>% filter(df$gene_id == result_signif_total[[1]][j])
  
  chromosome <- colsplit(temp$variant_id, "_", c("chromosome", "base_pair_location", "other_allele", "effect_allele", "name"))[[1]] 
  location <- colsplit(temp$variant_id, "_", c("chromosome", "base_pair_location", "other_allele", "effect_allele", "name"))[[2]]
  
  if (length(location) > 1){
    location <- location[1]
  }
  
  if (length(chromosome) > 1){
    chromosome <- chromosome[1]
  }
  
  if (chromosome == "chr6"){
    if (between(location, 28510120, 33480577)){
      filtered_MHC <- rbind(filtered_MHC, temp)
    } else {
      # Could be chromosome 6 but not MHC region
      filtered <- rbind(filtered, temp)
    }
  } else {
    filtered <- rbind(filtered, temp)
  }
}


saveRDS(result_signif_total, 'New/Pancreas/Data_MR2/MR_Result/allResultSignif.rds')
saveRDS(filtered, 'New/Spleen/Data_MR10/MR_Result/allResults_filtered.rds')
saveRDS(filtered_MHC, 'New/Spleen/Data_MR10/MR_Result/allResults_filtered_MHC_region.rds')
