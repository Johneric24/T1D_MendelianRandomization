# Header ----------------------------------------------------------------------

# AUTHOR: John Eric Hamilton
# PURPOSE: Performing colocalization study (inspired by code from Basile 02 BJ Colocalization),
#          looping through all significant genes
# NOTES: 
#   *Inputs:   - MR Results filtered and non-filtered (to choose one gene at a time)
#              - Exposition file 
#              - GWAS (Chiou Only)
#
#   *Outputs:  - Colocalization test statistics  
#              - SNP plot with p-value

# Loading Library and Data -----------------------------------------------------
library(arrow)
library(coloc)
library(susieR)
library(gridExtra)
library(tidyverse)
library(LDlinkR) #To get the LD
library('reshape') #Formating data


# Initialize result table
tab_res <- data.frame()

# Significant results that were not in MHC region
signif_total <- readRDS('New/Spleen/Data_MR8/MR_Result/allResults_filteredWithNames.rds')
signif_total <- signif_total[signif_total$chromosome == 6,]

# Loop -------------------------------------------------------------------------
# Taking the gene TH to start (eventually will loop through all of them)
for (i in 1:length(signif_total$geneName)){
  # Gene to coloc analyze information
  ensembleID <- signif_total$id.exposure[i]
  prot <- signif_total$geneName[i]
  chrI <- signif_total$chromosome[i]
  entry <- df[grep(ensembleID, df$gene_id),]
  position <-colsplit(entry$variant_id, "_", c("chromosome", "base_pair_location", "other_allele", "effect_allele", "name"))[2]
  SNPI <- paste0("chr", chrI, ":", position[1,1])
  
  # 1000000 Mb windows
  pos_min <- position[1,1] - 500000
  pos_max <- position[1,1] + 500000
  
  ### EXPO FILE ###
  expo_combined <- readRDS(sprintf('Coloc/Spleen/Spleen_allPairs_chr%s.rds', chrI))
  
  # cut the Exposure for the given base pair range
  expo_cut <- expo_combined[expo_combined$base_pair_location <= pos_max & expo_combined$base_pair_location >= pos_min, ]
  
  # Remove Indel
  expo_cut = expo_cut[expo_cut$other_allele%in%c('A','T','C','G') & expo_cut$effect_allele%in%c('A','T','C','G'),]
  
  ### OUTCOME FILE ###
  outc_interest <- readRDS(sprintf('Pancreas/Data/Proxy/ChiouWithSNP_chr%s.rds', chrI))
  
  proxied <- readRDS('Pancreas/Data/Proxy/outcomeToMerge.rds')
  outc_interest <- rbind(outc_interest, proxied %>% filter(sprintf("%s", i) == proxied$chromosome))
  
  
  # Remove Indel
  outc_interest=outc_interest[outc_interest$other_allele%in%c('A','T','C','G') & outc_interest$effect_allele%in%c('A','T','C','G'),]
  
  # Allele mineur rare (to remove some of the duplicates in the outcome)
  outc_interest <- outc_interest[outc_interest$effect_allele_frequency > 0.01,]
  
  # merge the 2 GWAS using harmonise_data() function 
  colnames(expo_cut) <- c("gene_id", "chromosome", "base_pair_location", "other_allele.exposure", "effect_allele.exposure", "tss_distance", "ma_samples", "ma_count", "eaf.exposure", "pval.exposure", "beta.exposure", "se.exposure", "SNP")
  expo_cut_Renamed <- expo_cut[expo_cut$gene_id == ensembleID,]
  expo_cut_Renamed$id.exposure <- ensembleID
  expo_cut_Renamed$exposure <- ensembleID
  
  
  colnames(outc_interest) <- c("variant_id", "pval.outcome", "chromosome", "base_pair_location", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "sample_size", "SNP")
  outc_interest_Renamed <- outc_interest
  outc_interest_Renamed$id.outcome <- "T1D"
  outc_interest_Renamed$outcome <- "T1D"
  
  harm_data <- harmonise_data(expo_cut_Renamed, outc_interest_Renamed)
  saveRDS(harm_data, sprintf('New/Spleen/Coloc_Toute_MR8/Genes/Coloc_Harm_%s_result.rds', prot))
  
  
  # Remove duplicated SNPs that still remain
  test <- duplicated(harm_data$SNP)
  SNP_toRemove <- harm_data[test,"SNP"]
  harm_data <- harm_data[!(harm_data$SNP %in% SNP_toRemove),]
  
  # Rearrange data
  tmp <- harm_data[, c(1, 29, 16, 17, 6, 28, 27, 25, 8, 7, 18, 15)]
  
  colnames(tmp) <- c("SNP", "gene_id", "chr", "position", "expo_eff", "expo_se", "expo_pv", "expo_N", "expo_maf",  
                     "outc_eff", "outc_se", "outc_pv")
  # varbeta = square(se)
  tmp$expo_varbeta <- (tmp$expo_se)^2
  tmp$outc_varbeta <- (tmp$outc_se)^2
  
  # COLOC ANALYSIS
  # 2 dataset for analyse Prot - ANM
  D1_expo <- list() 
  D2_outc <- list()
  
  # expo
  D1_expo$type <- "quant"
  D1_expo$snp <- tmp$SNP
  D1_expo$position <- tmp$position
  D1_expo$beta <- tmp$expo_eff
  D1_expo$varbeta <- tmp$expo_varbeta
  D1_expo$varbeta[is.na(D1_expo$varbeta)] <- 1
  D1_expo$sdY <- 1 # (prot level are scale)
  
  
  # outc
  D2_outc$type <- "cc"
  D2_outc$snp <- tmp$SNP
  D2_outc$position <- tmp$position
  D2_outc$beta <- tmp$outc_eff
  D2_outc$varbeta <- tmp$outc_varbeta
  D2_outc$sdY <- 1 # (FROM GWAS STUDY : 1.3 for AAM and 4 for ANM)
  
  
  saveRDS(D1_expo, sprintf('New/Spleen/Coloc_Toute_MR8/Genes/Coloc_Expo_%s_result.rds', prot))
  saveRDS(D2_outc, sprintf('New/Spleen/Coloc_Toute_MR8/Genes/Coloc_Outc_%s_result.rds', prot))
  
  # Run coloc analysis 
  res <- coloc.abf(dataset1 = D1_expo, dataset2 = D2_outc)
  
  # save the result
  tab_res <- rbind(tab_res, c(paste(prot), as.vector(res$summary)))
  
}

# Format and Save Results
colnames(tab_res) <- c("Gene Name", "nSNP", "H0", "H1", "H2", "H3", "H4")
tab_res$H0 <- as.numeric(tab_res$H0)
tab_res$H1 <- as.numeric(tab_res$H1)
tab_res$H2 <- as.numeric(tab_res$H2)
tab_res$H3 <- as.numeric(tab_res$H3)
tab_res$H4 <- as.numeric(tab_res$H4)

saveRDS(tab_res, 'New/Spleen/Coloc_Toute_MR8/Coloc_result_chrome6.rds')

# Filter for genes that passed colocalization 
final_res <- tab_res[tab_res$H4 > 0.8,] 
final_moreAnalysis <- tab_res[!(tab_res$H4 > 0.8),]

saveRDS(final_res,'New/Spleen/Coloc_Toute_MR8/passingColocGenes.rds')
saveRDS(final_moreAnalysis,'New/Spleen/Coloc_Toute_MR8/ColocGenesToProxy.rds')


write.table(final_res$`Gene Name`, 'New/Spleen/Coloc_Toute_MR8/PassingGeneNames.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)

final_res <- readRDS('New/Spleen/Coloc_Toute_MR8/passingColocGenes.rds')
final_moreAnalysis <- readRDS('New/Spleen/Coloc_Toute_MR8/ColocGenesToProxy.rds')
final_all <- rbind(final_res, final_moreAnalysis)

#Filter for genes that might need proxy (High H3 and very low H4)
final_moreAnalysis_Proxy <- final_moreAnalysis[final_moreAnalysis$H3 > 0.8,]

