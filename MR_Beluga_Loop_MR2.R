# Header ---------------------------------------------------------------------

# AUTHOR: John Eric Hamilton
# PURPOSE: Running the Mendelian Randomisation
# NOTES: 
#   *Inputs:   - Exposure and Outcome data 
#
#   *Outputs:  - Preprocessed_Data
#              - Harmonised_Data
#              - Resulting MR files 

# Load libraries and functions -------------------------------------------------
library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)

freadf <- function(df){return(as.data.frame(fread(df)))}

# Load datasets ----------------------------------------------------------------
df2 <- read.delim('T1DGC_saige_EUR.txt')
df <- read.delim("Whole_Blood.v8.independent_eqtls.txt.gz")

exposure <- df

#Filtering for the Fstat with a threshold of 10
#Finding R2 for all values
exposure$R2 <- (get_r_from_bsen(exposure$slope,exposure$slope_se,755))^2 # CHANGE FOR EACH TISSUE

#Finding the Fstat
exposure$Fstat <- (exposure$R2/1)/((1-exposure$R2)/(755-1-1))


#Filtering
exposure <- exposure[exposure$Fstat > 10,]
seuil <- 0.05/length(unique(exposure$gene_id))

# Run MR -----------------------------------------------------------------------
for (i in 1:22){
  
  # Pre-processing exposure
  ## Getting a single chromosome
  eQTL_Pancreas <- exposure %>% filter(grepl(sprintf("chr%s_",i), exposure$variant_id))
  eQTL_Pancreas_Reduced <- select(eQTL_Pancreas, gene_id, variant_id, maf, pval_nominal, slope, slope_se, ref_factor)
  
  # Changes MAF based on value of ref_factor (either -1 or 1)
  eQTL_Pancreas_Reduced[eQTL_Pancreas_Reduced$ref_factor == -1,]$maf <- 1- eQTL_Pancreas_Reduced[eQTL_Pancreas_Reduced$ref_factor == -1,]$maf
  
  ## Renaming and separating columns for harmonize_data() function 
  ### variant_id column separation and merging 
  v_id <- colsplit(eQTL_Pancreas_Reduced$variant_id, "_", c("chromosome", "base_pair_location", "other_allele", "effect_allele", "name"))
  v_id$SNP <- paste(v_id$chromosome, v_id$base_pair_location, sep=":")
  
  test_merge <- cbind(eQTL_Pancreas_Reduced[1], v_id[,c(6,3,4)], eQTL_Pancreas_Reduced[,3:6])
  
  # Remove Indel
  test_merge <- test_merge[test_merge$other_allele%in%c('A','T','C','G') & test_merge$effect_allele%in%c('A','T','C','G'),]
  
  ### Renaming columns for harmonize_data() function
  colnames(test_merge) <- c("gene_id", "SNP", "other_allele.exposure", "effect_allele.exposure", "eaf.exposure","pval.exposure", "beta.exposure", "se.exposure")
  eQTL_Pancreas_Reduced_Renamed <- test_merge 
  
  eQTL_Pancreas_Reduced_Renamed$id.exposure <- eQTL_Pancreas_Reduced_Renamed$gene_id
  eQTL_Pancreas_Reduced_Renamed$exposure <- eQTL_Pancreas_Reduced_Renamed$gene_id
  
  
  # Pre-processing 
  ## Getting a single chromosome
  GCST4023 <- df2 %>% filter(sprintf("%s", i) == df2$Chromosome)
  GCST4023$SNP_Original <- GCST4023$SNP_ID
  
  ## Merging Proxied Data
  proxied <- readRDS("Proxy_GWAS_MichalekEUR_Pancreas_val.rds")
  proxied <- select(proxied, -SNP)
  GCST4023 <- rbind(GCST4023, proxied %>% filter(sprintf("%s", i) == proxied$Chromosome))
  
  ## Getting the columns of interest
  GCST4023_Reduced <- select(GCST4023, c(1,2,3,4,5,7,16,17,18))
  
  ## Renaming and Separating columns 
  ### Making SNP column
  GCST4023_Reduced$SNP <- paste(sprintf("chr%s",GCST4023_Reduced$Chromosome), GCST4023_Reduced$Position, sep=":")
  
  # Remove Indel
  GCST4023_Reduced <- GCST4023_Reduced[GCST4023_Reduced$Allele1%in%c('A','T','C','G') & GCST4023_Reduced$Allele2%in%c('A','T','C','G'),]
  
  ### Renaming columns for harmonize_data() function
  colnames(GCST4023_Reduced) <- c("SNP_ID", "chromosome", "base_pair_location", "other_allele.outcome", "effect_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "SNP")
  GCST4023_Reduced_Renamed <- GCST4023_Reduced
  
  ### Adding outcome column
  GCST4023_Reduced_Renamed$outcome <- "T1D" 
  GCST4023_Reduced_Renamed$id.outcome <- "T1D"
  
  saveRDS(eQTL_Pancreas_Reduced_Renamed, sprintf('eQTL_Pancreas_chr%s.rds', i))
  saveRDS(GCST4023_Reduced_Renamed, sprintf('GCST4023_chr%s.rds', i))
  
  
  # Harmonise Data
  harmonised <- harmonise_data(exposure_dat = eQTL_Pancreas_Reduced_Renamed, outcome_dat = GCST4023_Reduced_Renamed, action = 2)
  
  saveRDS(harmonised, sprintf('Har_Chr%s.rds', i))
  
  # Mendelian Randomization
  result <- mr(harmonised)
  saveRDS(result, sprintf('Raw_MR_Chr%s.rds', i))
  result_signif <- result %>% filter(result$pval < seuil) 
  
  saveRDS(result, sprintf('MR_Chr%s.rds', i))
}