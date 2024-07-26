# Header ---------------------------------------------------------------------

# AUTHOR: John Eric Hamilton
# PURPOSE: Replication eQTLgen dataset
# NOTES: 
#   *Inputs:   - MR results 
#              - eQTLgen dataset
#
#   *Outputs:  - MR Results with exposure as the eQTLgen dataset
#              - Searched MAF files

# Load library -----------------------------------------------------------------
library(rsnps)
library(TwoSampleMR)

# Loading in the data and Preprocessing ----------------------------------------
df_new <- read.delim("2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")

allGenes <- read.delim("allGenes.txt", header = FALSE)

eQTL <- df_new[df_new$GeneSymbol%in%allGenes$V1,]

# Genes that were not found in the Expression data
write.table(allGenes[!allGenes$V1%in%eQTL$GeneSymbol,], 'genesToProxyNewEQTL.txt', row.names = FALSE, quote = FALSE, col.names = FALSE)


# Calculating the values we need from the zscore and pval given
eQTL$Beta <- eQTL$Zscore/abs(qnorm(eQTL$Pvalue/2))
eQTL$SE <- eQTL$Beta/eQTL$Zscore

# MAF search -------------------------------------------------------------------
SNP_done <- data.frame()
# Assessed allele is the effect allele 
# 0.35 maf default then find SNP for each using ncbi_snp_query ($maf_population[[1]]) use dbGaP_PopFreq
allSNPs <- eQTL$SNP

#outcomeToMerge3 is the updated eQTL file with EAF
outcomeToMerge3 <- eQTL
toSkip <- data.frame()

for(i in 1:nrow(outcomeToMerge3)){
  tryCatch({
    t=as.data.frame(ncbi_snp_query(as.character(outcomeToMerge3[i,'SNP'])))$maf_population[[1]]
  t=t[t$study=='dbGaP_PopFreq',]
  
  t2=outcomeToMerge3[i,]
  t3=t[(t$ref_seq==t2$AssessedAllele & t$Minor==t2$OtherAllele) | (t$Minor==t2$AssessedAllele & t$ref_seq==t2$OtherAllele),]
  
  if(t3$Minor==t2$AssessedAllele){
    outcomeToMerge3[i,'EAF']=t3$MAF
  }else{
    outcomeToMerge3[i,'EAF']=1-t3$MAF
  }
  print(t3$MAF)
  print(i)},
  error = function(e) {
    toSkip <- rbind(toSkip, outcomeToMerge3[i,'SNP'])
    }
  )
}
saveRDS(outcomeToMerge3, 'eQTLgen/eQTL_searchedMAF.rds')

# Retried MAF search so now the file is called outcomeToMerge_updated
outcomeToMerge_updated <- readRDS('eQTLgen/eQTL_searchedMAF_updated.rds')

# Clumping (LD) ----------------------------------------------------------------
# Formatting the exposure data for clumping
formexposure01 <- format_data(
  dat = outcomeToMerge_updated,
  type = "exposure",
  snps = NULL, #leave it as is
  header = TRUE, #leave it as is
  phenotype_col = "Phenotype", #leave it as is
  snp_col = "SNP",
  beta_col = "Beta",
  se_col = "SE",
  eaf_col = "EAF", 
  effect_allele_col = "AssessedAllele",
  other_allele_col = "OtherAllele",
  pval_col = "Pvalue",
  units_col = "units", #leave it as is
  ncase_col = "ncase", #leave it as is
  ncontrol_col = "ncontrol", #leave it as is
  samplesize_col = "NrSamples",
  gene_col = "GeneSymbol", 
  id_col = "Gene", 
  min_pval = 1e-200, #leave it as is
  z_col = "z", #leave it as is
  info_col = "info", #leave it as is
  chr_col = "SNPChr",
  pos_col = "GenePos",
  log_pval = FALSE) #leave it as is

# Setting to remove a reoccuring error
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

#Clumping
clumped_exposure01<-clump_data(
  dat = formexposure01, 
  clump_kb = 10000, #should be long enough
  clump_r2 = 0.001, #leave it as is
  clump_p1 = 1, #leave it as is
  clump_p2 = 1, #leave it as is
  pop = "EUR",
  bfile = NULL, #leave it as is
  plink_bin = NULL #leave it as is
)

#Identifying the exposure phenotype
clumped_exposure01$exposure <- clumped_exposure01$id.exposure


# MR ---------------------------------------------------------------------------
for (i in 1:22){
  
  # Pre-processing 
  ## Getting a single chromosome
  eQTL_WB <- clumped_exposure01 %>% filter(clumped_exposure01$chr.exposure == i)
  saveRDS(eQTL_WB, sprintf('eQTLgen/MAF/Data_MR4/Preprocessed_Data/eQTL_WB_chr%s.rds', i))
  
  GCST4023_Reduced_Renamed <- readRDS(sprintf('New/Pancreas/Data_MR4/Preprocessed_Data/GCST4023_chr%s.rds', i))
  
  #colnames(GCST4023_Reduced_Renamed)[c(1, 2)] <- c("SNP", "Location") # For Chiou GWAS
  colnames(GCST4023_Reduced_Renamed)[c(1,10)] <- c("SNP", "Location")  # For Michalek GWAS
  
  
  # Harmonise Data
  harmonised <- harmonise_data(exposure_dat = eQTL_WB, outcome_dat = GCST4023_Reduced_Renamed, action = 2)
  
  saveRDS(harmonised, sprintf('eQTLgen/MAF/Data_MR4/Harmonised_Data/Har_Chr%s.rds', i))
  
  # Mendelian Randomization
  result <- mr(harmonised)
  saveRDS(result, sprintf('eQTLgen/MAF/Data_MR4/MR_Result/Raw_MR_Chr%s.rds', i))
}


# Data Analysis ----------------------------------------------------------------
seuil <- 1

# Importing the results for chromosome 1 to 22 and merging all significant genes in one table
result_signif_total <- data.frame()
for (i in 1:22){
  result_temp <- readRDS(sprintf('eQTLgen/MAF/Data_MR4/MR_Result/Raw_MR_Chr%s.rds', i))
  
  if (dim(result_temp %>% filter(result_temp$pval < seuil))[1] > 0){
    result_signif <- result_temp %>% filter(result_temp$pval < seuil)
    result_signif$chromosome <- sprintf("%s", i)
    
    result_signif_total <- rbind(result_signif_total, result_signif)
  }
  
  #result_temp$chromosome <- sprintf("%s", i)
  #result_all <- rbind(result_all, result_temp)
}

saveRDS(result_signif_total,'eQTLgen/MAF/Data_MR4/MR_Result/allResult_allNonSignifAndSignif.rds')


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

saveRDS(result_signif_total, 'eQTLgen/MAF/Data_MR4/MR_Result/allResultSignifWithNames.rds')

# Checking Replication ---------------------------------------------------------
seuil_eQTL <- 0.05/50

MR1 <- readRDS('eQTLgen/Data_MR1/MR_Result/allResultSignifWithNames.rds')
MR1_filtered <- MR1[MR1$pval < seuil_eQTL,]
write.table(MR1_filtered[!duplicated(MR1_filtered$geneName),]$geneName, 'eQTLgen/ChiouGenes.txt', quote = FALSE, row.names=FALSE, col.names=FALSE)


MR2 <- readRDS('eQTLgen/Data_MR2/MR_Result/allResultSignifWithNames.rds')
MR2_filtered <- MR2[MR2$pval < seuil_eQTL,]
write.table(MR2_filtered[!duplicated(MR2_filtered$geneName),]$geneName, 'eQTLgen/MichalekEURGenes.txt', quote = FALSE, row.names=FALSE, col.names=FALSE)


MR4 <- readRDS('eQTLgen/Data_MR4/MR_Result/allResultSignifWithNames.rds')
MR4_filtered <- MR4[MR4$pval < seuil_eQTL,]
write.table(MR4_filtered[!duplicated(MR4_filtered$geneName),]$geneName, 'eQTLgen/MichalekTransethnicGenes.txt', quote = FALSE, row.names=FALSE, col.names=FALSE)

