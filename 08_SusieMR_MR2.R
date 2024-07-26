# Header ----------------------------------------------------------------------

# AUTHOR: John Eric Hamilton
# PURPOSE: Performing Susie colocalization study 
#
# NOTES: 
#   *Inputs:   - LD Matrix
#              - Harm Data from Exposition + Outcome of Michalek only 
#              - Plink BED file
#
#   *Outputs:  - Colocalization statistics  
#              

# Prepare Data to be used for co localization (RSID) ---------------------------  

# Significant results that were not in MHC region (That didnt pass coloc that we want to do Susie Coloc with)
signif_total <- readRDS('New/Pancreas/Data_MR2/MR_Result/allResults_filteredWithNames.rds')

final_moreAnalysis <- readRDS('New/Pancreas/Coloc_Toute_MR2/ColocGenesToProxy.rds')
signif_total <- signif_total[(signif_total$geneName %in% final_moreAnalysis$`Gene Name`),]

# To generate the new tmp files to be used with Plink
for (i in 1:length(signif_total$geneName)){
  
  # Get information of gene of interest
  # Taking the gene TH to start (eventually will loop through all of them)
  ensembleID <- signif_total$id.exposure[i]
  prot <- signif_total$geneName[i]
  chrI <- signif_total$chromosome[i]
  entry <- df[grep(ensembleID, df$gene_id),]
  position <-colsplit(entry$variant_id, "_", c("chromosome", "base_pair_location", "other_allele", "effect_allele", "name"))[2]
  SNPI <- paste0("chr", chrI, ":", position[1,1])
  
  
  # Get the exposition file
  expo_combined <- readRDS(sprintf('Coloc/Spleen/Spleen_allPairs_chr%s.rds', chrI))
  
  # 1000000 Mb windows
  pos_min <- position[1,1] - 500000
  pos_max <- position[1,1] + 500000
  
  # cut the Exposure for the given base pair range
  expo_cut <- expo_combined[expo_combined$base_pair_location <= pos_max & expo_combined$base_pair_location >= pos_min, ]
  
  # remove Indel
  expo_cut = expo_cut[expo_cut$other_allele%in%c('A','T','C','G') & expo_cut$effect_allele%in%c('A','T','C','G'),]
  
  ### OUTCOME FILE ###
  outc_interest <- readRDS(sprintf('Pancreas/Data_MR4/Proxy/MichalekTransethnicWithSNP_chr%s.rds', chrI))
  
  # Remove Indel
  outc_interest=outc_interest[outc_interest$Allele1%in%c('A','T','C','G') & outc_interest$EA%in%c('A','T','C','G'),]
  
  # Allele mineur rare (to remove some of the duplicates in the outcome)
  #outc_interest <- outc_interest[outc_interest$effect_allele_frequency > 0.01,]
  
  # merge the 2 GWAS using harmonise_data() function 
  colnames(expo_cut) <- c("gene_id", "chromosome", "base_pair_location", "other_allele.exposure", "effect_allele.exposure", "tss_distance", "ma_samples", "ma_count", "eaf.exposure", "pval.exposure", "beta.exposure", "se.exposure", "SNP")
  
  expo_cut_Renamed <- expo_cut[expo_cut$gene_id == ensembleID,]
  expo_cut_Renamed$id.exposure <- ensembleID
  expo_cut_Renamed$exposure <- ensembleID
  
  colnames(outc_interest) <- c("variant_id", "chromosome", "base_pair_location", "other_allele.outcome", "effect_allele.outcome", "EA", "eaf.outcome", "sample_size", "EAF_case", "Case_size", "EAF_control", "Control_size", "OR", "LCI", "UCI", "beta.outcome", "se.outcome", "pval.outcome", "SNP")
  
  outc_interest_Renamed <- outc_interest
  outc_interest_Renamed$id.outcome <- "T1D"
  outc_interest_Renamed$outcome <- "T1D"
  
  
  harm_data <- harmonise_data(expo_cut_Renamed, outc_interest_Renamed)
  
  # Remove duplicated SNPs that still remain
  test <- duplicated(harm_data$SNP)
  SNP_toRemove <- harm_data[test,"SNP"]
  harm_data <- harm_data[!(harm_data$SNP %in% SNP_toRemove),]
  
  # rearrange data
  tmp <- harm_data[, c(1, 29, 15, 16, 6, 36, 35, 33, 8, 7, 26, 27)]
  
  colnames(tmp) <- c("SNP", "gene_id", "chr", "position", "expo_eff", "expo_se", "expo_pv", "expo_N", "expo_maf",  
                     "outc_eff", "outc_se", "outc_pv")
  
  # varbeta = square(se)
  tmp$expo_varbeta <- (tmp$expo_se)^2
  tmp$outc_varbeta <- (tmp$outc_se)^2
  
  ### AT THIS POINT TMP FILES ARE GENERATED ###
  
  # Get RS ID using outcome file and add that column in tmp
  outc_RSID <- readRDS(sprintf('New/Spleen/Data_MR10/Preprocessed_Data/GCST4023_chr%s.rds', chrI))
  
  # Remove Indel
  outc_RSID=outc_RSID[outc_RSID$other_allele.outcome%in%c('A','T','C','G') & outc_RSID$effect_allele.outcome%in%c('A','T','C','G'),]
  
  # Identify Similar SNPs 
  outc_RSID_row <- outc_RSID[outc_RSID$SNP%in%tmp$SNP,]
  outc_RSID_row <- outc_RSID_row %>% select("SNP_ID", "SNP")
  outc_RSID_row <- unique(outc_RSID_row)
  
  # Merge by SNP to add the RSID column to the data file 
  tmp_new <- merge(harm_data, outc_RSID_row, by.x="SNP", by.y="SNP")
  saveRDS(tmp_new, sprintf("New/Spleen/Plink10/tmp/tmp_Chr%s_%s.rds", chrI, prot))
  write.table(tmp_new$SNP_ID, sprintf('New/Spleen/Plink10/testChr%s_%s.txt', chrI, prot), row.names = FALSE, col.names = FALSE, quote = FALSE) #Use in Plink
}

# Plink ------------------------------------------------------------------------
# Open Ubuntu (Linux Terminal)

# filepath /mnt/c/Users/John/Documents/CHU_Stage2024/Summer_2024/ExpressionGWAS/Plink
# command ./plink.exe --bfile ../Plink/Chr11/Chr11.EUR --extract testChr11.txt --recodeA --out susieChr11


# Function to Run Susie -----------------------------------------------

signif_total <- signif_total[!(duplicated(signif_total$geneName)),]
susieAllResults <- NULL

for (i in 1:length(signif_total$geneName)){
  ensembleID <- signif_total$id.exposure[i]
  prot <- signif_total$geneName[i]
  chrI <- signif_total$chromosome[i]
  entry <- df[grep(ensembleID, df$gene_id),]
  position <-colsplit(entry$variant_id, "_", c("chromosome", "base_pair_location", "other_allele", "effect_allele", "name"))[2]
  SNPI <- paste0("chr", chrI, ":", position[1,1])
  
  
  # And after estimate LD with R
  plinkLD_raw <- function(raw_files = sprintf("New/Spleen/Plink10/susieChr%s_%s.raw", chrI, prot)) {
    #' fn : imput mean for NA data
    #'
    #' @param x a matrix
    #'
    #' @export
    fn <- function(x){
      m. <- mean(x, na.rm = TRUE)
      x[is.na(x)] <- m.
      return(x)
    }
    t1 <- read.table(raw_files, header = T)
    ld1 <- as.matrix(t1[, c(7:ncol(t1))])
    # number of NA
    nNA <- sum(is.na(ld1))
    nNA <- nNA / length(as.vector(ld1)) * 100
    # remplace NA by colmeans
    ld1 <- fn(ld1)
    # add intercept (just a trick to avaid null ecart type)
    ld1 <- rbind(rep(0.5, ncol(ld1)),
                 ld1)
    # R value (correlation MATRIX)
    ld1 <- cor(ld1)
    # issues with the rsid (need remove _A, _T, ...)
    rs <- data.frame(rs = colnames(ld1))
    rs <- tidyr::separate(rs, rs, c("rs", "other"), "_")
    colnames(ld1) <- rs$rs
    row.names(ld1) <- rs$rs
    return(list(ld_r = ld1,
                out_plink = t1,
                nNA = paste0(round(nNA, 4), " % of NAs")))
  }
  
  # Define Function to remove NA values
  fn <- function(x){
    
    m. <- mean(x, na.rm = TRUE)
    x[is.na(x)] <- m.
    return(x)
  }
  
  # Loading the file from plink into RStudio
  tmps <- plinkLD_raw(sprintf("New/Spleen/Plink10/susieChr%s_%s.raw", chrI, prot))
  
  
  # Making the LD matrix ---------------------------------------------------------
  ld_r <- tmps$ld_r
  
  # For outc_plink
  outc_plink <- tmps$out_plink
  outc_plink <- data.frame(rsid = colnames(outc_plink)[7:ncol(outc_plink)])
  outc_plink <- separate(outc_plink, rsid, c("rsid", "eff_all"), "_")
  colnames(outc_plink) <- c("rsid", "ld_eff_all")
  outc_plink <- outc_plink[outc_plink$rsid %in% colnames(tmps$ld_r), ]
  
  tmp_new <- readRDS(sprintf("New/Spleen/Plink10/tmp/tmp_Chr%s_%s.rds", chrI, prot))
  # Merge outc_plink with tmp file that has the RSIDs
  tmp_new_susie <- merge.data.frame(outc_plink, tmp_new, by.x = 1, by.y=43)
  
  # remove mismatch between eff allele expo, ouctome and ld matrix (metabol = effect_allele)
  dat_susie <- tmp_new_susie[(tmp_new_susie$ld_eff_all == tmp_new_susie$effect_allele.exposure) & (tmp_new_susie$ld_eff_all ==   tmp_new_susie$effect_allele.outcome), ]
  
  # ld with the right set of rsid
  ld_r <- ld_r[dat_susie$rsid, dat_susie$rsid] #LD MATRIX
  
  saveRDS(dat_susie, sprintf('New/Spleen/Plink10/susieResults/dat_susieChr%s_%s.rds', chrI, prot))
  saveRDS(ld_r, sprintf('New/Spleen/Plink10/susieResults/LD_MatrixChr%s_%s.rds', chrI, prot))
  
  # calculate Z-score for exposure data and outcome data
  dat_susie$z_expo <- dat_susie$beta.exposure / dat_susie$se.exposure
  dat_susie$z_outc <- dat_susie$beta.outcome / dat_susie$se.outcome
  
  colnames(dat_susie) <- c("rsid", "ld_eff_all", "snp", "effect_allele", "other_allele", "effect_allele.outcome", "other_allele.outcome", "beta.exposure", "beta.outcome", "eaf.exposure", "eaf.outcome.", "remove", "palindromic", "ambiguous", "id.outcome","variant_id.x", "chromosome", "base_pair_location", "effect.allele" ,"sample_size", "EAF_case", "Case_size", "EAF_control", "Control_size", "OR", "LCI", "UCI", "se.outcome", "pval.outcome", "outcome", "gene_id", "chromosome.y", "base_pair_location.y", "tss_distance", "ma_samples", "ma_count", "pval.exposure",  "se.exposure", "id.exposure", "exposure", "action", "SNP_index", "mr_keep", "samplesize.outcome", "z_expo", "z_outc")
  
  # Check for NA
  ld_r_na <- fn(ld_r)
  
  # calculate lambda
  # concordance between the LD matrix and the z-scores (should be near 0)
  lambda_expo <- estimate_s_rss(dat_susie$z_expo, (ld_r_na)^2, n = mean(dat_susie$ma_samples))
  lambda_outc <- estimate_s_rss(dat_susie$z_outc, (ld_r_na)^2, n = mean(dat_susie$sample_size))
  
  # use susis_rss same as runsusie (runsusie add some names)
  fit_expo <- susie_rss(dat_susie$z_expo, (ld_r_na)^2, n = mean(dat_susie$ma_samples))
  fit_outc <- susie_rss(dat_susie$z_outc, (ld_r_na)^2, n = mean(dat_susie$sample_size))
  
  
  # Coloc susie ------------------------------------------------------------------
  res <- coloc.susie(fit_expo, fit_outc) 
  
  #View(res$summary)
  
  
  saveRDS(res, sprintf("New/Spleen/Plink10/susieResults/coloc_susie_Chr%s_%s.rds", chrI, prot))
  saveRDS(res$summary , sprintf("New/Spleen/Plink10/susieResults/susieResults_%s_%s.rds", chrI, prot))
  
  susieSignifResults <- res$summary[res$summary$PP.H4.abf > 0.8, ]
  susieAllResults <- rbind(susieAllResults, susieSignifResults)
  print("Next one...")
}

# Summary Table
susieAllResults_ <- data.frame()
for (i in 1:length(signif_total$geneName)){
  ensembleID <- signif_total$id.exposure[i]
  prot <- signif_total$geneName[i]
  chrI <- signif_total$chromosome[i]
  entry <- df[grep(ensembleID, df$gene_id),]
  position <-colsplit(entry$variant_id, "_", c("chromosome", "base_pair_location", "other_allele", "effect_allele", "name"))[2]
  SNPI <- paste0("chr", chrI, ":", position[1,1])
  
  
  res <- readRDS(sprintf("New/Spleen/Plink10/susieResults/coloc_susie_Chr%s_%s.rds", chrI, prot))
  
  if (length(res) != 1){
    susieSignifResults <- res$summary[res$summary$PP.H4.abf > 0.8, ]
    susieSignifResults$Gene <- paste(prot)
    susieAllResults_ <- rbind(susieAllResults_, susieSignifResults)
  }
}

saveRDS(susieAllResults_, "New/Spleen/Plink9/susieResults/allSignificantSusies.rds")


