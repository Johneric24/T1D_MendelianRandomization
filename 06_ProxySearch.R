# Header ----------------------------------------------------------------------

# AUTHOR: John Eric Hamilton
# PURPOSE: Proxy Search for SNPs that are in the exposure and are missing the outcome
# NOTES: 
#   *Inputs:   - Exposure (eQTL)
#              - Outcome (df3_GCST4023)
#              - Proxy Data from ProxyBatch.R
#
#   *Outputs:  - outcomeToMerge with the processed proxied SNPs found (to be added to the outcome) 
            

# Loading Data and looking for missing SNPs in outcome ---------------------------------
### GWAS Chiou 'df3_GCST4023' (outcome) vs eQTL 'df' (exposure)
eQTL <- data.frame()
for (i in 1:22){
  eQTL <- rbind(eQTL, readRDS(sprintf('Pancreas/Data/Preprocessed_Data/eQTL_Pancreas_chr%s.rds', i)))
}
proxy_MR1 <- eQTL[!eQTL$SNP%in%df3_GCST4023$SNP,]
saveRDS(proxy_MR1, 'Data/Proxy/missingFromOutcome.rds')


# Prepping Chiou GWAS for LD Proxy ------------------------------------------------------
# Adding SNP column
df3_GCST4023$SNP <- paste("chr", df3_GCST4023$chromosome, ":", df3_GCST4023$base_pair_location, sep="")


# Remove entry if there is multiple alleles in either the other_allele or effect_allele (indels)
df3_GCST4023 <- df3_GCST4023[df3_GCST4023$other_allele.outcome%in%c('A','T','C','G') & df3_GCST4023$effect_allele.outcome%in%c('A','T','C','G'),]
# Save
saveRDS(df3_GCST4023, 'Data/Proxy/ChiouWithSNP.rds')

# Prepping Michalek EUR GWAS for LD Proxy ------------------------------------------------------
# Adding SNP column
df2$SNP <- paste("chr", df2$Chromosome, ":", df2$Position, sep="")


# Remove entry if there is multiple alleles in either the other_allele (Allele1) or effect_allele (Allele2) (indels)
df2 <- df2[df2$Allele1%in%c('A','T','C','G') & df2$EA%in%c('A','T','C','G'),]

# Save
saveRDS(df2, 'Data_MR2/Proxy/MichalekWithSNP.rds')

df2 <- readRDS('Pancreas/Data_MR2/Proxy/MichalekWithSNP.rds')

proxy_MR2 <- eQTL[!(eQTL$SNP%in%df2$SNP),]

saveRDS(proxy_MR2, 'missingFromOutcome.rds')


# Prepping Michalek AMR GWAS for LD Proxy ------------------------------------------------------

dfAM <- dfAM[dfAM$Allele1%in%c('A','T','C','G') & dfAM$EA%in%c('A','T','C','G'),]

# Save
saveRDS(dfAM, 'Data_MR3/Proxy/MichalekAMRWithSNP.rds')

dfAM <- readRDS('Data_MR3/Proxy/MichalekAMRWithSNP.rds')

proxy_MR3 <- eQTL[!(eQTL$SNP%in%dfAM$SNP),]

saveRDS(proxy_MR3, 'Data_MR3/Proxy/missingFromOutcome.rds')


# Prepping Michalek Transethnic GWAS for LD Proxy ------------------------------
df4$SNP <- paste("chr", df4$Chromosome, ":", df4$Position, sep="")

# Remove entry if there is multiple alleles in either the other_allele (Allele1) or effect_allele (Allele2) (indels)
df4 <- df4[df4$Allele1%in%c('A','T','C','G') & df4$EA%in%c('A','T','C','G'),]

# Save
saveRDS(df4, 'Data_MR4/Proxy/MichalekTransethnicWithSNP.rds')

proxy_MR4 <- eQTL[!(eQTL$SNP%in%df4$SNP),]
saveRDS(proxy_MR4, 'Data_MR4/Proxy/missingFromOutcome.rds')


# Troubleshooting LD Proxy Errors ---------------------------------------------- 
proxy_MR1 <- readRDS('Data/Proxy/missingFromOutcome.rds')
df3_GCST4023 <- readRDS('Data/Proxy/ChiouWithSNP.rds')

# Removing the ones that were skipped
toSkip <- readRDS('Data/Proxy/toSkip.rds')
proxy_MR1_toSearch <- proxy_MR1[!(proxy_MR1$SNP%in%toSkip),]

# Removing the SNPs in the outcome that are duplicated by filtering the duplicates by p-value
dups <- df3_GCST4023 %>% arrange(p_value)

dups_checked <- duplicated(dups$SNP) # Returns a list of Booleans for the SNP column in the dups table
df3_GCST4023 <- dups[!dups_checked,] # Keeps non-duplicate entries
saveRDS(df3_GCST4023, 'Data/Proxy/dupsRemoved.rds')


# Running the proxy search -----------------------------------------------------
outcomeToMerge <- NULL
proxySearch <- read.delim('Pancreas/Data_MR2/combined_query_snp_list_grch38.txt', row.names = NULL)

df2$SNP_Original <- df2$SNP_ID

token <- "2dd0e33b428c"
for(snp_p in proxy_MR2$SNP){
  print(snp_p)
  #proxy search
  ldp=proxySearch[proxySearch$query_snp==snp_p,]
  #saveRDS(snp_p,"Data/Proxy/Proxy_GWAS_Chiou_Pancreas_snp.rds")
  if(sum(proxySearch$query_snp==snp_p)==0){
    next
  }
  
  ldp <- ldp[ldp$R2 > 0.8,]
  #ldp$Coord <- gsub("chr","",ldp$Coord)
  ldp=ldp[ldp$Coord%in%df2$SNP,]
  if(nrow(ldp)==0){
    next
    }
  ldp2=ldp[1,]
  toMerge=df2[df2$SNP==ldp2$Coord,]
  toMerge$SNP=ldp2$Coord
  toMerge$SNP_Original=ldp2$RS_Number
  if(ldp2$Distance!=0){
    ref=colsplit(ldp2$Correlated_Alleles,',',c('a','b'))
    alt=colsplit(ref$b,'=',c('exp','out'))
    alt[alt==TRUE]='T'
    ref=colsplit(ref$a,'=',c('exp','out'))
    ref[ref==TRUE]='T'
    if(toMerge$EA==alt$out){
      toMerge$EA=alt$exp
      toMerge$Allele1=ref$exp
    }else{
      toMerge$EA=ref$exp
      toMerge$Allele1=alt$exp
    }
  }
  print("Binding...")
  outcomeToMerge=rbind(outcomeToMerge,toMerge)
  saveRDS(outcomeToMerge,"Pancreas/Data_MR2/Proxy/Proxy_GWAS_MichalekEUR_Pancreas_val.rds")
}

df2_Corrected <- rbind(df2, outcomeToMerge)
df2_Corrected <- df2_Corrected[df2_Corrected$Allele1%in%c('A','T','C','G') & df2_Corrected$EA%in%c('A','T','C','G'),]


saveRDS(df2_Corrected, 'Pancreas/Data_MR2/Corrected/proxiedGWAS_MichalekEUR.rds')