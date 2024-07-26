# Header ---------------------------------------------------------------

# AUTHOR: John Eric Hamilton
# PURPOSE: Steiger Filtering 
# NOTES: 
#   *Inputs:   -  Harmonised Data files (rows for MR Result signif genes only )
#              
#
#   *Outputs:  - Steiger directionality TRUE/FALSE
#              - Steiger pval

# Load Harmonized Data Files ---------------------------------------------------
harm_total <- NULL
for (i in 1:22){
  harm_temp <- readRDS(sprintf('New/Spleen/Data_MR10/Harmonised_Data/Har_Chr%s.rds', i))
  harm_total <- rbind(harm_total, harm_temp)
}

# Filter Harmonized Data Files  
results <- readRDS('New/Spleen/Data_MR10/MR_Result/allResults_filteredWithNames.rds')
harm_total_Reduced <- harm_total[harm_total$id.exposure%in%results$id.exposure,]


# Add relevant columns for steiger filtering
harmdat11steiger <- harm_total_Reduced #creating a copy of the harmonized data since we need to modify it for the Steiger test of directionality
harmdat11steiger$units.outcome <- "log odds"

#harmdat11steiger$ncase.outcome <- 18942 # number of cases of T1D for Chiou
#harmdat11steiger$ncontrol.outcome <- 501638 # number of healthy controls for Chiou
#harmdat11steiger$prevalence.outcome <- (18942/(18942+501638)) # cases/total for Chiou

#harmdat11steiger$ncase.outcome <- 3428 # number of cases of T1D for Michalek EUR
#harmdat11steiger$ncontrol.outcome <- 3428 # number of healthy controls for Michalek EUR 
#harmdat11steiger$prevalence.outcome <- (3428/(3428+3428)) #cases/total for Michalek EUR

harmdat11steiger$ncase.outcome <- 3990 # number of cases of T1D for Michalek Transethnic
harmdat11steiger$ncontrol.outcome <- 4065 # number of healthy controls for Michalek Transethnic
harmdat11steiger$prevalence.outcome <- (3990/(3990+4065)) # cases/total for Michalek Transethnic


harmdat11steiger$effective_n.exposure <- 241 #Samples in GTEx (P: 328, WB: 755, Spleen: 241)
#harmdat11steiger$effective_n.outcome <- 4/(1/501638 + 1/18942) #effect size calculation using ncase and ncontrol from outcome Chiou

#harmdat11steiger$effective_n.outcome <- 4/(1/3428 + 1/3428) #effect size calculation using ncase and ncontrol from outcome Michalek EUR

harmdat11steiger$effective_n.outcome <- 4/(1/3990 + 1/4065) #effect size calculation using ncase and ncontrol from outcome Michalek Transethnic


  
# Adds the values of R^2 calculated 
harmdat11steiger$rsq.exposure <- (get_r_from_bsen(harm_total_Reduced$beta.exposure,harm_total_Reduced$se.exposure,241))^2 #(P: 328, WB: 755, Spleen: 241)
#harmdat11steiger$rsq.outcome <- (get_r_from_bsen(harm_total_Reduced$beta.outcome,harm_total_Reduced$se.outcome,(18942+501638)))^2

harmdat11steiger$rsq.outcome <- get_r_from_lor(
  lor = (harmdat11steiger$beta.outcome),
  af = harmdat11steiger$eaf.outcome,
  ncase = harmdat11steiger$ncase.outcome,
  ncontrol = harmdat11steiger$ncontrol.outcome,
  prevalence = harmdat11steiger$prevalence.outcome,
  model = "logit",
  correction = FALSE)^2


mr11_steig_filt <- steiger_filtering(harmdat11steiger)
View(mr11_steig_filt)

if (sum(mr11_steig_filt$steiger_dir) == nrow(harm_total_Reduced)){
  print("Steiger Passed")
  saveRDS(mr11_steig_filt, 'New/Steiger/steiger_MR10.rds')
}



# The directionality test is run as support to the results obtained from steiger_filtering.
# Gives a lot of False values for directionality test
mr11_direct_res <- directionality_test(harmdat11steiger)
