# Header ----------------------------------------------------------------------

# AUTHOR: John Eric Hamilton
# PURPOSE: Proxy Search for SNPs that are in the exposure and are missing the outcome
# NOTES: 
#   *Inputs:   - List of missing SNPs
#
#   *Outputs:  - Proxy Data 
# Load library -----------------------------------------------------------------
library("LDlinkR")

# Running the proxy search for single SNP --------------------------------------
snp_p <- "chr11:1341591"

token <- "2dd0e33b428c"

#proxy search
ldp=tryCatch(LDproxy(snp_p, token = token, genome_build = "grch38", pop = "ALL"))


# Running the proxy search for SNP batch ---------------------------------------
proxy_MR2 <- readRDS('missingFromOutcome.rds')

token <- "###"

#proxy search
LDproxy_batch(proxy_MR2$SNP, token = token, genome_build = "grch38", append = TRUE)