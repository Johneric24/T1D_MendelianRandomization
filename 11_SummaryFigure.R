# Header ---------------------------------------------------------------------

# AUTHOR: John Eric Hamilton
# PURPOSE: Making Barplots for final results figure
# NOTES: 
#   *Inputs:   - MR results  
#
#   *Outputs:  - Barplots that show OR and gene names 

# Load library -----------------------------------------------------------------
library(ggplot2)
library(ggpubr)

# Load Data Pancreas -----------------------------------------------------------
pancreas <- NULL
for (i in c(1,2,4)){
  # Get Data
  temp <- readRDS(sprintf('New/Pancreas/Data_MR%s/MR_Result/allResults_filteredWithNames.rds', i))
  if (i == 1){
    temp$GWAS <- "Chiou"
  } else if (i==2){
    temp$GWAS <- "Michalek EUR"
  } else if (i ==4){
    temp$GWAS <- "Michalek Transethnic"
  }
  
  temp <- temp[!duplicated(temp$geneName),]
  
  
  pancreas <- rbind(pancreas, temp)
}

# Only keep those with successful coloc/susie
passing_pancreas <- c(
  "NAA25",
  "PRKD2",
  "SH2B3",
  "METTL21A",
  "PTEN",
  "LIPJ",
  "DNAJC27-AS1",
  "NOTCH2",
  "SULT1A2",
  "CDC37P1",
  "CBFA2T3",
  "ST7L",
  "CLECL1P",
  "CTSH",
  "RPL7AP16",
  "NTN5",
  "SBK1",
  "LINC00240",
  "BTN2A2"
)
toPlot_pancreas <- pancreas[pancreas$geneName%in%passing_pancreas,]
toPlot_pancreas$Tissue <- "Pancreas"

toPlot_pancreas$geneName=as.character(toPlot_pancreas$geneName)
toPlot_pancreas$GWAS=as.character(toPlot_pancreas$GWAS)

# Order Beta
toPlot_pancreas[toPlot_pancreas$geneName == "LIPJ" & toPlot_pancreas$GWAS == "Chiou",]$geneName <- "LIPJ (C)"
toPlot_pancreas[toPlot_pancreas$geneName == "LIPJ" & toPlot_pancreas$GWAS == "Michalek EUR",]$geneName <- "LIPJ (M-Eur)"
toPlot_pancreas[toPlot_pancreas$geneName == "LIPJ" & toPlot_pancreas$GWAS == "Michalek Transethnic",]$geneName <- "LIPJ (M-T)"

toPlot_pancreas[toPlot_pancreas$geneName == "PTEN" & toPlot_pancreas$GWAS == "Chiou",]$geneName <- "PTEN (C)"
toPlot_pancreas[toPlot_pancreas$geneName == "PTEN" & toPlot_pancreas$GWAS != "Chiou",]$geneName <- "PTEN (M-T)"

toPlot_pancreas[toPlot_pancreas$geneName == "LINC00240" & toPlot_pancreas$GWAS == "Chiou",]$geneName <- "LINC00240 (C)"
toPlot_pancreas[toPlot_pancreas$geneName == "LINC00240" & toPlot_pancreas$GWAS == "Michalek EUR",]$geneName <- "LINC00240 (M-Eur)"
toPlot_pancreas[toPlot_pancreas$geneName == "LINC00240" & toPlot_pancreas$GWAS == "Michalek Transethnic",]$geneName <- "LINC00240 (M-T)"

toPlot_pancreas <- toPlot_pancreas[order(toPlot_pancreas$b),]

toPlot_pancreas$geneName <- factor(toPlot_pancreas$geneName,levels=(unique(toPlot_pancreas$geneName)))



# Identifying Repeats
#toPlot_pancreas$Tissue <- "Pancreas"
#toPlot_pancreas[toPlot_pancreas$geneName == "DNAJC27-AS1",]$Tissue <- "Pancreas + WB"
#toPlot_pancreas[toPlot_pancreas$geneName%in%c("CDC37P1", "CLECL1P", "CTSH","METTL21A"),]$Tissue <- "Pancreas + Spleen"


png(filename="New/finalGenesPancreas.png", width = 600, height = 480)
toPlot_pancreas$xmin=seq(0.55,24,1)
toPlot_pancreas$xmax=seq(1.45,24.5,1)

a <- ggplot(data=toPlot_pancreas, aes(y=exp(b),x=geneName, fill=GWAS))+
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 1, ymax = exp(b)),size=2)+
  scale_y_continuous(trans='log10')+
  #geom_bar(stat="identity", position=position_dodge(),linewidth=0.75)+
  scale_fill_manual(values=c("#56B4E9", "#CC79A7", "#009E73"), na.value = NA)+
  geom_errorbar(aes(ymin=exp(b-1.96*se), ymax=exp(b+1.96*se)), width=.4,linewidth=1, position=position_dodge(.9))+ ggtitle('Significant Genes Pancreas') + xlab('Gene Names') + ylab('Odd Ratio (OR)') + theme_bw() + 
  #facet_wrap(~GWAS)+
  theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust=1, size=13), plot.margin = margin(1,1,1,1, "cm"))
  
dev.off()


# Load Data Whole_Blood --------------------------------------------------------
WB <- NULL
for (i in c(5, 6)){
  # Get Data
  temp <- readRDS(sprintf('New/WB/Data_MR%s/MR_Result/allResults_filteredWithNames.rds', i))
  if (i == 5){
    temp$GWAS <- "Chiou"
  } else if (i == 6){
    temp$GWAS <- "Michalek EUR"
  }
  
  temp <- temp[!duplicated(temp$geneName),]
  
  
  WB <- rbind(WB, temp)
}

passing_WB <- c(
  "IL27",
  "SLC22A4",
  "CAMK4",
  "ISYNA1",
  "SEPTIN2",
  "PIK3R3",
  "LINC02390",
  "MAPK8IP1P1",
  "ACAP1",
  "SKAP2",
  "HOXA4",
  "GSDMB",
  "HOXA1",
  "HOXA-AS2",
  "ENSG00000236935_lncRNA",
  "CENPO",
  "DNAJC27",
  "ANKRD55",
  "RRP7BP",
  "IKZF3",
  "MTMR3",
  "CCDC88B",
  "CDC42SE2",
  "ATP5MG",
  "H2BC15"
)

toPlot_WB <- WB[WB$geneName%in%passing_WB,]
toPlot_WB$Tissue <- "Whole Blood"

# Order Beta
toPlot_WB <- toPlot_WB[order(toPlot_WB$b),]
toPlot_WB$geneName <- factor(toPlot_WB$geneName,levels=(unique(toPlot_WB$geneName)))

png(filename="New/finalGenesWB.png", width = 600, height = 480)

toPlot_WB$xmin=seq(0.55,24,1)
toPlot_WB$xmax=seq(1.45,24.5,1)

b <- ggplot(data=toPlot_WB, aes(y=exp(b),x=geneName, fill=GWAS))+
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 1, ymax = exp(b)),size=2)+
  scale_y_continuous(trans='log10')+
  #geom_bar(stat="identity", position=position_dodge(),linewidth=0.75)+
  scale_fill_manual(values=c("#56B4E9", "#CC79A7", "#009E73"), na.value = NA)+
  geom_errorbar(aes(ymin=exp(b-1.96*se), ymax=exp(b+1.96*se)), width=.4,linewidth=1, position=position_dodge(.9))+ ggtitle('Significant Genes WB') + xlab('Gene Names') + ylab('Odd Ratio (OR)') + theme_bw() + 
  #facet_wrap(~GWAS)+
  theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust=1, size=13), plot.margin = margin(1,1,1,1, "cm")) 
  #scale_y_continuous(limits = c(0, 5))

dev.off()

# Load Data Spleen -------------------------------------------------------------
Spleen <- NULL
for (i in c(8, 9, 10)){
  # Get Data
  temp <- readRDS(sprintf('New/Spleen/Data_MR%s/MR_Result/allResults_filteredWithNames.rds', i))
  if (i == 8){
    temp$GWAS <- "Chiou"
  } else if (i == 9){
    temp$GWAS <- "Michalek EUR"
  } else if (i == 10){
    temp$GWAS <- "Michalek Transethnic"
  }
  
  temp <- temp[!duplicated(temp$geneName),]
  
  
  Spleen <- rbind(Spleen, temp)
}

passing_spleen <- c("ACAP1",
                    "METTL21A",
                    "TYK2",
                    "SUOX",
                    "LINC01806",
                    "LINC02362",
                    "LINC02362",
                    "ENSG00000236935_lncRNA",
                    "ENSG00000262539_processed_pseudogene",
                    "CDC37P1",
                    "SYNGR1",
                    "ENSG00000280216_TEC",
                    "ANKRD55",
                    "ZPBP2",
                    "CLECL1P",
                    "CTSH",
                    "DOK6",
                    "RASGRP1",
                    "MTMR3",
                    "PFKL",
                    "CENPW",
                    "SLC17A3",
                    "ZSCAN31"
)

toPlot_spleen <- Spleen[Spleen$geneName%in%passing_spleen,]
toPlot_spleen$Tissue <- "Spleen"

toPlot_spleen$geneName=as.character(toPlot_spleen$geneName)
toPlot_spleen$GWAS=as.character(toPlot_spleen$GWAS)

# Order Beta
toPlot_spleen[toPlot_spleen$geneName == "METTL21A" & toPlot_spleen$GWAS == "Michalek EUR",]$geneName <- "METTL21A (M-Eur)"

toPlot_spleen[toPlot_spleen$geneName == "METTL21A" & toPlot_spleen$GWAS == "Michalek Transethnic",]$geneName <- "METTL21A (M-T)"


toPlot_spleen[toPlot_spleen$geneName == "SUOX" & toPlot_spleen$GWAS == "Chiou",]$geneName <- "SUOX (C)"

toPlot_spleen[toPlot_spleen$geneName == "SUOX" & toPlot_spleen$GWAS == "Michalek EUR",]$geneName <- "SUOX (M-Eur)"

toPlot_spleen[toPlot_spleen$geneName == "SUOX" & toPlot_spleen$GWAS == "Michalek Transethnic",]$geneName <- "SUOX (M-T)"



toPlot_spleen[toPlot_spleen$geneName == "LINC01806" & toPlot_spleen$GWAS == "Michalek EUR",]$geneName <- "LINC01806 (M-Eur)"

toPlot_spleen[toPlot_spleen$geneName == "LINC01806" & toPlot_spleen$GWAS == "Michalek Transethnic",]$geneName <- "LINC01806 (M-T)"

toPlot_spleen[toPlot_spleen$geneName == "ZSCAN31" & toPlot_spleen$GWAS == "Michalek EUR",]$geneName <- "ZSCAN31 (M-Eur)"

toPlot_spleen[toPlot_spleen$geneName == "ZSCAN31" & toPlot_spleen$GWAS == "Michalek Transethnic",]$geneName <- "ZSCAN31 (M-T)"




toPlot_spleen <- toPlot_spleen[order(toPlot_spleen$b),]
toPlot_spleen$geneName <- factor(toPlot_spleen$geneName,levels=(unique(toPlot_spleen$geneName)))

png(filename="New/finalGenesSpleen.png", width = 600, height = 480)

toPlot_spleen$xmin=seq(0.55,24,1)
toPlot_spleen$xmax=seq(1.45,24.5,1)

c <- ggplot(data=toPlot_spleen, aes(y=exp(b),x=geneName, fill=GWAS))+
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 1, ymax = exp(b)),size=2)+
  scale_y_continuous(trans='log10')+
  #geom_bar(stat="identity", position=position_dodge(),linewidth=0.75)+
  scale_fill_manual(values=c("#56B4E9", "#CC79A7", "#009E73"), na.value = NA)+
  geom_errorbar(aes(ymin=exp(b-1.96*se), ymax=exp(b+1.96*se)), width=.4,linewidth=1, position=position_dodge(.9))+ ggtitle('Significant Genes Pancreas') + xlab('Gene Names') + ylab('Odd Ratio (OR)') + theme_bw() + 
  #facet_wrap(~GWAS)+
  theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust=1, size=13), plot.margin = margin(1,1,1,1, "cm"))

dev.off()


saveRDS(toPlot_pancreas, 'New/finalGenesPancreas.rds')
saveRDS(toPlot_WB, 'New/finalGenesWB.rds')
saveRDS(toPlot_spleen, 'New/finalGenesSpleen.rds')

png(filename="New/allGenesFigure.png", width = 1300, height = 480)
ggarrange(a, b, c, ncol = 3, nrow = 1, common.legend=TRUE,legend='right',widths = c(1.5,1.5,1.5),align = "h")
# le widths ajuste la largeur de chaque figure. le align est celle qui permet d'aligner en bas les graphiques
dev.off()


# Load Data Replicated Whole_Blood -----------------------------------------------------------------------------
WB_replicated <- NULL
for (i in c(1, 2 ,4)){
  # Get Data
  temp <- readRDS(sprintf('eQTLgen/MAF/Data_MR%s/MR_Result/allResultSignifWithNames.rds', i))
  if (i == 1){
    temp$GWAS <- "Chiou"
  } else if (i==2){
    temp$GWAS <- "Michalek EUR"
  } else if (i ==4){
    temp$GWAS <- "Michalek Transethnic"
  }
  
  WB_replicated <- rbind(WB_replicated, temp)
}

passing_WB_replicated <- c("PIK3R3",
                    "SLC17A3",
                    "CENPW",
                    "ZSCAN31",
                    "SKAP2",
                    "CCDC88B",
                    "SH2B3",
                    "SUOX",
                    "CTSH",
                    "SULT1A2",
                    "ACAP1",
                    "PRKD2",
                    "NTN5",
                    "SYNGR1"
)

WB_replicated <- WB_replicated[WB_replicated$geneName%in%passing_WB_replicated,]

signif <- 0.05/50
WB_replicated_signif <- WB_replicated[WB_replicated$pval < signif,]

WB_replicated_signif$Tissue <- "Whole Blood"
WB_replicated_signif$geneName=as.character(WB_replicated_signif$geneName)
WB_replicated_signif$GWAS=as.character(WB_replicated_signif$GWAS)

WB_replicated_signif[WB_replicated_signif$geneName == "ZSCAN31" & WB_replicated_signif$GWAS == "Michalek EUR",]$geneName <- "ZSCAN31 (M-Eur)"

WB_replicated_signif[WB_replicated_signif$geneName == "ZSCAN31" & WB_replicated_signif$GWAS == "Michalek Transethnic",]$geneName <- "ZSCAN31 (M-T)"

WB_replicated_signif[WB_replicated_signif$geneName == "ZSCAN31" & WB_replicated_signif$GWAS == "Chiou",]$geneName <- "ZSCAN31 (c)"

# Check for duplicates after only for this one
WB_replicated_signif <- WB_replicated_signif[!(duplicated(WB_replicated_signif$geneName)),]

# Order Beta
WB_replicated_signif <- WB_replicated_signif[order(WB_replicated_signif$b),]
WB_replicated_signif$geneName <- factor(WB_replicated_signif$geneName,levels=(unique(WB_replicated_signif$geneName)))

png(filename="eQTLgen/finalGenesExpressionOnly_updated.png", width = 480, height = 480)

WB_replicated_signif$xmin=seq(0.55,16,1)
WB_replicated_signif$xmax=seq(1.45,16.5,1)

ggplot(data=WB_replicated_signif, aes(y=exp(b),x=geneName, fill=GWAS))+
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 1, ymax = exp(b)),size=1.5)+scale_y_continuous(trans='log10')+
  #geom_bar(stat="identity", position=position_dodge(),linewidth=0.75)+
  scale_fill_manual(values=c("#56B4E9", "#CC79A7", "#009E73"), na.value = NA)+
  geom_errorbar(aes(ymin=exp(b-1.96*se), ymax=exp(b+1.96*se)), width=.4,linewidth=1, position=position_dodge(.9))+ 
  theme_bw() + ggtitle('Significant Genes Replicated Whole Blood') + xlab('Gene Names') + ylab('Odd Ratio (OR)') +
  #facet_wrap(~GWAS)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=12))

dev.off()



# Common Genes with sQTL -------------------------------------------------------
toPlot_total <- rbind(toPlot_pancreas, toPlot_WB, toPlot_spleen)

common <- c(
  "ST7L",
  "SEPTIN2",
  "PRKD2",
  "SUOX",
  "CTSH",
  "NOTCH2",
  "ACAP1",
  "TYK2",
  "CCDC88B",
  "ACAP1"
)

common_genes <- toPlot_total[toPlot_total$geneName%in%common,]

# Replicated in different tissues-------------------------------------------------------------------
WB_spleen <- toPlot_total[toPlot_total$geneName%in%c("ANKRD55", "MTMR3", "ACAP1"),]


P_spleen <- toPlot_total[toPlot_total$geneName%in%c("METTL21A", "METTL21A (M-T)", "CLECL1P", "CDC37P1"),]
P_spleen[P_spleen$geneName == "METTL21A (M-T)" ,]$geneName <- "METTL21A"


png(filename="New/WB_Spleen.png", width = 480, height = 480)
WB_spleen$xmin=seq(0.55,2,1)
WB_spleen$xmax=seq(1.45,2.5,1)
ggplot(data=WB_spleen, aes(y=exp(b),x=geneName, fill=GWAS))+
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 1, ymax = exp(b)),size=1.5)+scale_y_continuous(trans='log10')+
  #geom_bar(stat="identity", position=position_dodge(),linewidth=0.75)+
  scale_fill_manual(values=c("#56B4E9", "#CC79A7", "#009E73"), na.value = NA)+
  geom_errorbar(aes(ymin=exp(b-1.96*se), ymax=exp(b+1.96*se)), width=.4,linewidth=1, position=position_dodge(.9))+ 
  theme_bw() + ggtitle('Significant Genes Shared Pancreas Spleen') + xlab('Gene Names') + ylab('Odd Ratio (OR)') +
  facet_wrap(~Tissue)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=12))

dev.off()

# Other type of graph --------------------------------------------------------
all <- rbind(WB_spleen, P_spleen)
all$xmin = 1
all$xmax = 1


png(filename="New/allRepeats.png", width = 480, height = 480)

all[all$Tissue == "Pancreas",]$xmin=0.55
all[all$Tissue == "Pancreas",]$xmax=1.55

all[all$Tissue == "Whole Blood",]$xmin=2.55
all[all$Tissue == "Whole Blood",]$xmax=3.55

all[all$Tissue == "Spleen",]$xmin=1.55
all[all$Tissue == "Spleen",]$xmax=2.55


ggplot(data=all, aes(y=exp(b),x=Tissue, fill=Tissue))+
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 1, ymax = exp(b)),size=1.5)+scale_y_continuous(trans='log10')+
  #geom_bar(stat="identity", position=position_dodge(),linewidth=0.75)+
  scale_fill_manual(values=c("#56B4E9", "#CC79A7", "#009E73"), na.value = NA)+
  geom_errorbar(aes(ymin=exp(b-1.96*se), ymax=exp(b+1.96*se)), width=.4,linewidth=1, position=position_dodge(.9))+ 
  theme_bw() + ggtitle('Significant Genes Replicated') + xlab('Gene Names') + ylab('Odd Ratio (OR)') +
  facet_wrap(~geneName)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=12))

dev.off()

# Splicing ----------------------------------------------------------------
toPlot_total <- readRDS('New/toPlot_total.rds')
WB_spleen <- toPlot_total[toPlot_total$geneName%in%c("SIRPG"),]
P_WB <- toPlot_total[toPlot_total$geneName%in%c("SULT1A1"),]
P_WB_spleen<- toPlot_total[toPlot_total$geneName%in%c("ST7L",
                                                      "SEPTIN2",
                                                      "PRKD2",
                                                      "SUOX",
                                                      "CTSH",
                                                      "GSDMB"),]

all <- rbind(WB_spleen, P_WB, P_WB_spleen)
all$xmin <- 1
all$xmax <- 1

png(filename="New/allSplicingRepeats.png", width = 480, height = 480)

all[all$Tissue == "Pancreas",]$xmin <- 0.55
all[all$Tissue == "Pancreas",]$xmax <- 1.55

all[all$Tissue == "Whole Blood",]$xmin=2.55
all[all$Tissue == "Whole Blood",]$xmax=3.55

all[all$Tissue == "Spleen",]$xmin=1.55
all[all$Tissue == "Spleen",]$xmax=2.55


ggplot(data=all, aes(y=exp(b),x=Tissue, fill=Tissue))+
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 1, ymax = exp(b)),size=1.5)+scale_y_continuous(trans='log10')+
  #geom_bar(stat="identity", position=position_dodge(),linewidth=0.75)+
  scale_fill_manual(values=c("#56B4E9", "#CC79A7", "#009E73"), na.value = NA)+
  geom_errorbar(aes(ymin=exp(b-1.96*se), ymax=exp(b+1.96*se)), width=.4,linewidth=1, position=position_dodge(.9))+ 
  theme_bw() + ggtitle('Significant Genes Replicated') + xlab('Gene Names') + ylab('Odd Ratio (OR)') +
  facet_wrap(~geneName)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=12))

dev.off()









