# Last Updated: 
# Author: Wesley Cai
# Purpose:

library(data.table)
library(ggplot2)
library(cBioPortalData)
library(maftools)
library(biomaRt)
library(reshape2)

# OG: Google_Drive/medschool/research/qin/tcga/210521_rna.R
mainwd <- "."
outputfolder <- "output/"
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

entrezid <- fread("../../dat/entrezid.dict.txt")
entrezid.goi <- entrezid[hgnc_symbol %in% c("KDM5B", "SETDB1", "CD8A", "CD3G", "CD8B", "IFNG", "TNF", "CXCL9", "CXCL10")]

cbio <- cBioPortal()
all.studies <- as.data.table(getStudies(cbio))

mol.profiles <- molecularProfiles(cbio, "skcm_tcga")
rna.zscore <- as.data.table(molecularData(cbio, "skcm_tcga_rna_seq_v2_mrna_median_Zscores", entrezGeneIds = entrezid.goi$entrezgene_id, sampleIds = allSamples(cbio, "skcm_tcga")$sampleId))
rna.rsem <- as.data.table(molecularData(cbio, "skcm_tcga_rna_seq_v2_mrna", entrezGeneIds = entrezid.goi$entrezgene_id, sampleIds = allSamples(cbio, "skcm_tcga")$sampleId))

colnames(rna.rsem) <- sub("^.*\\.", "", colnames(rna.rsem))
rna.rsem.cast <- as.data.table(dcast(rna.rsem, formula = patientId + sampleId ~ entrezGeneId))
rna.rsem.cast <- rna.rsem.cast[-grep("07", sampleId)]
# These patients are duplicated so choose met over primary
rna.rsem.cast <- rna.rsem.cast[-grep("TCGA-ER-A19T-01", sampleId)]
rna.rsem.cast <- rna.rsem.cast[-grep("TCGA-ER-A2NF-01", sampleId)]
which(duplicated(rna.rsem.cast$patientId))
rna.rsem.cast <- rna.rsem.cast[!duplicated(patientId)]
head(rna.rsem.cast)
rna.rsem.cast <- rna.rsem.cast[,.SD,.SDcols = c("sampleId", entrezid.goi$entrezgene_id)]
ncol(rna.rsem.cast) == nrow(entrezid.goi) + 1

rna.rsem.cast <- cbind(rna.rsem.cast[,.SD,.SDcols = c("sampleId")],
                       log2(rna.rsem.cast[,.SD,.SDcols = 2:ncol(rna.rsem.cast)] + 1))
colnames(rna.rsem.cast) <- c("sampleId", entrezid.goi$hgnc_symbol)

cor.res <- cor.test(rna.rsem.cast$KDM5B, rna.rsem.cast$SETDB1, alternative = "two.sided", method = "pearson")
p.value.lab <- signif(cor.res$p.value, 3)
rho.lab <- signif(cor.res$estimate, 3)
n.patients <- nrow(rna.rsem.cast)
mp <- ggplot(rna.rsem.cast, aes(KDM5B, SETDB1)) +
  geom_point() +
  geom_smooth(formula = "y ~ x", method = "lm", se = F, color = "red") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab(label = "SETDB1 log2(RSEM)") +
  xlab(label = "KDM5B log2(RSEM)") +
  annotate("text", x = 9, y = 12.3, label = paste0("P = ", p.value.lab), hjust = 0, size = 7) +
  annotate("text", x = 9, y = 12, label = paste0("r = ", rho.lab), hjust = 0, size = 7) +
  ggtitle(label = paste0("SKCM n = ", n.patients)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 18, color = "black"),
        panel.border = element_rect(fill=NA, colour = "black", size=2),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position = "none")
mp
ggsave(paste0("KDM5B_v_SETDB1.skcm.tcga.pdf"), mp, width = 5, height = 5.2)

cor.res <- cor.test(rna.rsem.cast$KDM5B, rna.rsem.cast$CD8A, alternative = "two.sided", method = "pearson")
p.value.lab <- signif(cor.res$p.value, 3)
rho.lab <- signif(cor.res$estimate, 3)
n.patients <- nrow(rna.rsem.cast)
mp <- ggplot(rna.rsem.cast, aes(KDM5B, CD8A)) +
  geom_point() +
  geom_smooth(formula = "y ~ x", method = "lm", se = F, color = "red") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab(label = "CD8A log2(RSEM)") +
  xlab(label = "KDM5B log2(RSEM)") +
  annotate("text", x = 8.9, y = 1.7, label = paste0("P = ", p.value.lab), hjust = 0, size = 7) +
  annotate("text", x = 8.9, y = 0.6, label = paste0("r = ", rho.lab), hjust = 0, size = 7) +
  ggtitle(label = paste0("SKCM n = ", n.patients)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 18, color = "black"),
        panel.border = element_rect(fill=NA, colour = "black", size=2),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position = "none")
mp
ggsave(paste0("KDM5B_v_CD8A.skcm.tcga.pdf"), mp, width = 5, height = 5)






rna.rsem.cast[KDM5B < quantile(KDM5B, c(0.5)), quant := "Low"]
rna.rsem.cast[KDM5B >= quantile(KDM5B, c(0.5)), quant := "High"]
rna.rsem.cast[, quant := factor(quant, c("Low", "High"))]
for(goi in c(entrezid.goi$hgnc_symbol)){
  t.test.res <- t.test(rna.rsem.cast[quant == "Low", get(goi)], rna.rsem.cast[quant == "High", get(goi)])
  p.value.lab <- formatC(signif(t.test.res$p.value, 3), format = "e", digits = 2)
  n.patients <- nrow(rna.rsem.cast)
  
  p.value.y.low <- quantile(rna.rsem.cast[quant == "Low", get(goi)], c(0.75))
  p.value.y.high <- quantile(rna.rsem.cast[quant == "High", get(goi)], c(0.75))
  p.value.y <- p.value.y.high
  if(p.value.y.low > p.value.y.high) p.value.y <- p.value.y.low
  
  mp <- ggplot(rna.rsem.cast, aes_string("quant", goi, fill = "quant")) +
    stat_boxplot(geom ='errorbar', width = 0.2, size = 0.8) + 
    geom_boxplot(notch=T, width = 0.5, size = 0.8, outlier.shape = NA) +
    ylab(label = paste0(goi, " log2(RSEM)")) +
    xlab(label = "KDM5B expression") +
    annotate("text", x = 1.1, y = p.value.y+1, label = paste0("P = ", p.value.lab), hjust = 0, size = 7, vjust = 0) +
    scale_fill_manual(values = c("dodgerblue2", "darkorange1")) +
    ggtitle(label = paste0("SKCM n = ", n.patients)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 18, hjust = 0.5),
          axis.title = element_text(size = 16, hjust = 0.5),
          axis.text = element_text(size = 18, color = "black"),
          panel.border = element_rect(fill=NA, colour = "black", size=2),
          axis.ticks = element_line(colour = "black", size = 1),
          legend.position = "none")
  mp
  ggsave(paste0("KDM5B_median_", goi,".skcm.tcga.pdf"), mp, width = 5, height = 5)
}
