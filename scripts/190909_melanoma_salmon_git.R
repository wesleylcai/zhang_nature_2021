# Last Updated: 
# Author: Wesley Cai
# Purpose:

library(data.table)
library(ggplot2)
library(DESeq2)
library(tximport)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# OG: Google_Drive/medschool/research/qin/melanoma/GSE78220_hugo/melanoma_pembro/salmon/190909_melanoma_salmon_edited.R

mainwd <- "."
outputfolder <- "output/"
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

clinical <- fread("../../dat/hugo.clinical.txt")
clinical[irRECIST == "Progressive_Disease",irR := gsub("Progressive_Disease", "PD", irRECIST)]
clinical[irRECIST == "Complete_Response",irR := gsub("Complete_Response", "CR", irRECIST)]
clinical[irRECIST == "Partial_Response",irR := gsub("Partial_Response", "PR", irRECIST)]
clinical.df <- as.data.frame(clinical)
row.names(clinical.df) <- clinical.df$Patient_ID

load(file = "../../dat/hugo.vst.raw.dt.RData")
vst.raw.dt.t <- transpose(vst.raw.dt, keep.names = "Patient_ID", make.names = 1)

#### Correlations ####
# # Find ERV correlating with KDM5B expression
# # PDCD1 = PD-1, CD274 = PD-L1
# cor.genes <- c("CD8A", "IFNG", "CCR5", "CXCL9", "CXCL10", "CXCL11", "IDO1", "PRF1", "GZMA")
# for(cor.gene in cor.genes){
#   dir.create(paste0(cor.gene, "/", cor.gene, "_v_erv"), showWarnings = FALSE, recursive = TRUE)
#   plotResponse(cor.gene, outfolder = cor.gene)
#   which(colnames(vst.raw.df.t) == cor.gene)
#   erv.cor.res <- data.table(erv = grep("ERVmap", vst.raw.dt$hgnc, value = TRUE))
#   # Filter for those ERVs concordant with KDM5B
#   erv.cor.res <- erv.cor.res[erv %in% kdm5b.erv.cor.res[kdm5b.cor.p < 0.05 & kdm5b.cor.r < 0,get("erv")]]
#   rown <- 1
#   for(rown in 1:nrow(erv.cor.res)){
#     erv.symbol <- erv.cor.res[rown,get("erv")]
#     erv.exp <- vst.raw.df.t[,.SD,.SDcols = c(cor.gene, erv.symbol)]
#     colnames(erv.exp) <- c("cor.gene", "erv")
#     cor.res <- cor.test(erv.exp$cor.gene, erv.exp$erv)
#     erv.cor.res[rown, cor.p := signif(cor.res$p.value,3)]
#     erv.cor.res[rown, cor.r := signif(cor.res$estimate,3)]
#     if(!is.na(cor.res$p.value)){
#       if(cor.res$p.value < 0.05 & abs(cor.res$estimate) > 0.2){
#         if(plotResponse(erv.symbol, pval.thr = 0.05, outfolder = cor.gene) < 0.05){
#           mp <- ggplot(erv.exp, aes_string("cor.gene", "erv")) +
#             geom_point() +
#             geom_smooth(method = "lm") +
#             ylab(label = erv.symbol) +
#             xlab(label = cor.gene) +
#             theme_bw()
#           ggsave(paste0(paste0(cor.gene, "/", cor.gene, "_v_erv/"), erv.symbol, "_v_", cor.gene, ".png"), mp, width = 4, height = 4)
#         }
#       }
#     }
#   }
#   save(erv.cor.res, file = paste0(cor.gene, "/", cor.gene, ".erv.cor.res.RData"))
# }
#### Correlations ####
load("../../dat/kdm5b.erv.cor.res.RData")

#### HEATMAPS ####
ervmap <- vst.raw.dt[grep("ERVmap|KDM5B", hgnc),]
ervmap[, rowvars := rowVars(as.matrix(ervmap[,.SD,.SDcols = grep("Pt", colnames(ervmap))]))]
ervmap <- ervmap[order(rowvars, decreasing = TRUE),]

# Comp 3
kdm5b.sigcor.erv <- kdm5b.erv.cor.res[kdm5b.cor.p < 0.05, get("erv")]
ervmap.kdm5bcor <- ervmap[hgnc %in% c(kdm5b.sigcor.erv, "KDM5B")]
ervmap.kdm5bcor.mat <- as.matrix(ervmap.kdm5bcor[,.SD,.SDcols = grep("Pt", colnames(ervmap.kdm5bcor))])
row.names(ervmap.kdm5bcor.mat) <- ervmap.kdm5bcor$hgnc
ervmap.kdm5bcor.mat.scaled <- t(scale(t(ervmap.kdm5bcor.mat)))
erv.sorted <- ervmap.kdm5bcor.mat.scaled[,names(sort(ervmap.kdm5bcor.mat.scaled["KDM5B",]))]
response.annote <- HeatmapAnnotation(response = clinical.df[colnames(erv.sorted),"irR"], col = list(response = c("PD" = "red", "PR" = "green", "CR" = "blue")))
hm <- Heatmap(erv.sorted, show_row_names = TRUE, cluster_columns = FALSE,
              col = colorRamp2(seq(-4,4,length.out = 11), rev(brewer.pal(11, "RdBu"))),
              top_annotation = response.annote, bottom_annotation = response.annote)
png("kdm5b_erv_cor_heatmap_kdm5b_sorted.png", height = 15, width = 8, units = "in", res = 300)
print(hm)
dev.off()
#### HEATMAPS ####

# Manual plot
goi <- "KDM5B"
dat <- merge(clinical[,.(Patient_ID,irRECIST)], vst.raw.dt.t[,.SD,.SDcols = c("Patient_ID", goi)], by = "Patient_ID")
pval <- tryCatch(signif(t.test(dat[irRECIST == "Complete_Response",get(goi)], dat[irRECIST == "Progressive_Disease",get(goi)])$p.value, 3),
                 error = function(e) return(1))
pos <- position_jitter(width = 0.2, seed = 2)
mp <- ggplot(dat, aes_string("irRECIST", goi, fill = "irRECIST")) +
  stat_boxplot(geom ='errorbar', width = 0.2, size = 0.5) + 
  geom_boxplot(width = 0.5, show.legend = FALSE) +
  geom_point(alpha = 1, position = pos, color = "black", show.legend = F) +
  annotate("text", label = paste0("CR vs PD p = ", pval), x = 2, y = max(dat[,get(goi)])+0.2) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 18),
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=0.75))
ggsave("KDM5B_figure1b.v4.pdf", width = 3, height = 2.2)

# Plot KDM5B expression in response groups
# https://stackoverflow.com/questions/55933524
library(showtext)
font_add("Arial", "/Applications/Microsoft Word.app/Contents/Resources/DFonts/arial.ttf")
showtext_auto()

scaleFUN <- function(x) sprintf("%.1f", x)
plotResponse <- function(goi, pval.thr = 1, outfolder = "."){
  dir.create(file.path(outfolder), recursive = TRUE, showWarnings = FALSE)
  dat <- merge(clinical[,.(Patient_ID,irRECIST)], vst.raw.dt.t[,.SD,.SDcols = c("Patient_ID", goi)], by = "Patient_ID")
  pval <- tryCatch(signif(t.test(dat[irRECIST == "Complete_Response",get(goi)], dat[irRECIST == "Progressive_Disease",get(goi)])$p.value, 3),
                   error = function(e) return(1))
  if(pval < pval.thr){
    labelscale <- max(dat[,get(goi)])+(max(dat[,get(goi)])-min(dat[,get(goi)]))/10
    pos <- position_jitter(width = 0.2, seed = 2)
    mp <- ggplot(dat, aes_string("irRECIST", goi, fill = "irRECIST")) +
      stat_boxplot(geom ='errorbar', width = 0.2, size = 0.5) + 
      geom_boxplot(width = 0.5, show.legend = FALSE) +
      geom_point(alpha = 1, position = pos, color = "black", show.legend = F) +
      annotate("text", label = paste0("CR vs PD p = ", pval), x = 2, y = labelscale) +
      theme_bw() +
      scale_y_continuous(labels = scaleFUN) +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 12, color = "black"),
            axis.title.y = element_text(size = 18),
            panel.grid = element_blank(),
            panel.border = element_rect(fill=NA, colour = "black", size=0.75),
            text=element_text(family="Arial"))
    ggsave(paste0(outfolder, "/", goi, "_response.pdf"), mp, width = 2.5, height = 3.5)
    mp <- ggplot(dat, aes_string("irRECIST", goi, fill = "irRECIST")) +
      geom_boxplot() +
      annotate("text", label = paste0("CR vs PD p = ", pval), x = 2, y = max(dat[,get(goi)])+0.2) +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank())
    ggsave(paste0(outfolder, "/FOR_LEGEND_response.pdf"), mp, width = 4, height = 3.5)
  }
  return(pval)
}
gois <- c("ERVmap_2637", "KDM5A", "KDM5C", "KDM5D", "KDM1A", "DNMT1", "DNMT3B", "EZH2", "SETDB1")
for(goi in gois){
  plotResponse(goi, outfolder = ".")
}

#### Correlation plots ####
gene1 <- "KDM5B"
gene2 <- "ERVmap_2637"
vst.raw.dt.t.subset <- vst.raw.dt.t[,.SD,.SDcols = c(gene1,gene2)]
cor.res <- cor.test(vst.raw.dt.t.subset[,get(gene1)], vst.raw.dt.t.subset[,get(gene2)])
mp <- ggplot(vst.raw.dt.t.subset, aes_string(gene1, gene2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  annotate("text", 12.2, 10.4, hjust = 0, label = paste0("r2 = ", signif(cor.res$estimate,3), "\np = ", signif(cor.res$p.value,3))) +
  theme_bw() +
  scale_y_continuous(labels = scaleFUN) +
  scale_x_continuous(labels = scaleFUN) +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)),
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=0.75),
        text=element_text(family="Arial"))
ggsave(paste0(gene1, "_v_", gene2, ".cor.pdf"), mp, width = 3.5, height = 3.5)
#### Correlation plots ####
