# Last Updated: 
# Author: Wesley Cai
# Purpose:

library(data.table)
library(ggplot2)
library(DESeq2)
library(reshape2)
library(openxlsx)
library(ggrepel)
library(scales)

# OG: Google_Drive/medschool/research/qin/atac/analysis/210513_atac_diffbind_sizefactors_update.R
mainwd <- "."
outputfolder <- "output/"
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

#### Load files ####
# File can be found at GSE161065
rep.atac <- fread("../../dat/repmask.atac.bed")
res.dt <- fread("../../dat/atac_results_library.sizefactor.txt")
#### Load files ####

#### For Motif analysis ####

# More accessible KO vs WT -IFNG
# Less accessible KO vs WT -IFNG
# More accessible KO vs WT +IFNG
# Less accessible KO vs WT +IFNG
# Common IFNG induced
# Common IFNG inhibited
# KO specific IFNG induced
# WT specific IFNG induced
motif.input.list <- list()

res.dt[, c("chr", "start", "end") := as.data.table(do.call(rbind, strsplit(peakid, "-")))]
motif.input.list[["kovwt.noifng.up"]] <- res.dt[ko4H3_v_YR.l2fc > 0 & ko4H3_v_YR.padj < 0.05, .(chr, start, end, peakid)]
motif.input.list[["kovwt.noifng.down"]] <- res.dt[ko4H3_v_YR.l2fc < 0 & ko4H3_v_YR.padj < 0.05, .(chr, start, end, peakid)]
motif.input.list[["kovwt.ifng.up"]] <- res.dt[ko4H3.IFNG_v_YR.IFNG.l2fc > 0 & ko4H3.IFNG_v_YR.IFNG.padj < 0.05, .(chr, start, end, peakid)]
motif.input.list[["kovwt.ifng.down"]] <- res.dt[ko4H3.IFNG_v_YR.IFNG.l2fc < 0 & ko4H3.IFNG_v_YR.IFNG.padj < 0.05, .(chr, start, end, peakid)]

# More accessible under all repetitive elements -IFNG
# More accessible under all repetitive elements +IFNG
# More accessible under all repetitive elements -IFNG overlap with KDM5B peaks
# More accessible under all repetitive elements +IFNG overlap with KDM5B peaks

overlaps <- fread("../../dat/allYRandall4H3-ATAC.repmask.YR-KDM5B.closest.bed")
overlaps <- merge(overlaps[,.(peakid = V4, repmask.name = V8, rep.dist = V9, kdm5b.dist = V14)], res.dt, by = "peakid", all.x = T)
overlaps <- overlaps[!duplicated(peakid)]

overlaps.rep <- overlaps[rep.dist == 0]
motif.input.list[["overlap.rep.kovwt.noifng.up"]] <- overlaps.rep[ko4H3_v_YR.l2fc > 0 & ko4H3_v_YR.padj < 0.05, .(chr, start, end, peakid)]
motif.input.list[["overlap.rep.kovwt.ifng.up"]] <- overlaps.rep[ko4H3.IFNG_v_YR.IFNG.l2fc > 0 & ko4H3.IFNG_v_YR.IFNG.padj < 0.05, .(chr, start, end, peakid)]

overlaps.rep.kdm5b <- overlaps[rep.dist == 0 & kdm5b.dist == 0]
motif.input.list[["overlap.rep.kdm5b.kovwt.noifng.up"]] <- overlaps.rep.kdm5b[ko4H3_v_YR.l2fc > 0 & ko4H3_v_YR.padj < 0.05, .(chr, start, end, peakid)]
motif.input.list[["overlap.rep.kdm5b.kovwt.ifng.up"]] <- overlaps.rep.kdm5b[ko4H3.IFNG_v_YR.IFNG.l2fc > 0 & ko4H3.IFNG_v_YR.IFNG.padj < 0.05, .(chr, start, end, peakid)]

# Extraneous, these would have duplicate peakids because using repmask.loc rna
rep.atac.rna.merge[, c("chr", "start", "end") := as.data.table(do.call(rbind, strsplit(V10, "-")))]
motif.input.list[["overlap.rep.kovwt.rna.up.noifng.up"]] <- rep.atac.rna.merge[rna.l2fc.multi > 0 & ko4H3_v_YR.l2fc > 0 & ko4H3_v_YR.padj < 0.05, .(chr, start, end, V10)]
motif.input.list[["overlap.rep.kovwt.rna.up.noifng.down"]] <- rep.atac.rna.merge[rna.l2fc.multi > 0 & ko4H3_v_YR.l2fc < 0 & ko4H3_v_YR.padj < 0.05, .(chr, start, end, V10)]
motif.input.list[["overlap.rep.kovwt.rna.up.ifng.up"]] <- rep.atac.rna.merge[rna.l2fc.multi > 0 & ko4H3.IFNG_v_YR.IFNG.l2fc > 0 & ko4H3.IFNG_v_YR.IFNG.padj < 0.05, .(chr, start, end, V10)]
motif.input.list[["overlap.rep.kovwt.rna.down.ifng.down"]] <- rep.atac.rna.merge[rna.l2fc.multi > 0 & ko4H3.IFNG_v_YR.IFNG.l2fc < 0 & ko4H3.IFNG_v_YR.IFNG.padj < 0.05, .(chr, start, end, V10)]

# Print
for(comp in names(motif.input.list)){
  tmp <- motif.input.list[[comp]]
  tmp <- tmp[!duplicated(tmp)]
  fwrite(tmp, paste0(comp, ".bed"), col.names = F, sep = "\t")
}

#### For Motif analysis ####

#### Combined analyses ####
# RE in increased peaks
res.dt <- fread("../../dat/atac_results_library.sizefactor.txt")
res.dt[, ko4H3_v_YR.neglogp := -log10(ko4H3_v_YR.padj)]
res.dt[, ko4H3.IFNG_v_YR.IFNG.neglogp := -log(ko4H3.IFNG_v_YR.IFNG.padj)]
repmask.dict <- fread("../../dat/repmask.dfam.final.dict.txt")

overlaps <- fread("../../dat/allYRandall4H3-ATAC.repmask.YR-KDM5B.closest.bed")
overlaps <- merge(overlaps[,.(peakid = V4, repmask.name = V8, repmask.name.loci = paste0(V8, "_", V5, "-", V6, "-", V7), rep.dist = V9, kdm5b.dist = V14)],
                  res.dt, by = "peakid", all.x = T)
overlaps <- merge(overlaps, repmask.dict, by = "repmask.name")
overlaps <- overlaps[order(ko4H3_v_YR.padj)]
overlaps.mmvl30 <- overlaps[!duplicated(peakid) & rep.dist == 0 & ko4H3_v_YR.padj < 0.05 & repmask.name == "MMVL30-int"]
overlaps.unique.peaks <- overlaps[!duplicated(peakid), .(peakid, ko4H3_v_YR.l2fc, ko4H3_v_YR.padj, ko4H3.IFNG_v_YR.IFNG.l2fc, ko4H3.IFNG_v_YR.IFNG.padj, ko4H3_v_YR.neglogp, ko4H3.IFNG_v_YR.IFNG.neglogp)]

# Final figure for interferon
library(showtext)
font_add("Arial", "/Applications/Microsoft Word.app/Contents/Resources/DFonts/arial.ttf")
showtext_auto()
up.peaks <- nrow(overlaps.unique.peaks[ko4H3.IFNG_v_YR.IFNG.padj < 0.05 & ko4H3.IFNG_v_YR.IFNG.l2fc > 0])
down.peaks <- nrow(overlaps.unique.peaks[ko4H3.IFNG_v_YR.IFNG.padj < 0.05 & ko4H3.IFNG_v_YR.IFNG.l2fc < 0])
mmvl30.up.peaks <- nrow(overlaps.mmvl30[ko4H3.IFNG_v_YR.IFNG.padj < 0.05 & ko4H3.IFNG_v_YR.IFNG.l2fc > 0])
mmvl30.down.peaks <- nrow(overlaps.mmvl30[ko4H3.IFNG_v_YR.IFNG.padj < 0.05 & ko4H3.IFNG_v_YR.IFNG.l2fc < 0])
mp <- ggplot(overlaps, aes(ko4H3.IFNG_v_YR.IFNG.l2fc, ko4H3.IFNG_v_YR.IFNG.neglogp)) +
  geom_vline(xintercept = 0, linetype="dashed", color = "gray") +
  geom_point(data = subset(overlaps, ko4H3.IFNG_v_YR.IFNG.padj >= 0.05), color = "gray", size = 0.2) +
  geom_point(data = subset(overlaps, ko4H3.IFNG_v_YR.IFNG.l2fc < 0 & ko4H3.IFNG_v_YR.IFNG.padj < 0.05), color = "#0000FF", size = 0.2) +
  geom_point(data = subset(overlaps, ko4H3.IFNG_v_YR.IFNG.l2fc > 0 & ko4H3.IFNG_v_YR.IFNG.padj < 0.05), color = "#F94040", size = 0.2) +
  # geom_point(data = subset(overlaps, !is.na(display.name)), color = "black") +
  geom_point(data = subset(overlaps.mmvl30, rep.dist == 0 & ko4H3.IFNG_v_YR.IFNG.padj < 0.05 & repmask.name == "MMVL30-int"), color = "black", size = 1, pch = 21, stroke = 1) +
  # geom_label_repel(aes(label = display.name), max.overlaps = 2000, force = 50, size = 2.5, nudge_x = -1.5, nudge_y = 5, segment.size = 0.25) +
  annotate("text", x = -5, y = 145, label = paste0("n = ", comma_format()(up.peaks)), color = "#F94040", hjust = 0, size = 4.6) +
  annotate("text", x = -5, y = 135, label = paste0("n = ", comma_format()(down.peaks)), color = "#0000FF", hjust = 0, size = 4.6) +
  # annotate("text", x = -3.1, y = 104, label = paste0("MMVL30 n = ", mmvl30.up.peaks), color = "red", hjust = 0) +
  # annotate("text", x = -3.1, y = 96, label = paste0("MMVL30 n = ", mmvl30.down.peaks), color = "blue", hjust = 0) +
  xlim(c(-5, 7.5)) +
  ylab(label = expression(-Log[10]*"(p value)")) +
  xlab(label = expression(Log[2]*"(KO IFN"*gamma*" vs WT IFN"*gamma*")")) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16),
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1),
        text=element_text(family="Arial"))
ggsave("ko4H3.IFNG_v_YR.IFNG.volcano.png", mp, width = 3.5, height = 4, dpi = 1200)

# Boxplots
overlaps.unique.rep <- overlaps[rep.dist == 0]
overlaps.unique.rep.type <- overlaps.unique.rep[!duplicated(overlaps.unique.rep[,.(peakid, final.type)])]
overlaps.unique.rep.type[, final.type := factor(final.type, c("LTR", "LINE", "SINE", "DNA"))]
mp <- ggplot(overlaps.unique.rep.type[ko4H3.IFNG_v_YR.IFNG.padj < 0.05 & final.type %in% c("LTR", "DNA", "SINE", "LINE")], aes(final.type, ko4H3.IFNG_v_YR.IFNG.l2fc)) +
  geom_hline(yintercept = 0, linetype="dashed", color = "gray") +
  geom_boxplot(outlier.shape = NA, width = 0.5, aes(fill = final.type), show.legend = F) +
  geom_jitter(width = 0.1, size = 0.1, alpha = 0.5) +
  ylab(expression(atop("ATAC-seq", paste(Log[2]*"(KO IFN"*gamma*" vs WT IFN"*gamma*")")))) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1),
        text=element_text(family="Arial"))
ggsave("ko4H3.IFNG_v_YR.IFNG.boxplot.repeat.type.png", mp, width = 4, height = 4, dpi = 1200)

#### combined analyses ####


