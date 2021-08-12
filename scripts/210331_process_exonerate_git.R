# Last Updated: 
# Author: Wesley Cai
# Purpose:

library(data.table)
library(ggplot2)
library(reshape2)
library(spgs)
library(scales)
library(ggrepel)

# OG: Google_Drive/medschool/research/qin/rna/stringtie/210331_process_exonerate.R
mainwd <- "."
outputfolder <- "exonerate/"
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

tx.repmask.overlap.file <- "../../dat/yummer_sgCtrl_sgKDM5B.transcript.mm10.repmask.s.bed"
tx.deseq.file <- "../../dat/yummer_sgCtrl_sgKDM5B.stringtie.exon.results.txt"

res <- fread("../../dat/exonerate.output.tab")
res <- res[V3 != "revcomp"]
tx.deseq <- fread(tx.deseq.file)
tx.repmask.overlap <- fread(tx.repmask.overlap.file)
tx.repmask.overlap <- tx.repmask.overlap[,.(V1, V2, V3, V4, V11, V12, V13, V14, repelement = paste0(V11, ":", V12, "-", V13))]

res.repmask <- merge(res, tx.repmask.overlap[,.(repelement, str.1 = V4, first.repelement = V14)], by.x = "V1", by.y = "repelement")
res.repmask <- merge(res.repmask, tx.repmask.overlap[,.(repelement, str.2 = V4, second.repelement = V14)], by.x = "V2", by.y = "repelement")
res.repmask <- res.repmask[str.1 == str.2]
res.repmask.deseq <- merge(res.repmask, tx.deseq, by.x = "str.1", by.y = "symbol")
res.repmask.deseq[, neglogp := -log(sgKDM5B_v_sgctrl.padj)]
nrow(res.repmask.deseq[sgKDM5B_v_sgctrl.l2fc > 0])
nrow(res.repmask.deseq[sgKDM5B_v_sgctrl.l2fc < 0])
fwrite(res.repmask.deseq, "invertrepeat.repmask.deseq.txt", sep = "\t")
res.repmask.deseq <- fread("invertrepeat.repmask.deseq.txt")
res.repmask.deseq <- res.repmask.deseq[order(sgKDM5B_v_sgctrl.padj)]
res.repmask.deseq[1:10, label := paste0(first.repelement, ";", second.repelement)]

up.elements <- nrow(res.repmask.deseq[sgKDM5B_v_sgctrl.l2fc > 0 & sgKDM5B_v_sgctrl.padj < 0.05])
down.elements <- nrow(res.repmask.deseq[sgKDM5B_v_sgctrl.l2fc < 0 & sgKDM5B_v_sgctrl.padj < 0.05])
mp <- ggplot(res.repmask.deseq, aes(sgKDM5B_v_sgctrl.l2fc, neglogp)) +
  geom_vline(xintercept = 0, linetype="dashed", color = "gray") +
  geom_point(data = subset(res.repmask.deseq, sgKDM5B_v_sgctrl.padj >= 0.05), color = "gray", size = 0.5) +
  geom_point(data = subset(res.repmask.deseq, sgKDM5B_v_sgctrl.l2fc < 0 & sgKDM5B_v_sgctrl.padj < 0.05), color = "blue", size = 0.5) +
  geom_point(data = subset(res.repmask.deseq, sgKDM5B_v_sgctrl.l2fc > 0 & sgKDM5B_v_sgctrl.padj < 0.05), color = "red", size = 0.5) +
  # geom_point(data = subset(overlaps, !is.na(display.name)), color = "black") +
  geom_label_repel(aes(label = label), max.overlaps = 2000, force = 100, size = 2.5, nudge_x = 1, nudge_y = 5, segment.size = 0.25, min.segment.length = 0.1) +
  annotate("text", x = -4, y = 120, label = paste0("n = ", comma_format()(up.elements)), color = "red", hjust = 0) +
  annotate("text", x = -4, y = 112, label = paste0("n = ", comma_format()(down.elements)), color = "blue", hjust = 0) +
  # annotate("text", x = -3.1, y = 104, label = paste0("MMVL30 n = ", mmvl30.up.peaks), color = "red", hjust = 0) +
  # annotate("text", x = -3.1, y = 96, label = paste0("MMVL30 n = ", mmvl30.down.peaks), color = "blue", hjust = 0) +
  # xlim(c(-5, 7.5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
mp
ggsave("invertrepeat.repmask.volcano.pdf", mp, width = 4, height = 4)
ggsave("invertrepeat.repmask.volcano.png", mp, width = 4, height = 4, dpi = 1200)
