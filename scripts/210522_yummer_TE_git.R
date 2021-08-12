# Last Updated: 
# Author: Wesley Cai
# Purpose:

library(data.table)
library(ggplot2)
library(ggrepel)
library(DESeq2)

# Colors: ERV1, ERVK, ERVL, LINE, SINE: "#FF64B0", "#F8766D", "#DE8C00", "#7CAE00", "#00BFC4"
# OG file: Google_Drive/medschool/research/qin/rna/TE/210522_yummer_TE_edited.R

mainwd <- "."
outputfolder <- "output/"
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

repmask.dict <- fread("../../dat/repmask.dfam.final.dict.txt")

#### non-loci ####
load("../../dat/yummer.res.merged.nonloci.RData")
res.merged <- as.data.table(res.merged)
res.merged.subset <- merge(res.merged[,.(repmask.name = row, l2fc.yummer.sgKDM5B_v_sgCtrl = log2FoldChange.sgKDM5B_vs_sgCtrl,
                                         padj.yummer.sgKDM5B_v_sgCtrl = padj.sgKDM5B_vs_sgCtrl)], repmask.dict, by = "repmask.name")
res.merged.subset.te <- res.merged.subset[final.type %in% c("LTR", "LINE", "SINE", "DNA")]
fwrite(res.merged.subset.te, "yummer.sgKDM5B_v_sgCtrl.results.TEonly.txt", sep = "\t")
#### non-loci ####

#### loci ####
load("../../dat/yummer.loc.res.merged.RData")
res.merged.subset <- as.data.table(res.merged)
res.merged.subset[, repmask.name := gsub("_chr.*", "", row)]
res.merged.subset <- res.merged.subset[,.(repmask.name, repmask.name.loci = row, rna.l2fc = log2FoldChange.sgKDM5B_vs_sgCtrl,
                                          rna.padj = padj.sgKDM5B_vs_sgCtrl, rna.neglogp = neglogp.sgKDM5B_vs_sgCtrl)]
res.merged.subset <- res.merged.subset[order(rna.padj)]
res.merged.subset <- res.merged.subset[complete.cases(res.merged.subset)]
res.merged.subset[, repmask.name.loci.num := paste0(repmask.name, "_loci_", seq(1, .N)), by = .(repmask.name)]
res.merged.subset <- merge(res.merged.subset, repmask.dict, by = "repmask.name")
unique(res.merged.subset$final.type)
unique(res.merged.subset$final.subtype)
unique(res.merged.subset$ste.type)
unique(res.merged.subset$ste.subtype)
res.merged.subset.te <- res.merged.subset[final.type %in% c("LTR", "LINE", "SINE", "DNA")]

res.merged.subset.te <- res.merged.subset.te[order(rna.padj)]
res.merged.subset.te[1:10, label := repmask.name]
res.merged.subset.te[!is.na(label), color := "black"]
res.merged.subset.te[!is.na(label) & repmask.name == "MMVL30-int", color := "red"]
mp <- ggplot(res.merged.subset.te, aes(rna.l2fc, rna.neglogp)) +
  geom_vline(xintercept = 0, linetype="dashed", color = "gray") +
  geom_point(data = subset(res.merged.subset.te, rna.padj >= 0.05), color = "gray", size = 0.4) +
  geom_point(data = subset(res.merged.subset.te, rna.l2fc < 0 & rna.padj < 0.05), color = "#0000FF", size = 0.4) +
  geom_point(data = subset(res.merged.subset.te, rna.l2fc > 0 & rna.padj < 0.05), color = "#F94040", size = 0.4) +
  geom_label_repel(aes(label = label, color = color), max.overlaps = 2000, force = 20, size = 3,
                   nudge_y = -10, segment.size = 0.25, min.segment.length = 0.5) +
  xlim(c(-10,10)) +
  scale_color_manual(values = c("#000000", "#F94040")) +
  guides(color = guide_none()) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.title.x = element_blank(),
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        panel.border = element_rect(fill=NA, colour = "black", size=1),
        axis.ticks = element_line(colour = "black", size = 1))
ggsave("yummer.sgKDM5B_v_sgCtrl.volcano.update.png", mp, width = 4.2, height = 4.5, dpi = 1200)

res.merged.subset.te.box <- res.merged.subset.te
res.merged.subset.te.box[final.subtype %in% c("ERV1", "ERVK", "ERVL"), te.label := final.subtype]
res.merged.subset.te.box[final.subtype %in% c("ERVL-MaLR"), te.label := "ERVL"]
res.merged.subset.te.box[final.type %in% c("LINE", "SINE"), te.label := final.type]
res.merged.subset.te.box[,te.label := factor(te.label, c("ERV1", "ERVK", "ERVL", "LINE", "SINE"))]
res.merged.subset.te.box <- res.merged.subset.te.box[!is.na(te.label)]
pos <- position_jitter(width = 0.2, seed = 2)
mp <- ggplot(res.merged.subset.te.box[rna.padj < 0.05], aes(te.label, rna.l2fc, fill = te.label, size = rna.neglogp)) +
  geom_hline(yintercept = 0, linetype="dashed", color = "gray") +
  geom_violin(width = 1, lwd = 0.5) +
  # geom_violin(width = 1, lwd = 0.5, color = NA) +
  geom_point(alpha = 0.3, position = pos, color = "black") +
  # geom_point(alpha = 0.1, position = pos, color = "black") +
  # geom_point(alpha = 0, position = pos, color = "black") +
  geom_point(data = subset(res.merged.subset.te.box, !is.na(label)), aes(size = rna.neglogp), color = "black", position = pos) +
  geom_point(data = subset(res.merged.subset.te.box, !is.na(label)), aes(size = rna.neglogp), fill = "red", color = "black", position = pos, pch = 21, show.legend = F) +
  geom_label_repel(aes(label = label), max.overlaps = 2000, force = 50,
                   size = 5, segment.size = 0.5, fill = "white",
                   position = pos, show.legend = F) +
  guides(fill = guide_legend(override.aes = list(shape = NA)),
         alpha = guide_none()) +
  scale_fill_manual(name = "RE", values = c("#FF64B0", "#F8766D", "#DE8C00", 
                                            "#7CAE00", "#00BFC4")) +
  scale_size(name = "-log P") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        panel.border = element_rect(fill=NA, colour = "black", size=1),
        axis.ticks = element_line(colour = "black", size = 1))
# ERV1, ERVK, ERVL, LINE, SINE: "#FF64B0", "#F8766D", "#DE8C00", "#7CAE00", "#00BFC4"
ggsave("yummer.sgKDM5B_v_sgCtrl.boxplot.update.png", mp, width = 7.3, height = 5.5, dpi = 1200)
#### loci ####