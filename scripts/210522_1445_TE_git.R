# Last Updated: 
# Author: Wesley Cai
# Purpose:

library(data.table)
library(ggplot2)
library(ggrepel)

# Colors: ERV1, ERVK, ERVL, LINE, SINE: "#FF64B0", "#F8766D", "#DE8C00", "#7CAE00", "#00BFC4"

# OG: Google_Drive/medschool/research/qin/rna/TE/210522_1445_TE.R
mainwd <- "."
outputfolder <- "output/"
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

repmask.dict <- fread("../../dat/repmask.dfam.final.dict.txt")

#### Non-loci ####
load("../../dat/1445.res.merged.repmask.RData")
res.merged.subset <- as.data.table(res.merged)
res.merged.subset <- res.merged.subset[,.(repmask.name = row,
                                          padj.sg_vs_sgCtrl, neglogp.sg_vs_sgCtrl, l2fc.sg_vs_sgCtrl)]
res.merged.subset <- res.merged.subset[order(padj.sg_vs_sgCtrl)]
fwrite(res.merged.subset, "1445.sgKDM5B_v_sgCtrl.results.txt", sep = "\t")
res.merged.subset <- fread("1445.sgKDM5B_v_sgCtrl.results.txt")

# res.merged.subset <- res.merged.subset[complete.cases(res.merged.subset)]
res.merged.subset <- merge(res.merged.subset, repmask.dict, by = "repmask.name")
unique(res.merged.subset$final.type)
unique(res.merged.subset$final.subtype)
unique(res.merged.subset$ste.type)
unique(res.merged.subset$ste.subtype)
res.merged.subset.te <- res.merged.subset[final.type %in% c("LTR", "LINE", "SINE", "DNA")]

res.merged.subset.te <- res.merged.subset.te[order(padj.sg_vs_sgCtrl)]
fwrite(res.merged.subset.te[,.(repmask.name, l2fc.1445.sgKDM5B_v_sgCtrl = l2fc.sg_vs_sgCtrl, padj.1445.sgKDM5B_v_sgCtrl = padj.sg_vs_sgCtrl, 
                               final.type, final.subtype)], "1445.sgKDM5B_v_sgCtrl.results.TEonly.txt", sep = "\t")

res.merged.subset.te <- res.merged.subset.te[complete.cases(res.merged.subset.te)]
res.merged.subset.te[1:1, label := repmask.name]
mp <- ggplot(res.merged.subset.te[padj.sg_vs_sgCtrl < 1], aes(l2fc.sg_vs_sgCtrl, neglogp.sg_vs_sgCtrl)) +
  geom_vline(xintercept = 0, linetype="dashed", color = "gray") +
  geom_point(data = subset(res.merged.subset.te, padj.sg_vs_sgCtrl >= 0.05), color = "gray", size = 0.2) +
  geom_point(data = subset(res.merged.subset.te, l2fc.sg_vs_sgCtrl < 0 & padj.sg_vs_sgCtrl < 0.05), color = "blue", size = 0.2) +
  geom_point(data = subset(res.merged.subset.te, l2fc.sg_vs_sgCtrl > 0 & padj.sg_vs_sgCtrl < 0.05), color = "red", size = 0.2) +
  geom_label_repel(aes(label = label), max.overlaps = 2000, force = 10, size = 2.5, nudge_y = 1, segment.size = 0.25, min.segment.length = 0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
mp
ggsave("1445.sgKDM5B_v_sgCtrl.nonloci.volcano.pdf", mp, width = 4, height = 4)

res.merged.subset.te.box <- res.merged.subset.te
res.merged.subset.te.box[final.subtype %in% c("ERV1", "ERVK", "ERVL"), te.label := final.subtype]
res.merged.subset.te.box[final.type %in% c("LINE", "SINE"), te.label := final.type]
res.merged.subset.te.box[,te.label := factor(te.label, c("ERV1", "ERVK", "ERVL", "LINE", "SINE"))]
res.merged.subset.te.box <- res.merged.subset.te.box[!is.na(te.label)]
mp <- ggplot(res.merged.subset.te.box, aes(te.label, l2fc.sg_vs_sgCtrl, fill = te.label)) +
  geom_hline(yintercept = 0, linetype="dashed", color = "gray") +
  geom_violin(width = 1, lwd = 0.5) +
  # geom_jitter(aes(size = neglogp.sg_vs_sgCtrl, alpha = 1-padj.sg_vs_sgCtrl), width = 0.05) +
  geom_jitter(aes(size = neglogp.sg_vs_sgCtrl), width = 0.05) +
  geom_label_repel(aes(label = label), max.overlaps = 2000, force = 50,
                   size = 4, nudge_x = 1, segment.size = 0.25, fill = "white",
                   show.legend = F) +
  guides(fill = guide_legend(override.aes = list(shape = NA)),
         alpha = guide_none()) +
  scale_fill_manual(name = "TE", values = c("#FF64B0", "#F8766D", "#DE8C00", 
                                            "#7CAE00", "#00BFC4")) +
  scale_size(name = "-log P") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        panel.border = element_rect(fill=NA, colour = "black", size=1),
        axis.ticks = element_line(colour = "black", size = 1))
mp
# ERV1, ERVK, ERVL, LINE, SINE: "#FF64B0", "#F8766D", "#DE8C00", "#7CAE00", "#00BFC4"

ggsave("1445.sgKDM5B_v_sgCtrl.nonloci.boxplot.fullalpha.pdf", mp, width = 6, height = 5)
#### Non-loci ####
library(scales)
show_col(hue_pal()(4))
show_col(hue_pal()(12))

#### Loci ####
# File can be found at GSE161065: Superseries
load("../../dat/1445.res.merged.repmask.loc.RData")
res.merged.subset <- as.data.table(res.merged)
res.merged.subset[, repmask.name := gsub("_chr.*", "", row)]
res.merged.subset <- res.merged.subset[,.(repmask.name, repmask.name.loci = row,
                                          padj.1445.sgKDM5B_v_sgCtrl = padj.sg_vs_sgCtrl,
                                          neglogp.sg_vs_sgCtrl,
                                          l2fc.1445.sgKDM5B_v_sgCtrl = l2fc.sg_vs_sgCtrl)]
res.merged.subset <- res.merged.subset[complete.cases(res.merged.subset)]
res.merged.subset <- res.merged.subset[order(padj.sg_vs_sgCtrl)]
res.merged.subset[, repmask.name.loci.num := paste0(repmask.name, "_loci_", seq(1, .N)), by = .(repmask.name)]
res.merged.subset <- merge(res.merged.subset, repmask.dict, by = "repmask.name")
res.merged.subset.te <- res.merged.subset[final.type %in% c("LTR", "LINE", "SINE", "DNA")]

res.merged.subset.te <- res.merged.subset.te[order(padj.sg_vs_sgCtrl)]
res.merged.subset.te[1:1, label := repmask.name.loci.num]
mp <- ggplot(res.merged.subset.te[padj.sg_vs_sgCtrl < 1], aes(l2fc.sg_vs_sgCtrl, neglogp.sg_vs_sgCtrl)) +
  geom_vline(xintercept = 0, linetype="dashed", color = "gray") +
  geom_point(data = subset(res.merged.subset.te, padj.sg_vs_sgCtrl >= 0.05), color = "gray", size = 0.2) +
  geom_point(data = subset(res.merged.subset.te, l2fc.sg_vs_sgCtrl < 0 & padj.sg_vs_sgCtrl < 0.05), color = "blue", size = 0.2) +
  geom_point(data = subset(res.merged.subset.te, l2fc.sg_vs_sgCtrl > 0 & padj.sg_vs_sgCtrl < 0.05), color = "red", size = 0.2) +
  geom_label_repel(aes(label = label), max.overlaps = 2000, force = 50, size = 2.5, nudge_x = -1.5, nudge_y = 5, segment.size = 0.25) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
mp
ggsave("1445.sgKDM5B_v_sgCtrl.volcano.pdf", mp, width = 4, height = 4)

res.merged.subset.te.box <- res.merged.subset.te
res.merged.subset.te.box[final.subtype %in% c("ERV1", "ERVK", "ERVL"), te.label := final.subtype]
res.merged.subset.te.box[final.type %in% c("LINE", "SINE"), te.label := final.type]
res.merged.subset.te.box[,te.label := factor(te.label, c("ERV1", "ERVK", "ERVL", "LINE", "SINE"))]
res.merged.subset.te.box <- res.merged.subset.te.box[!is.na(te.label)]
mp <- ggplot(res.merged.subset.te.box, aes(te.label, l2fc.sg_vs_sgCtrl, fill = te.label)) +
  geom_hline(yintercept = 0, linetype="dashed", color = "gray") +
  # geom_violin(width = 1) +
  geom_boxplot(width = 0.5) +
  geom_jitter(aes(size = neglogp.sg_vs_sgCtrl, alpha = 1-padj.sg_vs_sgCtrl), width = 0.05) +
  geom_label_repel(aes(label = label), max.overlaps = 2000, force = 50,
                   size = 4, nudge_x = 1, segment.size = 0.25, fill = "white",
                   show.legend = F) +
  guides(fill = guide_legend(override.aes = list(shape = NA)),
         alpha = guide_none()) +
  scale_fill_discrete(name = "TE") +
  scale_size(name = "-log P") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        panel.border = element_rect(fill=NA, colour = "black", size=1),
        axis.ticks = element_line(colour = "black", size = 1))
mp
ggsave("1445.sgKDM5B_v_sgCtrl.boxplot.pdf", mp, width = 6, height = 5)
#### Loci ####
