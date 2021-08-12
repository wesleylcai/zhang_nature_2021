# Last Updated: 
# Author: Wesley Cai
# Purpose:

library(ggplot2)
library(reshape2)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(matrixStats)

# Colors: ERV1, ERVK, ERVL, LINE, SINE: "#FF64B0", "#F8766D", "#DE8C00", "#7CAE00", "#00BFC4"

# OG: Google_Drive/medschool/research/qin/chip/analysis/210524_1445_heatmaps/210524_1445_heatmaps
mainwd <- "."
outputfolder <- "output/"
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

source("common/metaplot_heatmap_functions.R")

# File can be found at GSE161065
repmask <- fread("../../dat/mm10.repmask.s.bed")
repmask.dict <- fread("../../dat/repmask.dfam.final.dict.txt")

load("../../dat/1445.res.merged.repmask.RData")
res.merged.subset <- as.data.table(res.merged)
res.merged.subset <- res.merged.subset[,.(repmask.name = row,
                                          padj.sg_vs_sgCtrl, neglogp.sg_vs_sgCtrl, l2fc.sg_vs_sgCtrl)]
res.merged.subset[, neglogp.sg_vs_sgCtrl := -log10(padj.sg_vs_sgCtrl)]
res.merged.subset <- res.merged.subset[order(padj.sg_vs_sgCtrl)]
res.merged.subset <- res.merged.subset[!is.na(l2fc.sg_vs_sgCtrl)]
res.merged.subset <- merge(res.merged.subset, repmask.dict, by = "repmask.name")
res.merged.subset[final.subtype %in% c("ERV1", "ERVK", "ERVL"), te.label := final.subtype]
res.merged.subset[final.type %in% c("LINE", "SINE"), te.label := final.type]
res.merged.subset <- res.merged.subset[!is.na(te.label)]

#### Subset data for heatmap ####
# File can be found at GSE161062: 1445 ChIP-Seq
load("../../dat/1445.rep.signal.subtract.RData")
signal.matrix <- as.data.table(do.call(cbind, rep.signal.subtract))
signal.matrix.pos <- signal.matrix[,.(position = `7SLRNA.position`)]
signal.matrix <- transpose(signal.matrix[,.SD,.SDcols = !grepl("position", colnames(test))], keep.names = "repmask.name")
signal.matrix[, repmask.name := gsub("\\.Signal", "", repmask.name)]
head(signal.matrix[,.SD,.SDcols = 1:5])
signal.matrix.res <- merge(signal.matrix, res.merged.subset[,.(repmask.name, l2fc.sg_vs_sgCtrl, neglogp.sg_vs_sgCtrl, te.label)], by = "repmask.name")

signal.matrix.subset <- signal.matrix.res[,.SD,.SDcols = c("repmask.name",
                                                       "l2fc.sg_vs_sgCtrl",
                                                       "neglogp.sg_vs_sgCtrl",
                                                       "te.label",
                                                       paste0("V", 3000:7000))]
head(signal.matrix.subset[,.SD,.SDcols = 1:5])
signal.matrix.subset[, signalsums := rowSums(signal.matrix.subset[,.SD,.SDcols = paste0("V", 3000:7000)])]
signal.matrix.subset <- signal.matrix.subset[order(signalsums, decreasing = T),]
head(signal.matrix.subset[,.SD,.SDcols = 1:5])
fwrite(signal.matrix.subset[,.SD,.SDcols = 1:4], "heatmap.order.txt", sep = "\t")
#### Subset data for heatmap ####

#### Draw heatmap ####
annote.dt <- signal.matrix.subset[,.SD,.SDcols = c("repmask.name",
                                                   "l2fc.sg_vs_sgCtrl",
                                                   "neglogp.sg_vs_sgCtrl",
                                                   "te.label")]
colnames(annote.dt) <- c("name", "L2FC", "nlogP", "TE_class")

ha = rowAnnotation(df = annote.dt[,.SD,.SDcols = c("L2FC", "nlogP", "TE_class")],
                   annotation_width = unit(c(5, 5, 5), "mm"),
                   col = list(TE_class = c("SINE" = "#00BFC4", "ERVK" = "#F8766D", "LINE" = "#7CAE00", "ERV1" = "#FF64B0", "ERVL" = "#DE8C00"),
                              # nlogP = colorRamp2(c(0, max(annote.dt$nlogP, na.rm = TRUE)), c("white", "dodgerblue4")),
                              nlogP = colorRamp2(c(0, max(annote.dt$nlogP, na.rm = TRUE)), c("white", "black")),
                              L2FC = colorRamp2(c(-max(abs(annote.dt$L2FC), na.rm = TRUE), 0, max(abs(annote.dt$L2FC), na.rm = TRUE)), c("green3", "white", "red2"))
                   )
)

#We can do significant ones, the ones I want to highlight particularly  are mmervk10d3, x2-Line mmervk10c,rltr6b_Mm, rltr6-int 
ht <- Heatmap(as.matrix(signal.matrix.subset[,.SD,.SDcols = paste0("V", 3000:7000)]), name = "Signal", column_title = "KDM5B", cluster_rows = FALSE, cluster_columns = FALSE,
              show_column_names = FALSE, show_row_names = FALSE, col = colorRamp2(c(0, 20), c("white", "red")), use_raster = T)

repmask.goi <- annote.dt[(L2FC > 0 & nlogP > -log10(0.05)) | name %in% c("MMERVK10D3_I-int", "MMERVK10C-int",
                                                            "RLTR6B_Mm", "RLTR6-int"), get("name")]
labelSubset <- which(annote.dt$name %in% repmask.goi)
labels <- annote.dt[labelSubset, get("name")]
rowlab <- rowAnnotation(link = anno_mark(at = labelSubset, labels = labels)) #, width = unit(1, "cm") + max_text_width(labels))

ht_list <- ha + ht + rowlab

pdf("1445.KDM5B.heatmap.pdf", width = 6, height = 10)
draw(ht_list, show_heatmap_legend = T, show_annotation_legend = F)
decorate_heatmap_body(c("Signal"), {
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
})
dev.off()
#### Draw heatmap ####


