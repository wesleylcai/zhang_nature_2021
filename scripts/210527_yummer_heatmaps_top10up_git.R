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

# OG: Google_Drive/medschool/research/qin/chip/analysis/210527_yummer_heatmaps/210527_yummer_heatmaps_top10up.R
mainwd <- "."
outputfolder <- "output/"
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

source("common/metaplot_heatmap_functions.R")
# File can be found at GSE161065: Superseries
repmask <- fread("../../dat/mm10.repmask.s.bed")
repmask.dict <- fread("../../dat/repmask.dfam.final.dict.txt")
load("../../dat/yummer.res.merged.nonloci.subset.RData")
# File can be found at GSE175542: Yummer ChIP-seq
load("../../dat/yummer.rep.signals.subtractinput.rpkm.RData")

#### Subset data for heatmap ####
res.merged.subset <- res.merged.subset[order(padj.sg_vs_sgCtrl)]
bigwigs.dt <- data.table(track.name = c("KDM5B", "SETDB1", "H3K9me3", "H3K4me3", "H3K4me2", "Input"))
signal.matrix <- as.data.table(do.call(cbind, rep.signals.subtractinput.rpkm))
signal.matrix.pos <- signal.matrix[,.(position = `7SLRNA.position`)]
signal.matrix <- transpose(signal.matrix[,.SD,.SDcols = !grepl("position", colnames(signal.matrix))], keep.names = "repmask.name")
test <- signal.matrix[,.SD,.SDcols = 1:5]
head(test)
signal.matrices <- list()
target <- "KDM5B"
for(target in bigwigs.dt$track.name){
  tmp.ma <- signal.matrix[grepl(paste0(".", target), repmask.name)]
  tmp.ma[, repmask.name := gsub(paste0(".", target), "", repmask.name)]
  head(tmp.ma[,.SD,.SDcols = 1:5])
  if(target == "KDM5B"){
    tmp.ma <- merge(tmp.ma, res.merged.subset[,.(repmask.name, l2fc.sg_vs_sgCtrl, neglogp.sg_vs_sgCtrl, te.label)], by = "repmask.name")
    tmp.ma <- tmp.ma[,.SD,.SDcols = c("repmask.name",
                                      "l2fc.sg_vs_sgCtrl",
                                      "neglogp.sg_vs_sgCtrl",
                                      "te.label",
                                      paste0("V", 3000:7000))]
    head(tmp.ma[,.SD,.SDcols = 1:5])
  } else {
    tmp.ma <- tmp.ma[,.SD,.SDcols = c("repmask.name",
                                      paste0("V", 3000:7000))]
  }
  signal.matrices[[target]] <- tmp.ma
}
head(signal.matrices[["KDM5B"]][,.SD,.SDcols = c("repmask.name", paste0("V", 3000:3010))])
head(signal.matrices[["SETDB1"]][,.SD,.SDcols = c("repmask.name", paste0("V", 3000:3010))])
#### Subset data for heatmap ####

#### Draw heatmap ####
kdm5b.sorted <- signal.matrices[["KDM5B"]]
kdm5b.sorted[, signalsums := rowSums(kdm5b.sorted[,.SD,.SDcols = paste0("V", 3000:7000)])]
unique(kdm5b.sorted$te.label)
kdm5b.sorted[, te.label := factor(te.label, c("ERV1", "ERVK", "ERVL", "LINE", "SINE"))]

kdm5b.sorted <- kdm5b.sorted[order(signalsums, decreasing = T),]
kdm5b.sorted <- kdm5b.sorted[order(te.label),]
sort.order <- kdm5b.sorted$repmask.name
# plotSorted(sort.order = sort.order, sort.name = "TE_sorted") #### SEE BELOW

kdm5b.sorted <- kdm5b.sorted[order(signalsums, decreasing = T),]
sort.order <- kdm5b.sorted$repmask.name
# plotSorted(sort.order = sort.order, sort.name = "KDM5B_sorted") #### SEE BELOW

#### Cannot make a function because of the "draw" function so need to change variables and run all code below

sort.name <- "TE_sorted"
sort.name <- "KDM5B_sorted"

signal.matrices.sorted <- list()
for(target in names(signal.matrices)){
  tmp.ma <- signal.matrices[[target]]
  setkey(tmp.ma, "repmask.name")
  signal.matrices.sorted[[target]] <- tmp.ma[sort.order,]
}
head(signal.matrices.sorted[["KDM5B"]][,.SD,.SDcols = c("repmask.name", paste0("V", 3000:3010))])
head(signal.matrices.sorted[["SETDB1"]][,.SD,.SDcols = c("repmask.name", paste0("V", 3000:3010))])

annote.dt <- signal.matrices.sorted[["KDM5B"]][,.SD,.SDcols = c("repmask.name",
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
median(as.matrix(kdm5b.sorted[,.SD,.SDcols = paste0("V", 3000:7000)]))

matrices.dt <- data.table(mat = grep("Input", names(signal.matrices.sorted), invert = T, value = T),
                          color = c("red", "dodgerblue3", "forestgreen", "darkorchid", "black"))
# For some reason this doesn't work
# ht.list <- list()
# for(rown in 1:nrow(matrices.dt)){
#   target.value <- matrices.dt[rown, get("mat")]
#   color.value <- matrices.dt[rown, get("color")]
#   ht.list[[target.value]] <- Heatmap(as.matrix(signal.matrices.sorted[[target]][,.SD,.SDcols = paste0("V", 3000:7000)]), name = target.value, column_title = target.value, cluster_rows = FALSE, cluster_columns = FALSE,
#                                      show_column_names = FALSE, show_row_names = FALSE, col = colorRamp2(c(0, 5), c("white", color.value)), use_raster = T)
# }

start.value <- 3000
end.value <- 7000
# start.value <- 4500
# end.value <- 5500
ht.kdm5b <- Heatmap(as.matrix(signal.matrices.sorted[["KDM5B"]][,.SD,.SDcols = paste0("V", start.value:end.value)]), name = "KDM5B", column_title = "KDM5B", cluster_rows = FALSE, cluster_columns = FALSE,
                    show_column_names = FALSE, show_row_names = FALSE, col = colorRamp2(c(0, 5), c("white", "red")), use_raster = T)
ht.setdb1 <- Heatmap(as.matrix(signal.matrices.sorted[["SETDB1"]][,.SD,.SDcols = paste0("V", start.value:end.value)]), name = "SETDB1", column_title = "SETDB1", cluster_rows = FALSE, cluster_columns = FALSE,
                     show_column_names = FALSE, show_row_names = FALSE, col = colorRamp2(c(0, 5), c("white", "dodgerblue3")), use_raster = T)
ht.h3k9me3 <- Heatmap(as.matrix(signal.matrices.sorted[["H3K9me3"]][,.SD,.SDcols = paste0("V", start.value:end.value)]), name = "H3K9me3", column_title = "H3K9me3", cluster_rows = FALSE, cluster_columns = FALSE,
                      show_column_names = FALSE, show_row_names = FALSE, col = colorRamp2(c(0, 5), c("white", "forestgreen")), use_raster = T)
ht.h3k4me3 <- Heatmap(as.matrix(signal.matrices.sorted[["H3K4me3"]][,.SD,.SDcols = paste0("V", start.value:end.value)]), name = "H3K4me3", column_title = "H3K4me3", cluster_rows = FALSE, cluster_columns = FALSE,
                      show_column_names = FALSE, show_row_names = FALSE, col = colorRamp2(c(0, 5), c("white", "darkorchid")), use_raster = T)
ht.h3k4me2 <- Heatmap(as.matrix(signal.matrices.sorted[["H3K4me2"]][,.SD,.SDcols = paste0("V", start.value:end.value)]), name = "H3K4me2", column_title = "H3K4me2", cluster_rows = FALSE, cluster_columns = FALSE,
                      show_column_names = FALSE, show_row_names = FALSE, col = colorRamp2(c(0, 5), c("white", "black")), use_raster = T)


# repmask.goi <- annote.dt[(L2FC > 0 & nlogP > -log10(5e-10)) | name %in% c("MMERVK10D3_I-int", "MMERVK10C-int",
#                                                                           "RLTR6B_Mm", "RLTR6-int"), get("name")]

## Top 10 up
res.merged.subset <- res.merged.subset[order(padj.sg_vs_sgCtrl)]
up.dt <- res.merged.subset[l2fc.sg_vs_sgCtrl > 0,]
top10 <- up.dt[1:10, get("repmask.name")]
# repmask.goi <- annote.dt[name %in% c(top10, "MMERVK10D3_I-int", "MMERVK10C-int", "RLTR6B_Mm", "RLTR6-int"), get("name")]
repmask.goi <- annote.dt[name %in% c(top10), get("name")]

labelSubset <- which(annote.dt$name %in% repmask.goi)
labels <- annote.dt[labelSubset, get("name")]
rowlab <- rowAnnotation(link = anno_mark(at = labelSubset, labels = labels)) #, width = unit(1, "cm") + max_text_width(labels))

ht_list <- ha + ht.kdm5b + rowlab
pdf(paste0("yummer.KDM5B.heatmap.top10only.labelINCREASED.", sort.name,".pdf"), width = 6, height = 10)
draw(ht_list, show_heatmap_legend = T, show_annotation_legend = F)
decorate_heatmap_body(c("KDM5B"), {
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
})
dev.off()

ht_list <- ha + ht.kdm5b + ht.setdb1 + ht.h3k9me3 + rowlab
pdf(paste0("yummer.KDM5B.SETDB1.H3K9me3.heatmap.top10only.labelINCREASED.", sort.name,".pdf"), width = 9, height = 10)
draw(ht_list, show_heatmap_legend = T, show_annotation_legend = F)
decorate_heatmap_body(c("KDM5B"), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))})
decorate_heatmap_body(c("SETDB1"), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))})
decorate_heatmap_body(c("H3K9me3"), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))})
dev.off()

ht_list <- ha + ht.kdm5b + ht.setdb1 + ht.h3k9me3 + ht.h3k4me3 + ht.h3k4me2 + rowlab
pdf(paste0("yummer.KDM5B.SETDB1.H3K9me3.H3K4me3.H3K4me2.heatmap.top10only.labelINCREASED.", sort.name,".pdf"), width = 11, height = 10)
draw(ht_list, show_heatmap_legend = T, show_annotation_legend = F)
decorate_heatmap_body(c("KDM5B"), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))})
decorate_heatmap_body(c("SETDB1"), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))})
decorate_heatmap_body(c("H3K9me3"), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))})
decorate_heatmap_body(c("H3K4me3"), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))})
decorate_heatmap_body(c("H3K4me2"), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))})
dev.off()

## Top 10 plus
res.merged.subset <- res.merged.subset[order(padj.sg_vs_sgCtrl)]
up.dt <- res.merged.subset[l2fc.sg_vs_sgCtrl > 0,]
top10 <- up.dt[1:10, get("repmask.name")]
repmask.goi <- annote.dt[name %in% c(top10, "MMERVK10D3_I-int", "MMERVK10C-int", "RLTR6B_Mm", "RLTR6-int"), get("name")]
# repmask.goi <- annote.dt[name %in% c(top10), get("name")]

labelSubset <- which(annote.dt$name %in% repmask.goi)
labels <- annote.dt[labelSubset, get("name")]
rowlab <- rowAnnotation(link = anno_mark(at = labelSubset, labels = labels)) #, width = unit(1, "cm") + max_text_width(labels))

ht_list <- ha + ht.kdm5b + rowlab
pdf(paste0("yummer.KDM5B.heatmap.top10plus.labelINCREASED.", sort.name,".pdf"), width = 6, height = 10)
draw(ht_list, show_heatmap_legend = T, show_annotation_legend = F)
decorate_heatmap_body(c("KDM5B"), {
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
})
dev.off()

ht_list <- ha + ht.kdm5b + ht.setdb1 + ht.h3k9me3 + rowlab
pdf(paste0("yummer.KDM5B.SETDB1.H3K9me3.heatmap.top10plus.labelINCREASED.", sort.name,".pdf"), width = 9, height = 10)
draw(ht_list, show_heatmap_legend = T, show_annotation_legend = F)
decorate_heatmap_body(c("KDM5B"), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))})
decorate_heatmap_body(c("SETDB1"), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))})
decorate_heatmap_body(c("H3K9me3"), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))})
dev.off()

ht_list <- ha + ht.kdm5b + ht.setdb1 + ht.h3k9me3 + ht.h3k4me3 + ht.h3k4me2 + rowlab
pdf(paste0("yummer.KDM5B.SETDB1.H3K9me3.H3K4me3.H3K4me2.heatmap.top10plus.labelINCREASED.", sort.name,".pdf"), width = 11, height = 10)
draw(ht_list, show_heatmap_legend = T, show_annotation_legend = F)
decorate_heatmap_body(c("KDM5B"), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))})
decorate_heatmap_body(c("SETDB1"), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))})
decorate_heatmap_body(c("H3K9me3"), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))})
decorate_heatmap_body(c("H3K4me3"), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))})
decorate_heatmap_body(c("H3K4me2"), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))})
dev.off()

## Top 10 down
res.merged.subset <- res.merged.subset[order(padj.sg_vs_sgCtrl)]
up.dt <- res.merged.subset[l2fc.sg_vs_sgCtrl < 0,]
top10 <- up.dt[1:10, get("repmask.name")]
# repmask.goi <- annote.dt[name %in% c(top10, "MMERVK10D3_I-int", "MMERVK10C-int", "RLTR6B_Mm", "RLTR6-int"), get("name")]
repmask.goi <- annote.dt[name %in% c(top10), get("name")]

labelSubset <- which(annote.dt$name %in% repmask.goi)
labels <- annote.dt[labelSubset, get("name")]
rowlab <- rowAnnotation(link = anno_mark(at = labelSubset, labels = labels)) #, width = unit(1, "cm") + max_text_width(labels))

ht_list <- ha + ht.kdm5b + rowlab
pdf(paste0("yummer.KDM5B.heatmap.top10only.labelDECREASED.", sort.name,".pdf"), width = 6, height = 10)
draw(ht_list, show_heatmap_legend = T, show_annotation_legend = F)
decorate_heatmap_body(c("KDM5B"), {
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
})
dev.off()
#### Draw heatmap ####


