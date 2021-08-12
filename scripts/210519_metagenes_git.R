# Last Updated: 
# Author: Wesley Cai
# Purpose:

library(data.table)
library(ggplot2)
library(reshape2)

# OG: Google_Drive/medschool/research/qin/chip/analysis/210519_metagenes/210519_metagenes.R
mainwd <- "."
outputfolder <- "output/"
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

source("../common/metaplot_heatmap_functions.R")
repmask.rna.dict <- fread("../../dat/mm10.repmask.rna.dict.txt")

# File can be found at GSE161065
bedfile <- fread("../../dat/mm10.repmask.s.bed")
bedfile[, repmask.name := paste0(V4, "_", V1, "-", V2, "-", V3)]
bedfile.merge <- merge(bedfile, repmask.rna.dict, by = "repmask.name")

# Bigwig files can be found at GSE175542
bwfiles <- data.table(path = dir("bigwig/", pattern = "*bw", full.names = T))
bwfiles <- bwfiles[grep("merge", path)]
bwfiles[, line := gsub("-.*", "", basename(path))]
bwfiles[, target := gsub("^.*-", "", gsub("_.*", "", basename(path)))]
bwfiles[, scale := gsub("^.*\\.", "", gsub("\\.bw", "", basename(path)))]
bwfiles$target <- factor(bwfiles$target, c("KDM5B", "SETDB1", "H3K9me3", "H3K4me3", "H3K4me2", "ATAC", "Input", "RNA"))
bwfiles <- bwfiles[order(line, decreasing = T)]
bwfiles <- bwfiles[order(target, decreasing = F)]

repelement <- "MMVL30-int"
repelement <- "RLTR6_Mm"
bedfile.merge.temp <- bedfile.merge[grepl(repelement, repmask.name) & rna.padj.multi < 0.05]
condition.message <- "rna.padj.multi < 0.05"
bedfile.merge.temp <- bedfile.merge[grepl(repelement, repmask.name) & rna.padj.unique < 0.05]
condition.message <- "rna.padj.unique < 0.05"
bedfile.merge.temp <- bedfile.merge[grepl(repelement, repmask.name) & rna.l2fc.multi > 0]
condition.message <- "rna.l2fc.multi > 0"
bedfile.merge.temp <- bedfile.merge[grepl(repelement, repmask.name) & rna.l2fc.unique > 0]
condition.message <- "rna.l2fc.unique > 0"
bedfile.merge.temp <- bedfile.merge[grepl(repelement, repmask.name)]
condition.message <- "all loci"

mmvl30.align <- fread("../../dat/mmvl30.alignment.final.bed")
bedfile.merge.temp <- mmvl30.align
condition.message <- "MMVL30 aligned"

mmvl30.stringtie.bed <- fread("test.bed")[V10 != "."]
unique(mmvl30.stringtie.bed$V10)
median(mmvl30.stringtie.bed$V9 - mmvl30.stringtie.bed$V8)
bedfile.merge.temp <- mmvl30.stringtie.bed[!duplicated(V10),.(V1 = V7, V2 = V8, V3 = V9, V4 = V10, V5 = V11, V6 = V6)]
condition.message <- "MMVL30 all stringtie"

condition.message <- paste0(condition.message, " n = ", nrow(bedfile.merge.temp))
bedfile.merge.temp.file <- "repmask.temporary.repelement.bed"
fwrite(bedfile.merge.temp[,.(V1, V2, V3, V4, V5, V6)], bedfile.merge.temp.file, sep = "\t", col.names = F)
out.file <- "tmp.txt"

# Get median size
median(bedfile.merge.temp$V3 - bedfile.merge.temp$V2)
# RLTR6-Mm: median 590, 2000:600:2000
# MMVL30-int: median 1921, 5000:2000:5000
# MMVL30-int: stringtie 5000, 5000:5000:5000
parameter <- "5000:5000:5000"
metaGene("KDM5B", "RPKM", "Input", "RPKM", repelement, condition.message, parameter)
metaGene("KDM5B", "scale", "Input", "RPKM", repelement, condition.message, parameter)
metaGene("SETDB1", "RPKM", "Input", "RPKM", repelement, condition.message, parameter)
metaGene("SETDB1", "scale", "Input", "RPKM", repelement, condition.message, parameter)
metaGene("H3K9me3", "RPKM", "Input", "RPKM", repelement, condition.message, parameter)
metaGene("H3K9me3", "scale", "Input", "RPKM", repelement, condition.message, parameter)
metaGene("H3K4me3", "RPKM", "Input", "RPKM", repelement, condition.message, parameter)
metaGene("H3K4me3", "scale", "Input", "RPKM", repelement, condition.message, parameter)
metaGene("H3K4me2", "RPKM", "Input", "RPKM", repelement, condition.message, parameter)
metaGene("H3K4me2", "scale", "Input", "RPKM", repelement, condition.message, parameter)
metaGene("RNA", "scale", repelement = repelement, condition.message = condition.message, parameter = parameter)

metaGene("KDM5B", "RPKM", repelement = repelement, condition.message = condition.message, parameter = parameter)
metaGene("KDM5B", "scale", repelement = repelement, condition.message = condition.message, parameter = parameter)
metaGene("SETDB1", "RPKM", repelement = repelement, condition.message = condition.message, parameter = parameter)
metaGene("SETDB1", "scale", repelement = repelement, condition.message = condition.message, parameter = parameter)
metaGene("H3K9me3", "RPKM", repelement = repelement, condition.message = condition.message, parameter = parameter)
metaGene("H3K9me3", "scale", repelement = repelement, condition.message = condition.message, parameter = parameter)
metaGene("H3K4me3", "RPKM", repelement = repelement, condition.message = condition.message, parameter = parameter)
metaGene("H3K4me3", "scale", repelement = repelement, condition.message = condition.message, parameter = parameter)
metaGene("H3K4me2", "RPKM", repelement = repelement, condition.message = condition.message, parameter = parameter)
metaGene("H3K4me2", "scale", repelement = repelement, condition.message = condition.message, parameter = parameter)
metaGene("RNA", "scale", repelement = repelement, condition.message = condition.message, parameter = parameter)

target1 <- "H3K9me3"
scale1 <- "scale"
target2 <- NULL
scale2 <- NULL
metaGene <- function(target1, scale1, target2 = NULL, scale2 = NULL, repelement, condition.message, parameter = "5000:5000:5000"){
  upstream <- as.numeric(strsplit(parameter, ":")[[1]][1])
  meta <- as.numeric(strsplit(parameter, ":")[[1]][2])
  downstream <- as.numeric(strsplit(parameter, ":")[[1]][3])
  
  if(is.null(target2)){
    tracks <- c(bwfiles[target %in% c(target1) & scale == scale1, get("path")])
  } else {
    tracks <- c(bwfiles[target %in% c(target1) & scale == scale1, get("path")],
                bwfiles[target %in% c(target2) & scale == scale2, get("path")])
  }
  tracks.collapse <- paste(tracks, collapse = ",")
  tracks.names <- gsub("\\.msd\\.redup\\.bam\\.center|\\.bw|_merge", "", basename(tracks))
  tracks.collapse
  
  system(paste0("source ~/.bashrc; bwtool aggregate ", parameter," ", bedfile.merge.temp.file, " ", tracks.collapse, " ", out.file, " -meta-scale-all;",
                " grep -v '#' ", out.file, " > test.txt;"))
  test <- fread("test.txt")
  colnames(test) <- c("Position", tracks.names)
  test.melt <- melt(test, id.vars = "Position", measure.vars = tracks.names, variable.name = "Bigwig")
  
  mp <- ggplot(data = test.melt, aes(Position, value, color = Bigwig)) +
    geom_line(size = 0.5) +
    theme_bw() +
    ylab(label = "ChIP-seq Signal") +
    scale_x_continuous(expand=c(0,0), breaks = c(-upstream, 0, meta, meta+downstream),
                       labels = c(upstream, "Start", "End", downstream)) +
    ggtitle(paste0(repelement, " ", condition.message)) +
    theme(axis.line = element_line(colour = "black", size = 0.25),
          axis.text.x = element_text(size = 12),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())
  mp
  ggsave(paste(c(repelement, "metagene", target1, scale1, target2, scale2, "pdf"), collapse = "."), mp, width = 6, height = 4)
  system("source ~/.bashrc; rm test.txt; rm tmp.txt")
}
