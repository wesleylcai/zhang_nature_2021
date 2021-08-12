# Last Updated: 
# Author: Wesley Cai
# Purpose:

library(data.table)
library(ggplot2)
library(ggvenn)

# OG: Google_Drive/medschool/research/qin/chip/analysis/210516_overlaps/210518_overlaps_H3K9me3_KDM5B_SETDB1_updated.R

mainwd <- "."
outputfolder <- "output/"
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

# File can be found at GSE161065: Superseries
repmask.rna.dict <- fread("../../dat/mm10.repmask.rna.dict.txt")

library(showtext)
font_add("Arial", "/Applications/Microsoft Word.app/Contents/Resources/DFonts/arial.ttf")
showtext_auto()

#### TE stats ####
multi.repmask.file <- "../../dat/YR-H3K9me3.YR-KDM5B.YR-SETDB1.common.repmask.bed"
multi.repmask <- fread(multi.repmask.file)
multi.repmask[, repmask.name := paste0(V8, "_", V5, "-", V6, "-", V7)]
multi.repmask.merge <- merge(multi.repmask, repmask.rna.dict, by = "repmask.name")
multi.repmask.merge.te <- multi.repmask.merge[final.type %in% c("DNA", "LTR", "SINE", "LINE")]
multi.repmask.merge.te[, final.type := factor(final.type, c("LTR", "LINE", "SINE", "DNA"))]

mp <- ggplot(multi.repmask.merge.te[rna.padj.unique < 0.05], aes(final.type, rna.l2fc.unique)) +
  geom_hline(yintercept = 0, linetype="dashed", color = "gray") +
  geom_boxplot(width = 0.5, aes(fill = final.type), show.legend = F) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  ylab(label = expression("RNA "*Log[2]*"(FC)")) +
  # ylim(c(-8, 15)) +
  scale_fill_discrete(name = "RE type") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1),
        text=element_text(family="Arial"))
ggsave("repmask.YR-H3K9me3.YR-KDM5B.YR-SETDB1.regions.png", mp, width = 3.5, height = 3.5, dpi = 1200)

p.table <- data.table(final.type = c("DNA", "LTR", "LINE", "SINE"))
multi.repmask.merge.te.sig <- multi.repmask.merge.te[rna.padj.unique < 0.05]
for(rown in 1:nrow(p.table)){
  final.type.var <- p.table[rown, get("final.type")]
  p.table[rown, sig.up := nrow(multi.repmask.merge.te.sig[final.type == final.type.var & rna.l2fc.unique > 0])]
  p.table[rown, sig.down := nrow(multi.repmask.merge.te.sig[final.type == final.type.var & rna.l2fc.unique < 0])]
  
  p.table[rown, LTR.mann.whitney.p.value := wilcox.test(multi.repmask.merge.te.sig[final.type == "LTR", get("rna.l2fc.unique")],
                                                        multi.repmask.merge.te.sig[final.type == final.type.var, get("rna.l2fc.unique")])$p.value]
  p.table[rown, LTR.t.test.p.value := t.test(multi.repmask.merge.te.sig[final.type == "LTR", get("rna.l2fc.unique")],
                                             multi.repmask.merge.te.sig[final.type == final.type.var, get("rna.l2fc.unique")])$p.value]
  
}
fwrite(p.table, "repmask.YR-H3K9me3.YR-KDM5B.YR-SETDB1.pvalue.table.txt", sep = "\t")
#### TE stats ####

#### Venn diagrams ####
yr.setdb1.file <- "../../dat/YR-SETDB1_merge.bed"
yr.kdm5b.file <- "../../dat/YR-KDM5B_merge.bed"
yr.h3k9me3.file <- "../../dat/YR-H3K9me3_merge.bed"
mult.int.bed.file <- "temp.mult.int.bed"

system(paste("source ~/.bashrc; cat", yr.h3k9me3.file, yr.kdm5b.file, yr.setdb1.file, "| bedtools sort -i - | bedtools merge -i - > temp.all.peaks.bed;", sep =" "))
system(paste("source ~/.bashrc; bedtools intersect -wa -wb -a temp.all.peaks.bed -b", yr.h3k9me3.file, yr.kdm5b.file, yr.setdb1.file, ">", mult.int.bed.file, sep =" "))

mult.int <- fread(mult.int.bed.file)

mult.int[,peakid := paste(V1, V2, V3, sep = "-")]
mult.int[, h3k9me3 := 0]
mult.int[V4 == 1, h3k9me3 := 1]
mult.int[, kdm5b := 0]
mult.int[V4 == 2, kdm5b := 1]
mult.int[, setdb1 := 0]
mult.int[V4 == 3, setdb1 := 1]

mult.int.unique <- mult.int[,.(chr = unique(V1), start = unique(V2), end = unique(V3),
                               H3K9me3 = sum(h3k9me3), KDM5B = sum(kdm5b), SETDB1 = sum(setdb1)), by = peakid]

venn.final <- list()
nm <- "KDM5B"
for(nm in c("KDM5B", "SETDB1", "H3K9me3")){
  venn.final[[nm]] <- mult.int.unique[get(nm) > 0, get("peakid")]
}
mp <- ggvenn(venn.final[c(1,2,3)], text_size = 3.5, fill_color = c("red", "dodgerblue2", "lawngreen"),
             stroke_size = 0.5, set_name_size = 4)
mp
ggsave("YR-KDM5B.YR-SETDB1.YR-H3K9me3.venn.pdf", mp, width = 5, height = 5)
mp <- ggvenn(venn.final[c(1,2,3)], text_size = 0, fill_color = c("red", "dodgerblue2", "lawngreen"),
             stroke_size = 0.5, set_name_size = 4)
mp
ggsave("YR-KDM5B.YR-SETDB1.YR-H3K9me3.venn.nolabel.pdf", mp, width = 5, height = 5)


png("KDM5B.H3K9me3.venn.scale.png", width = 3, height = 3, units = "in", res = 300)
venn.plot <- draw.pairwise.venn(
  area1 = 56.5+31.2,
  area2 = 12.3+31.2,
  cross.area = 31.2,
  category = c("KDM5B", "H3K9me3"),
  fill = c("red", "green"),
  alpha = c(0.6, 0.6),
  lty = 1,
  cex = 0,
  cat.cex = 1,
  cat.pos = c(210, 130),
  cat.dist = 0.09,
  cat.just = list(c(-1, -1), c(1, 1)),
  ext.pos = 30,
  ext.dist = -0.05,
  ext.length = 0.85,
  ext.line.lwd = 2,
  ext.line.lty = "dashed")
dev.off()
#### Venn diagrams ####