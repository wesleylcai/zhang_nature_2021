# Last Updated: 
# Author: Wesley Cai
# Purpose:

library(data.table)
library(ggplot2)
library(reshape2)
library(janitor)

# OG: Google_Drive/medschool/research/qin/chip/analysis/210520_signal_correlation/210520_signal_correlation_updated.R

mainwd <- "."
outputfolder <- "output/"
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

counts <- fread("../../dat/YR-KDM5B.counts.txt")
colnames(counts) <- gsub("^([0-9])", "ko\\1", gsub("-", ".", gsub("\\.msd.*", "", basename(colnames(counts)))))
counts.summary <- fread("../../dat/YR-KDM5B.counts.txt.summary")
colnames(counts.summary) <- gsub("^([0-9])", "ko\\1", gsub("-", ".", gsub("\\.msd.*", "", basename(colnames(counts.summary)))))

#### Get RPKM ####
# counts.summary.total <- adorn_totals(counts.summary, where = "row")
# fwrite(counts.summary.total, "totals.test.txt", sep = "\t")
counts.summary <- counts.summary[,.SD,.SDcols = 2:ncol(counts.summary)]
samples <- unique(gsub("_R.*", "", colnames(counts.summary)))
for(sam in samples){
  if(!grepl("1G8|815|old", sam)){
    counts.summary[, eval(paste0(sam, "_merge")) := get(paste0(sam, "_Rep1")) + get(paste0(sam, "_Rep2"))]
  }
}
counts.summary.total <- colSums(counts.summary)/1e6

counts.ma <- as.matrix(counts[,.SD,.SDcols = 7:ncol(counts)])
# counts.ma.rpm <- sweep(counts.ma, 2, counts.summary.total, "/")
# counts.ma.rpkm <- sweep(counts.ma.rpm, 1, counts$Length, "/")
# chip.rpkm <- as.data.table(cbind(counts[,.(Geneid)], counts.ma.rpkm))


## Normal RPKM
counts.ma.plus1 <- counts.ma + 1
counts.ma.rpm.plus1 <- sweep(counts.ma.plus1, 2, counts.summary.total[1:ncol(counts.ma.plus1)], "/")
counts.ma.rpkm.plus1 <- sweep(counts.ma.rpm.plus1, 1, counts$Length, "/")
chip.log2.rpkm.plus1 <- as.data.table(cbind(counts[,.(Geneid)], log2(counts.ma.rpkm.plus1)))
fwrite(chip.log2.rpkm.plus1, "chip.log2rpkmplus1.replicates.txt", sep = "\t")
## Normal RPKM


## Merged/Input RPKM ##
counts.ma.dt <- as.data.table(counts.ma)
# Only keep 4H3 and YR
counts.ma.dt <- counts.ma.dt[,.SD,.SDcols = grepl("4H3|YR", colnames(counts.ma.dt)) &
                                           !(grepl("old", colnames(counts.ma.dt)))]
# Combine replicates
samples <- unique(gsub("_R.*", "", colnames(counts.ma.dt)))
for(sam in samples){
  if(!grepl("1G8|815|old", sam)){
    counts.ma.dt[, eval(paste0(sam, "_merge")) := get(paste0(sam, "_Rep1")) + get(paste0(sam, "_Rep2"))]
  }
}
# Combine replicates for total number of reads
counts.ma.dt <- counts.ma.dt[,.SD,.SDcols = grepl("merge", colnames(counts.ma.dt))]
counts.total.dt <- data.table(sample = names(counts.summary.total), total = counts.summary.total)
counts.total.df <- as.data.frame(counts.total.dt)
row.names(counts.total.df) <- counts.total.dt[, get("sample")]
counts.total.df.subset <- counts.total.df[colnames(counts.ma.dt),]
# Double-check row names of totals match colnames of read counts
row.names(counts.total.df.subset)
colnames(counts.ma.dt)
# Get RPKM
counts.ma <- as.matrix(counts.ma.dt)
counts.ma.plus1 <- counts.ma + 1
counts.ma.rpm.plus1 <- sweep(counts.ma.plus1, 2, counts.total.df.subset$total, "/")
counts.ma.rpkm.plus1 <- sweep(counts.ma.rpm.plus1, 1, counts$Length, "/")
class(counts.ma.rpkm.plus1)
# Get ratios of Target RPKM/Input RPKM
counts.ma.rpkm.plus1.dt <- as.data.table(counts.ma.rpkm.plus1)
samples <- unique(gsub("_merge.*", "", colnames(counts.ma.rpkm.plus1.dt)))
for(sam in samples){
  line <- gsub("\\..*", "", sam)
  counts.ma.rpkm.plus1.dt[, eval(paste(sam, line, "Input", sep = ".")) := get(paste0(sam, "_merge"))/get(paste0(line, ".Input_merge"))]
}
# Only keep ratios and then log2
counts.ma.rpkm.plus1.dt <- counts.ma.rpkm.plus1.dt[,.SD,.SDcols =
                                                     grepl("\\..*\\..*\\.Input", colnames(counts.ma.rpkm.plus1.dt)) &
                                                     !(grepl("Input.*Input", colnames(counts.ma.rpkm.plus1.dt)))]

chip.log2.rpkm.plus1.dt <- as.data.table(cbind(counts[,.(Geneid)], log2(counts.ma.rpkm.plus1.dt)))
fwrite(chip.log2.rpkm.plus1.dt, "chip.log2ratiostargetinput.merge.txt", sep = "\t")
chip.log2.rpkm.plus1.dt <- fread("chip.log2ratiostargetinput.merge.txt")
## Merged/Input RPKM ##

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#### Plot Merged/Input RPKM #####
tmp <- chip.log2.rpkm.plus1.dt
tmp.combn <- as.data.table(t(combn(colnames(tmp[,.SD,.SDcols = 2:ncol(tmp)]), 2)))
colnames(tmp.combn) <- c("comp1", "comp2")
comp1 <- "YR.KDM5B.YR.Input"
comp2 <- "YR.SETDB1.YR.Input"
comp2 <- "YR.H3K9me3.YR.Input"
rown <- 1

# Only do a subset
tmp.combn <- data.table(comp1 = c("YR.KDM5B.YR.Input", "YR.KDM5B.YR.Input"), comp2 = c("YR.SETDB1.YR.Input", "YR.H3K9me3.YR.Input"))

library(showtext)
font_add("Arial", "/Applications/Microsoft Word.app/Contents/Resources/DFonts/arial.ttf")
showtext_auto()

## MODIFIED TO ONLY KEEP POSITIVE PEAKS (variable numbers)
for(rown in 1:nrow(tmp.combn)){
  comp1 <- tmp.combn[rown, get("comp1")]
  comp2 <- tmp.combn[rown, get("comp2")]
  tmp.positives <- tmp[get(comp1) > 0 & get(comp2) > 0]
  cor.res <- cor.test(tmp.positives[, get(comp1)], tmp.positives[, get(comp2)], method = "pearson", alternative = "two.sided")
  cor.res$p.value
  cor.res$estimate
  tmp.combn[rown, estimate := cor.res$estimate]
  cor.res.p.value <- cor.res$p.value
  if(cor.res.p.value == 0) cor.res.p.value <- "p < 2.2e-16"
  tmp.combn[rown, p.value := cor.res.p.value]
  # https://stackoverflow.com/questions/36669095/y-axis-wont-start-at-0-in-ggplot
  tmp.positives[, dens := get_density(tmp.positives[, get(comp1)], tmp.positives[, get(comp2)], h = c(0.5,0.5), n = 200)]
  nrow(tmp.positives)
  mp <- ggplot(tmp.positives, aes_string(comp1, comp2, color = "dens")) +
    geom_point(size = 0.1) +
    scale_color_viridis_c() +
    geom_smooth(formula = "y ~ x", method = "lm", se = F, color = "red") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 5), breaks = 1:4) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 5), breaks = 1:4) +
    theme_bw() +
    theme(axis.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 18, color = "black"),
          panel.border = element_rect(fill=NA, colour = "black", size=1.4),
          axis.ticks = element_line(colour = "black", size = 1),
          legend.position = "none",
          text=element_text(family="Arial"))
  mp
  ggsave(paste0(comp1, "_v_", comp2, ".positives.png"), mp, width = 4, height = 4, dpi = 1200)
}
fwrite(tmp.combn[comp2 != "dens"], "logratiorpkm.cor.pvalues.txt", sep = "\t")
#### Plot Merged/Input RPKM #####

#### Plot replicate RPKM ####
chip.log2.rpkm.plus1[, input.merge := YR.Input_Rep1 + YR.Input_Rep2 + ko4H3.Input_Rep1 + ko4H3.Input_Rep2]
chip.log2.rpkm.plus1 <- chip.log2.rpkm.plus1[order(input.merge, decreasing = T)]
input.cutoff <- quantile(chip.log2.rpkm.plus1$input.merge, probs = c(0.99))
tmp <- chip.log2.rpkm.plus1[input.merge < input.cutoff,.SD,.SDcols = grepl("_Rep", colnames(chip.log2.rpkm.plus1))]
tmp.sampled <- tmp[sample(1:nrow(tmp), 5e3)]

tmp.combn <- as.data.table(t(combn(colnames(tmp[,.SD,.SDcols = 2:ncol(tmp)]), 2)))
colnames(tmp.combn) <- c("comp1", "comp2")
for(rown in 1:nrow(tmp.combn)){
  comp1 <- tmp.combn[rown, get("comp1")]
  comp2 <- tmp.combn[rown, get("comp2")]
  cor.res <- cor.test(tmp[, get(comp1)], tmp[, get(comp2)], method = "pearson", alternative = "two.sided")
  tmp.combn[rown, estimate := cor.res$estimate]
  cor.res.p.value <- cor.res$p.value
  # message(cor.res.p.value == 0)
  # if(cor.res.p.value == 0){cor.res.p.value <- "p < 2.2e-16"}
  tmp.combn[rown, p.value := cor.res.p.value]
  
  if(grepl("ko4H3|YR", comp1) & grepl("ko4H3|YR", comp2)){
    mp <- ggplot(tmp.sampled, aes_string(comp1, comp2)) +
      geom_point(size = 0.1) +
      geom_smooth(formula = "y ~ x", method = "lm", se = F, color = "red") +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 18, color = "black"),
            panel.border = element_rect(fill=NA, colour = "black", size=2),
            axis.ticks = element_line(colour = "black", size = 1),
            legend.position = "none")
    ggsave(paste0(comp1, "_v_", comp2, ".5Ksample.pdf"), mp, width = 5, height = 5)
  }
}
fwrite(tmp.combn[comp2 != "dens"], "logrpkmreplicates.cor.pvalues.txt", sep = "\t")
#### Plot replicate RPKM ####
