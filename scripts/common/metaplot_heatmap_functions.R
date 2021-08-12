
# Functions
getMeltedSizeDistribution <- function(bedList, filter){
  sizeTable <- data.frame(size = 1:filter)
  for(i in names(bedList)){
    size <- getSizes(bedList[[i]], filter+1)
    tmp <- data.frame(table(size))
    colnames(tmp) <- c("size", i)
    tmp[,1] <- as.numeric(as.character(tmp[,1]))
    sizeTable <- merge(sizeTable, tmp, by = "size", all = TRUE)
  }
  #colnames(sizeTable)[2:ncol(sizeTable)] <- names(bedList)
  melt(sizeTable, id.vars = "size")
}
getMeltedDistanceDistribution <- function(bedList, filter){
  disTable <- data.frame(distance = 1:filter)
  for(i in names(bedList)){
    dis <- getDistances(bedList[[i]], filter+1)
    tmp <- data.frame(table(dis))
    colnames(tmp) <- c("distance", i)
    tmp[,1] <- as.numeric(as.character(tmp[,1]))
    disTable <- merge(disTable, tmp, by = "distance", all = TRUE)
  }
  colnames(disTable)[2:ncol(disTable)] <- names(bedList)
  melt(disTable, id.vars = "distance")
}

# plotProfile <- function(bed, tracks, names, cPalette, returnDataFrame = FALSE, centered = FALSE, dataframeOutput, interval = c(-1000:-1,1:1000)){
plotProfile <- function(bed, tracks, names, returnDataFrame = FALSE, centered = FALSE, dataframeOutput, interval = c(-1000:-1,1:1000)){
    if(nrow(bed) > 100000){
    warning("Large bedfile!")
  } else {
    signals <- list()
    write.table(bed, "temp.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    signals[["position"]] <- data.frame(position = interval)
    for(i in tracks){
      stdout <- system(paste0("source ~/.bashrc; bwtool aggregate ", abs(min(interval)), ":", abs(max(interval)), " temp.bed ", i, " /dev/stdout"), intern = TRUE)
      con <- textConnection(paste(stdout, collapse = '\n'))
      signals[[i]] <- read.table(con , header = FALSE, stringsAsFactors = FALSE, sep = "\t")[,2,drop = FALSE]
      colnames(signals[[i]]) <- i
      close(con)
    }
    res <- do.call(cbind, signals)
    if(centered){
      modified <- scale(res[,2:ncol(res)], scale = FALSE)
      ylab <- "Mean Centered Signal"
    } else {
      modified <- res[,2:ncol(res)]
      ylab <- "Mean Signal"
    }
    #modified <- scale(res[,2:ncol(res)], center = FALSE, scale = FALSE)
    #modified <- scale(res[,2:ncol(res)], scale = FALSE)
    colnames(modified) <- names
    
    if(returnDataFrame){
      return(melt(cbind(res[,1,drop = FALSE], modified), id.vars = "position"))
    }
    fwrite(cbind(res[,1,drop = FALSE], modified), dataframeOutput, sep = "\t", row.names = FALSE)
    myplot <- ggplot2::ggplot(melt(cbind(res[,1,drop = FALSE], modified), id.vars = "position"), ggplot2::aes(position, value, color = variable)) + ggplot2::geom_line() +
      ggplot2::theme_bw() + 
      # ggplot2::scale_color_manual(values = cPalette) +
      ggplot2::xlab(label = "Relative Position (bp)") +
      ggplot2::ylab(label = ylab) +
      ggplot2::geom_vline(xintercept=0, size = 0.5, linetype = 3) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(size = 14, hjust = 0.5), 
                     axis.title = ggplot2::element_text(size = 12), 
                     axis.text = ggplot2::element_text(size = 10))
    #myplot
    
    # Clean
    system("source ~/.bashrc; rm temp.bed;")
    
    return(myplot)
    

  }
  

}

getMatrix <- function(bed, tracks, window = c(-2000,2000), colN = 100, sums = FALSE){
  # message(nrow(bed))
  bed <- as.data.table(bed)
  if(ncol(bed) < 4){
    colnames(bed) <- c("V1", "V2", "V3")
    bed[, V4 := paste(V1, V2, V3, sep = "-")]
  } else {
    colnames(bed) <- c("V1", "V2", "V3", "V4")
  }
  # message(nrow(bed))
  # 
  require(data.table)
  
  if(nrow(bed) > 1000000){
    warning("Large bedfile!")
  } else {
    #colN <- 100 #Has to be divible into window!!!
    tiled <- (abs(min(window)) + abs(max(window)))/colN
    fileConn<-file("bigwiglist.tmp.lst")
    writeLines(tracks, fileConn)
    close(fileConn)

    write.table(as.data.frame(bed)[,1:4], "temp.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    
    system(paste0("source ~/.bashrc; bwtool matrix ", abs(min(window)), ":", abs(max(window)), " -tiled-averages=", tiled, " temp.bed bigwiglist.tmp.lst heatmap.output.tmp.txt;"))
    df <- fread("heatmap.output.tmp.txt")
    if(sums){
      df.sums <- data.frame(matrix(nrow = nrow(df)))
      df.means <- data.frame(matrix(nrow = nrow(df)))
      for(i in 1:length(tracks)){
        df.sums[,i] <- rowSums(df[,.SD, .SDcols = c(((i-1)*colN+1):(i*colN))], na.rm = TRUE)
        df.means[,i] <- rowMeans(df[,.SD, .SDcols = c(((i-1)*colN+1):(i*colN))], na.rm = TRUE)
      }
      colnames(df.sums) <- paste0("sum.", 1:length(tracks))
      colnames(df.means) <- paste0("mean.", 1:length(tracks))
      
      df <- as.data.table(cbind(df.sums, df.means, as.data.frame(bed)[,1:4], ".", ".", df))
      # Can't sort by sums because that messes up order
      # for(col in rev(colnames(df.sums))){
      #   df <- df[order(get(col), decreasing = TRUE),]
      # }
    } else {
      df <- cbind(as.data.frame(bed)[,1:4], ".", ".", df)
    }
  
    # Clean
    system("source ~/.bashrc; rm bigwiglist.tmp.lst;")
    
    df <- as.data.table(df)
    hm.attr <- list()
    hm.attr[["sampleN"]] <- length(tracks)
    hm.attr[["sampleN"]] <- length(tracks)
    hm.attr[["colN"]] <- colN
    hm.attr[["tiled"]] <- tiled
    hm.attr[["upstream"]] <- abs(min(window))
    hm.attr[["downstream"]] <- abs(max(window))
    hm.attr[["group_boundaries"]] <- paste0(0, ",", nrow(df))
    hm.attr[["sample_boundaries"]] <- paste0(0:length(tracks)*colN, collapse = ",")
    
    attr(df, which = "heatmapAttr") <- hm.attr
    return(df)

  }
}

sortBySums <- function(df, sumsCols){
  df <- as.data.table(df)
  for(col in sumsCols){
    df <- df[order(get(col), decreasing = TRUE),]
  }
  return(df)
}

sortByMeans <- function(df, meanCols){
  df <- as.data.table(df)
  for(col in meanCols){
    df <- df[order(get(col), decreasing = TRUE),]
  }
  return(df)
}

#sampleNames <- c("GATA6", "ATAC", "H3K27ac", "H3K4me3", "Input")
#outFile <- "test.txt"
writeHeatmapMatrix <- function(df, outFile, heatmapAttr = NULL, sampleNames){
  require(data.table)
  if(is.null(heatmapAttr)){
    heatmapAttr <- attributes(df)$heatmapAttr 
  }
  if(length(sampleNames) != heatmapAttr[["sampleN"]]){
    message("Error, sample number doesn't match names.")
  } else {
    sampleN <- heatmapAttr[["sampleN"]]
    colN <- heatmapAttr[["colN"]]
    tiled <- heatmapAttr[["tiled"]]
    upstream <- heatmapAttr[["upstream"]]
    downstream <- heatmapAttr[["downstream"]]
    sample_boundaries <- heatmapAttr[["sample_boundaries"]]
    group_boundaries <- heatmapAttr[["group_boundaries"]]
    sample_labels <- gsub(" ", "", toString(shQuote(sampleNames)))

    headerLine <- gsub("'", "\\\"", paste0("@{",shQuote("verbose"),":true,",
                                           shQuote("scale"),":1,",shQuote("skip zeros"),
                                           ":false,",shQuote("nan after end"),":false,",
                                           shQuote("sort using"),":",shQuote("sum"),",",
                                           shQuote("unscaled 5 prime"),":0,",shQuote("body"),
                                           ":0,",shQuote("sample_labels"),":[",sample_labels,"],",
                                           shQuote("downstream"),":",downstream,",",shQuote("unscaled 3 prime"),
                                           ":0,",shQuote("group_labels"),":[",shQuote("peaks"),"],",
                                           shQuote("bin size"),":",tiled,",",shQuote("upstream"),":",upstream,",",
                                           shQuote("group_boundaries"),":[", group_boundaries,"],",
                                           shQuote("sample_boundaries"),":[", sample_boundaries,"],",
                                           shQuote("missing data as zero"),":true,",shQuote("ref point"),
                                           ":",shQuote("center"),",",shQuote("min threshold"),":null,",
                                           shQuote("sort regions"),":",shQuote("no"),",",
                                           shQuote("proc number"),":20,",shQuote("bin avg type"),":",
                                           shQuote("mean"),",",shQuote("max threshold"),":null}\n"))
    cat(headerLine, file = outFile)
    #fwrite(as.data.table(headerLine), outFile, append = TRUE, col.names = FALSE)
    fwrite(df, outFile, append = TRUE, col.names = FALSE, sep = "\t", quote = FALSE, na = 0)
    system(paste0("source ~/.bashrc; gzip < ", outFile, " > ", paste0(outFile, ".gz")))
  }
}

getDistances <- function(df, filter){
  distances <- cbind(df[1:(nrow(df)-1), 2:3], df[2:nrow(df), 2:3])
  #colnames(distances) <- c("start1", "end1", "start2", "end2")
  #distances$dis <- distances$start2 - distances$end1
  dis <- distances[,3] - distances[,2]
  dis[dis < filter]
  #nums <- distances$dis[distances$dis < filter]
  #nums
}

getSizes <- function(df, filter){
  sizes <- df[,3]-df[,2]
  sizes[sizes < filter]
}

# bedMerge <- function(df, size){
#   df <- data.frame(df)
#   output <- matrix(NA, nrow = nrow(df), ncol = 3)
#   entry <- c(df[1,1], df[1,2], df[1,3])
#   for(i in 1:(nrow(df))){
#     if(i == nrow(df)){
#       break
#     }
#     if(df[i+1,2] - entry[3] <= size){
#       entry[3] <- df[i+1,3]
#     } else {
#       output[i,] <- entry
#       entry <- c(df[i+1,1], df[i+1,2], df[i+1,3])
#     }
#   }
#   output[i,] <- entry
#   output[complete.cases(output),]
# }

bedFilter <- function(df, size){
  df <- data.frame(df)
  df[which(df[,3] - df[,2] > size),]
}

bedMerge2 <- function(df, size){
  write.table(df[,1:3], "temp.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  stdout <- system(paste0("bedtools merge -d ", size, " -i temp.bed"), intern = TRUE)
  con <- textConnection(paste(stdout, collapse = '\n'))
  output <- read.table(con , header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  close(con)
  return(output)
}

bedIntersect <- function(df1, df2){
  write.table(df1[,1:3], "temp1.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(df2[,1:3], "temp2.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  stdout <- system(paste0("bedtools intersect -a temp1.bed -b temp2.bed"), intern = TRUE)
  con <- textConnection(paste(stdout, collapse = '\n'))
  output <- read.table(con , header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  close(con)
  return(output)
}

bedElementOf <- function(a, b, length){
  write.table(a[,1:3], "temp1.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(b[,1:3], "temp2.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  stdout <- system(paste0("bedops --element-of ", length, " temp1.bed temp2.bed"), intern = TRUE)
  con <- textConnection(paste(stdout, collapse = '\n'))
  output <- read.table(con , header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  close(con)
  return(output)
}
