#!/usr/bin/env Rscript

# Relationship between input parameters and the ones used here
inputargs = commandArgs(trailingOnly=TRUE)
outputdir = inputargs[1]
inputed_ctrls = inputargs[2]
grnainfofiles = inputargs[3:length(inputargs)]

# Open pdf
pdf(paste(outputdir, "outputs/Log_counts.pdf", sep = ""))
par(mfrow=c(2,2))


#############################
## Plot useful information ##
#############################

for (gfile in grnainfofiles[1:(length(grnainfofiles))]) {

  # Check the test number
  savename = gsub(".*results_MAGeCK_", "", gfile)
  savename = gsub(".sgrna_summary.txt", "", savename)

  # Import the gRNA info files
  ginfo = read.table(gfile, header = T, stringsAsFactors = F)

  # Import hits
  hitspath = gsub("results_MAGeCK_.*", "", gfile)
  hits = read.table(
    paste(hitspath, "hits.all_MAGeCK_", savename, ".txt", sep = ""),
    header = T, stringsAsFactors = F)
  hitstwo = read.table(
    paste(hitspath, "hits.all_PBNPA_", savename, ".txt", sep = ""),
    header = T, stringsAsFactors = F)

  # Import the normalized count files
  countspath = gsub("sgrna_summary", "normalized", gfile)
  counts = read.table(countspath, header = T, stringsAsFactors = F)
  counts_second = counts


  ## CREATE 1st PLOT ##
  plot(ginfo$LFC, -log10(ginfo$FDR),
       main = paste("Volcano plot, test", savename, sep = ""),
       pch = 19, col = "gray",
       xlab = "Log2 fold-change",
       ylab = "-Log10(adjusted p-value)")
  abline(v=0, col = "dimgray", lty = 5)
  used_colors =  rainbow(nrow(hits)+nrow(hitstwo))
  for (j in 1:nrow(hits)){
    points(ginfo$LFC[ginfo$Gene %in% hits$group_id[j]],
           -log10(ginfo$FDR)[ginfo$Gene %in% hits$group_id[j]],
           pch = 19, col = used_colors[j], cex = 0.7)
  }
  for (k in 1:nrow(hitstwo)){
    points(ginfo$LFC[ginfo$Gene %in% hitstwo$Gene[k]],
           -log10(ginfo$FDR)[ginfo$Gene %in% hitstwo$Gene[k]],
           pch = 19, col = used_colors[k+j], cex = 0.7)
  }
  # Add legend
  legend(x = "topleft", bty = "n", pch = 19,
         legend = c("gRNAs of selected genes"),
         col = c("red"), x.intersp= 0.7, cex = 0.5)


  ## CREATE 2nd PLOT ##
  # Set graph limits
  minval = min(ginfo$LFC)
  maxval = max(ginfo$LFC)

  # Sort the dataframe by the log2 fold change
  ginfoordered = ginfo[as.numeric(order(ginfo$LFC)),]

  # Plot figures
  # Open an empty plot that allows to add lines afterwards
  plot(1,minval, type = "n",
       xlim = c(1, nrow(ginfoordered)),
       ylim = c(minval, maxval),
       xlab = "",
       axes = F,
       ylab = "Log2 fold-change",
       main = paste("Log2 fold-change, test ", savename, sep =""))
  # Plot a line following all of the gRNAs
  lines(seq(1,nrow(ginfoordered)), ginfoordered$LFC)
  # Add 0 y-axis
  abline(0,0, col = "dimgray", lty = 5)
  # Write y-axis
  axis(side=2,
       at = seq(minval, maxval, length.out = 5),
       labels = round(seq(minval, maxval, length.out = 5), 1))
  # Write x-label
  mtext(side=1, line=0, "Ranked gRNAs")
  # Add position of guides targeting hits to the plot
  for (j in 1:nrow(hits)){
    points(seq(1,nrow(ginfoordered))[ginfoordered$Gene %in% hits$group_id[j]],
           ginfoordered$LFC[ginfoordered$Gene %in% hits$group_id[j]],
           pch = 19, col = used_colors[j], cex = 0.7)
  }
  for (k in 1:nrow(hitstwo)){
    points(seq(1,nrow(ginfoordered))[ginfoordered$Gene %in% hitstwo$Gene[k]],
           ginfoordered$LFC[ginfoordered$Gene %in% hitstwo$Gene[k]],
           pch = 19, col = used_colors[j+k], cex = 0.7)
  }
  # Add legend
  legend(x = "topleft", bty = "n", pch = 19,
         legend = c("gRNAs of selected genes"),
         col = c("red"), x.intersp= 0.7, cex = 0.5)


  ## CREATE 3rd PLOT ##
  # Replace 0s by 1E-10 from counts file to avoid problems with log2
  for (jcol in 3:ncol(counts)) {
    counts[counts[,jcol] == 0,jcol] = 1
  }

  # Determine how many samples appear in the counts table
  numcols = ncol(counts)
  numsamples = numcols - 2

  # Determine how many samples are treatments or controls
  numtreated = numsamples/2
  dayzero = 3:(numtreated+2)
  daytreated = (numtreated+3):numcols

  # Compute the log2
  counts[,3:numcols] = log2(counts[,3:numcols])

  # Separate controls from treatments
  if (length(dayzero) > 1) {
    logcontrol = rowSums(counts[,dayzero])
    logtreated = rowSums(counts[,daytreated])
  } else {
    logcontrol = counts[,dayzero]
    logtreated = counts[,daytreated]
  }
  names(logcontrol) = counts$Gene
  names(logtreated) = counts$Gene

  # Plot log2 counts
  plot(logcontrol, logtreated,
       main = paste("Log2 counts, test", savename, sep = ""),
       pch = 19, col = "dimgray",
       xlab = "Log2 counts, untreated",
       ylab = "Log2 counts, treated")
  for (j in 1:nrow(hits)){
    points(logcontrol[names(logcontrol) %in% hits$group_id[j]],
           logtreated[names(logtreated) %in% hits$group_id[j]],
           pch = 19, col = used_colors[j], cex = 0.7)
  }
  for (k in 1:nrow(hitstwo)){
    points(logcontrol[names(logcontrol) %in% hitstwo$Gene[k]],
           logtreated[names(logtreated) %in% hitstwo$Gene[k]],
           pch = 19, col = used_colors[k+j], cex = 0.7)
  }
  # Add legend
  legend(x = "topleft", bty = "n", pch = 19,
         legend = c("gRNAs of selected genes"),
         col = c("red"), x.intersp= 0.7, cex = 0.5)


  ## CREATE 4th PLOT ##
  # Define how many colors will be used
  used_colors =  rainbow(ncol(counts_second)-2)

  # Open an empty plot that allows to add lines afterwards
  plot(0,0, type = "n",
       xlim = c(1, nrow(counts_second)),
       ylim = c(0, 100),
       xlab = "gRNAs",
       ylab = "Percentage of total counts (%)",
       main = paste("Cumulative counts, test ", savename, sep = ""))

  # Add lines for all the different samples
  used_legend = numeric()
  for (jcol in 3:ncol(counts_second)) {
    # Compute and plot cumulative distribution
    cumul_ctr = cumsum(sort(counts_second[,jcol], decreasing = T))
    cumul_dist = (cumul_ctr/sum(counts_second[,jcol]))*100
    lines(seq(1,length(cumul_dist), by = 1), cumul_dist,
          col = used_colors[jcol-2])
    used_legend[jcol-2] = colnames(counts_second)[jcol]
  }
  legend(x = "topleft", bty = "n", pch = 19,
         legend = used_legend,
         col = used_colors,
         x.intersp=0.5,
         inset = c(0, -0.1), xpd = T,
         horiz = T, cex = 0.5)
}

message = dev.off()

##########
## DONE ##
##########
