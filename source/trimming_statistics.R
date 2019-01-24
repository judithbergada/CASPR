#!/usr/bin/env Rscript

# Relationship between input parameters and the ones used here
files_total = commandArgs(trailingOnly=TRUE)

# Open pdf to save plots
outputdir = files_total[2]
name_pdf = paste(outputdir, "/outputs/Trimming_statistics.pdf", sep = "")
pdf(name_pdf)

#################################
## Plot statistics of trimming ##
#################################

# Determine how many nucleotides have been cut before trimming
i=1
cuted = numeric(1)
for (k in unlist(strsplit(files_total[1], split="\n"))) {
  cuted[i] = as.numeric(k)
  i = i + 1
}

# Prepare information of results
results = data.frame(numeric(2))
rownames(results) = c("Total reads", "Trimmed reads")

counter = 1
for (filename in files_total[3:length(files_total)]) {

  # Obtain SampleName
  name=gsub(".*trim_stat_", "", filename)
  name=gsub(".txt", "", name)

  # Read file
  info = readLines(filename)

  # Obtain number of input reads
  readsnumb = grep("(Total read pairs processed|Total reads processed)",
                    info, value = TRUE)
  assign(paste("readsnumb_", name, sep = ""),
         as.numeric(
           gsub(",", "",
                unlist(regmatches(
                  readsnumb, gregexpr(
                    "[[:digit:]]+(,?[[:digit:]]+)+", readsnumb))))))

  # Obtain percentage of trimmed reads
  trimmed = grep("with adapter", info, value = TRUE)
  assign(paste("trimmed_", name, sep = ""),
         max(as.numeric(
           gsub("%", "",
                unlist(regmatches(
                  trimmed, gregexpr(
                    "[[:digit:]]+(\\.[[:digit:]]+)?%", trimmed)))))))

  # Create dataframe with the results
  results[,counter] = c(get(paste("readsnumb_",name, sep = "")),
                        get(paste("trimmed_",name, sep = "")))
  colnames(results)[counter] = name

  # Create dataframe with the overview of removed sequences
  idx = grep("Overview of removed sequences", info) + 1
  if (any(grepl("Overview of removed sequences", info)[idx-1])) {
    data_to_plot = info[idx:length(info)]
    data_to_plot = read.table(
                    textConnection(data_to_plot), sep = "\t", header = T)

    # Create plots with information of the trimming position
    dat2 = unlist(apply(data_to_plot, 1, function(x) rep(x[1], x[2])))
    dat2 = dat2 + cuted[counter]
    hist(dat2,
         breaks = (max(dat2)-min(dat2)) + 1,
         freq = F,
         bty = "l",
         col = rgb(1,0,0,0.5),
         border = rgb(1,0,0,1),
         xlim = c(min(dat2), max(dat2)),
         xlab = "Position of the sgRNAs within the reads (bp)",
         ylab = "Percentage of reads",
         main = paste("Trimming statistics,", name, sep = " "))
  }

  # Update counter
  counter = counter + 1
}

# Barplot with percentage of trimmed reads per sample
par(xpd=NA, oma=c(0,0,0,2))
xvals = barplot(data.matrix(results[2,]),
              col = "cornflowerblue",
              border = "white",
              space = 0.25,
              axisnames = F,
              xlab = "",
              ylab = "Percentage (%)",
              ylim = c(0,100))
title("Trimmed reads", line = 3)
axis(side = 1, at = xvals, labels = colnames(results),
    cex.axis = 0.6, line = -1, tick = F, las = 2)

# Add number of reads of each sample in the plot
y2 = results[1,] / max(results[1,]) * 100 # Rescale into the bargraph
lines(xvals, y2, lwd = 4, type = "b")
axis(side = 4,
      at = seq(0, 100, length.out = 6),
      labels= round(seq(0, max(results[1,])/1E6, length.out = 6), 3))
mtext("Total reads (Millions)", side=4, line=2.5)
legend(x = "top",
      legend = c("Total reads", rownames(results)[2]),
      bty = "n",
      col = c("black", "cornflowerblue"),
      lty = c(1, 0),
      lwd = c(2, 0),
      pch = c(1, 15),
      horiz = T,
      x.intersp = c(0.8, 0),
      inset = -0.1)
par(xpd=FALSE)
grid(col = "lightgray",
      nx = ncol(results),
      ny = 5)

# Close pdf to save plots
message = dev.off()


##########
## DONE ##
##########
