#!/usr/bin/env Rscript

# Relationship between input parameters and the ones used here
files_total = commandArgs(trailingOnly=TRUE)

#######################################
## Plot statistics of unmapped reads ##
#######################################

# Open pdf
outputdir = files_total[1]
name_pdf = paste(outputdir, "outputs/Alignment_unmapped.pdf", sep = "")
pdf(name_pdf)

# Create dataframe to save the results
results = data.frame(numeric(4))
rownames(results) = c(
  "Total reads", "Unique", "Non-unique    ", "Unmapped")

# Analysis: from the previously unmapped reads, check how many
# of them align to 1 of the gRNAs + 5bp of the second gRNA.
counter = 1
for (filename in files_total[2:(((length(files_total)-1)/2)+1)]) {

  # Obtain SampleName
  name=gsub(".*Statistics_unmapped_pgrna_", "", filename)
  name=gsub(".txt", "", name)

  # Read file
  info = readLines(filename)

  # Obtain number of input reads
  readsnumb = grep("Number of input reads", info, value = TRUE)
  assign(paste("readsnumb_", name, sep = ""),
         as.numeric(
           unlist(regmatches(
             readsnumb, gregexpr(
               "[[:digit:]]+", readsnumb)))))

  # Obtain percentage of unique reads
  unique = grep("Uniquely mapped reads %", info, value = TRUE)
  assign(paste("unique_", name, sep = ""),
         as.numeric(
           unlist(regmatches(
             unique, gregexpr(
               "[[:digit:]]+\\.[[:digit:]]+", unique)))))

  # Obtain percentage of non-unique reads
  multiple = grep("% of reads mapped to", info, value = TRUE)
  assign(paste("multiple_", name, sep = ""),
         as.numeric(
           unlist(
             regmatches(
               multiple, gregexpr(
                 "[[:digit:]]+\\.[[:digit:]]+", multiple)))))

  # Obtain percentage of unmapped reads
  unmapped = grep("% of reads unmapped", info, value = TRUE)
  assign(paste("unmapped_", name, sep = ""),
         as.numeric(
           unlist(
             regmatches(
               unmapped, gregexpr(
                 "[[:digit:]]+\\.[[:digit:]]+", unmapped)))))

  # Update dataframe with the results
  results[,counter] = c(get(paste("readsnumb_",name, sep = "")),
                        get(paste("unique_",name, sep = "")),
                        sum(get(paste("multiple_",name, sep = ""))),
                        sum(get(paste("unmapped_",name, sep = ""))))
  colnames(results)[counter] = name

  # Update counter
  counter = counter + 1
}

# Stacked Bar Plot with Colors and Legend
colors_data = c("aquamarine", "coral", "cornflowerblue")
par(xpd=NA, oma=c(0,0,0,2))
xvals = barplot(data.matrix(results[2:4,]),
              col = colors_data,
              border = "white",
              space = 0.25,
              axisnames = F,
              xlab = "",
              ylab = "Percentage (%)")
title("Unmapped reads mapping 1gRNA+5bp of the second", line = 3)
axis(side = 1, at = xvals, labels = colnames(results),
      line = -1, tick = F, las = 2)

# Add number of reads of each sample
y2 = results[1,] / max(results[1,]) * 100 # Rescale into the bargraph
lines(xvals, y2, lwd = 4, type = "b")
axis(side = 4,
      at = seq(0, 100, length.out = 6),
      labels=round(seq(0, max(results[1,])/1E6, length.out = 6), 3))
mtext("Total reads (Millions)", side=4, line=2.5)
legend(x = "top",
      legend = c("Total reads", rownames(results)[2:4]),
      bty = "n",
      col = c("black", colors_data),
      lty = c(1, 0, 0, 0),
      lwd = c(2, 0, 0, 0),
      pch = c(1, 15, 15, 15),
      horiz = T,
      x.intersp = c(0.8, 0, 0, 0),
      inset = -0.1)

# Create dataframe to save the results
results = data.frame(numeric(4))
rownames(results) = c(
  "Total reads", "Unique", "Non-unique    ", "Unmapped")

# Analysis: from the previously unmapped reads, check how many
# of them align ONLY to 1 gRNA
counter = 1
for (filename in files_total[
      (((length(files_total)-1)/2)+2):length(files_total)]) {

  # Obtain SampleName
  name=gsub(".*Statistics_unmapped_sgrna_", "", filename)
  name=gsub(".txt", "", name)

  # Read file
  info = readLines(filename)

  # Obtain number of input reads
  readsnumb = grep("Number of input reads", info, value = TRUE)
  assign(paste("readsnumb_", name, sep = ""),
         as.numeric(
           unlist(regmatches(
             readsnumb, gregexpr(
               "[[:digit:]]+", readsnumb)))))

  # Obtain percentage of unique reads
  unique = grep("Uniquely mapped reads %", info, value = TRUE)
  assign(paste("unique_", name, sep = ""),
         as.numeric(
           unlist(regmatches(
             unique, gregexpr(
               "[[:digit:]]+\\.[[:digit:]]+", unique)))))

  # Obtain percentage of non-unique reads
  multiple = grep("% of reads mapped to", info, value = TRUE)
  assign(paste("multiple_", name, sep = ""),
         as.numeric(
           unlist(
             regmatches(
               multiple, gregexpr(
                 "[[:digit:]]+\\.[[:digit:]]+", multiple)))))

  # Obtain percentage of unmapped reads
  unmapped = grep("% of reads unmapped", info, value = TRUE)
  assign(paste("unmapped_", name, sep = ""),
         as.numeric(
           unlist(
             regmatches(
               unmapped, gregexpr(
                 "[[:digit:]]+\\.[[:digit:]]+", unmapped)))))

  # Update dataframe with the results
  results[,counter] = c(get(paste("readsnumb_",name, sep = "")),
                        get(paste("unique_",name, sep = "")),
                        sum(get(paste("multiple_",name, sep = ""))),
                        sum(get(paste("unmapped_",name, sep = ""))))
  colnames(results)[counter] = name

  # Update counter
  counter = counter + 1
}

# Stacked Bar Plot with Colors and Legend
colors_data = c("aquamarine", "coral", "cornflowerblue")
par(xpd=NA, oma=c(0,0,0,2))
xvals = barplot(data.matrix(results[2:4,]),
              col = colors_data,
              border = "white",
              space = 0.25,
              axisnames = F,
              xlab = "",
              ylab = "Percentage (%)")
title("Unmapped reads mapping only 1 gRNA", line = 3)
axis(side = 1, at = xvals, labels = colnames(results),
      line = -1, tick = F, las = 2)

# Add number of reads of each sample
y2 = results[1,] / max(results[1,]) * 100 # Rescale into the bargraph
lines(xvals, y2, lwd = 4, type = "b")
axis(side = 4,
      at = seq(0, 100, length.out = 6),
      labels=round(seq(0, max(results[1,])/1E6, length.out = 6), 3))
mtext("Total reads (Millions)", side=4, line=2.5)
legend(x = "top",
      legend = c("Total reads", rownames(results)[2:4]),
      bty = "n",
      col = c("black", colors_data),
      lty = c(1, 0, 0, 0),
      lwd = c(2, 0, 0, 0),
      pch = c(1, 15, 15, 15),
      horiz = T,
      x.intersp = c(0.8, 0, 0, 0),
      inset = -0.1)

message = dev.off()

##########
## DONE ##
##########
