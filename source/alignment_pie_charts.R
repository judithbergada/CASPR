#!/usr/bin/env Rscript

# Relationship between input parameters and the ones used here
inputdata = commandArgs(trailingOnly=TRUE)

# Basic variables
name = inputdata[1]
q = inputdata[2]
b = as.numeric(inputdata[3])
tot_len = as.numeric(inputdata[4])
# Graph 1
mapped = as.numeric(inputdata[5])
unmapped = as.numeric(inputdata[6])
# Graph 2
mappedm0 = as.numeric(inputdata[7])
mappedm3 = as.numeric(inputdata[8])
short = as.numeric(inputdata[9])
mappedmore3 = as.numeric(inputdata[10])
# Graph 3
mapped20 = as.numeric(inputdata[11])
nonunique20 = as.numeric(inputdata[12])
unmapped20 = as.numeric(inputdata[13])
# Graph 4
repeatedg = as.numeric(inputdata[14])
recombing = as.numeric(inputdata[15])
others = as.numeric(inputdata[16])

# Open pdf
name_pdf = paste(q, "/intermediate/Alignment_stat2_", name, ".pdf", sep = "")
pdf(name_pdf)
par(mfrow = c(2,2))

######################
## Plot pie chart 1 ##
######################

# Arrange slices and labels
slices <- c(mapped, unmapped)
lbls <- c(paste("Mapped ", b, "bp", sep =""),
          paste("Unmapped ", b, "bp", sep = ""))
pct <- round(slices/sum(slices)*100)
# Plot
pie(slices, labels = paste(pct, "%", sep = ""),
    col=rainbow(length(lbls)),
    main = paste("Mapping at least ", b, "bp", sep = ""),
    cex = 0.8)
par(xpd=TRUE)
legend(x = "bottom", bty = "n", pch = 19,
       legend = lbls, ncol = 2, inset = -0.1,
       col = rainbow(length(lbls)), cex = 0.7)
par(xpd=FALSE)

######################
## Plot pie chart 2 ##
######################

# Arrange slices and labels
slices <- c(mappedm0, mappedm3, short, mappedmore3)
lbls <- c(paste("Map ", tot_len, "bp, m = 0", sep = ""),
          paste("Map ", tot_len, "bp, 1<=m<=3", sep = ""),
          paste("Read size < ", tot_len-3, "bp", sep = ""),
          paste("Map ", tot_len, "bp, m > 3", sep = ""))
pct <- round(slices/sum(slices)*100)
# Plot
pie(slices, labels = paste(pct, "%", sep = ""),
    col=rainbow(length(lbls)),
    main = paste("Map at least ", b, "bp", sep = ""),
    cex = 0.8)
par(xpd=TRUE)
legend(x = "bottom", bty = "n", pch = 19,
       legend = lbls, ncol = 2, inset = -0.2,
       col = rainbow(length(lbls)), cex = 0.7)
par(xpd=FALSE)

######################
## Plot pie chart 3 ##
######################

# Arrange slices and labels
slices <- c(mapped20, nonunique20, unmapped20)
lbls <- c("Map 1 guide",
          "Map 1 guide multiple times",
          "Map < 1 guide")
pct <- round(slices/sum(slices)*100)
# Plot
pie(slices, labels = paste(pct, "%", sep = ""),
    col=rainbow(length(lbls)),
    main = paste("Mapping unmapped reads to ",
                 round(tot_len/2), "bp", sep = ""),
    cex = 0.8)
par(xpd=TRUE)
legend(x = "bottom", bty = "n", pch = 19,
       legend = lbls, ncol = 2, inset = -0.2,
       col = rainbow(length(lbls)), cex = 0.7)
par(xpd=FALSE)

######################
## Plot pie chart 4 ##
######################

# Arrange slices and labels
slices <- c(repeatedg, recombing, others)
lbls <- c("Repeated sgRNA", "Recombination", "Other")
pct <- round(slices/sum(slices)*100)
# Plot
pie(slices,labels = paste(pct, "%", sep = ""),
    col=rainbow(length(lbls)),
    main="Map a guide multiple times", cex = 0.8)
par(xpd=TRUE)
legend(x = "bottom", bty = "n", pch = 19,
       legend = lbls, ncol = 2, inset = -0.2,
       col = rainbow(length(lbls)), cex = 0.7)
par(xpd=FALSE)

# Add title to the page and close pdf
mtext(paste("Fastq files: ", name, sep=""),
      outer=TRUE, cex=1, line=-1.5)
message = dev.off()

##########
## DONE ##
##########
