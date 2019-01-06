#!/usr/bin/env Rscript

# Relationship between input parameters and the ones used here
inputargs = commandArgs(trailingOnly=TRUE)
fdr = as.numeric(inputargs[1])
outputdir = inputargs[2]
filenames = inputargs[3:length(inputargs)]

# Open pdf
name_pdf = paste(outputdir, "/outputs/Comp_MAGeCK_PBNPA.pdf", sep = "")
pdf(name_pdf)


############################
## Write useful functions ##
############################

# FishersMethod to combine p-values
fishersMethod = function(x) {
  pchisq(-2 * sum(log(x)), df=2*length(x), lower=FALSE)
}

#Create circles for the diagrams
circle <- function(x, y, r, ...) {
  ang <- seq(0, 2*pi, length = 100)
  xx <- x + r * cos(ang)
  yy <- y + r * sin(ang)
  polygon(xx, yy, ...)
}

# Function to create Venn Diagrams
venndia <- function(A, B, test){
  # Compute intersections and unions needed
  uniqueA <- setdiff(A, B)
  uniqueB <- setdiff(B, A)
  intersAB <- intersect(A, B)

  # Compute length of the previous values
  nA <- length(uniqueA)
  nB <- length(uniqueB)
  nAB <- length(intersAB)

  # Plot the circles
  par(mar=c(0, 0, 0, 0))
  plot(-10, -10, ylim=c(0, 9), xlim=c(0, 9), axes=FALSE)
  circle(x=3, y=4, r=3, col=rgb(0,0.9,0.1,.4), border=NA)
  circle(x=6, y=4, r=3, col=rgb(0.2,0.1,0.7,.3), border=NA)

  # Add text
  text(x=c(1.5, 7.4), y=c(5.8, 5.8), c("MAGeCK", "PBNPA"), cex=1.5, col="black")
  text(x=4.5, y=8, paste("Venn Diagram of hits, test", test), cex=1.5)
  if (nAB >= 3){
    text(x=4.5, y=0.5,
         paste("Some common hits: ", intersAB[1], ", ",
               intersAB[2], ", ", intersAB[3], sep = ""), cex=1)
  } else if (nAB == 2){
    text(x=4.5, y=0.5, paste(
      "Common hits: ", intersAB[1], ", ", intersAB[2], sep = ""), cex=1)
  } else if (nAB == 1){
    text(x=4.5, y=0.5, paste("Common hits:", intersAB[1]), cex=1)
  }
  text(
    x=c(2, 7, 4.5), y=c(4, 4, 4),
    c(nA, nB, nAB),
    cex=2)
}

##########################################
## Comparison of MAGeCK & PBNPA outputs ##
##########################################

# Check how many tests have been performed
numtests = length(filenames)/2

# Compare results of MAGeCK and PBNPA for each of the tests
for (i in 1:numtests){
  # Load data
  mageck = read.table(
    filenames[i], header = T, stringsAsFactors = F, dec =".")
  pbnpa = read.table(
    filenames[numtests+i], header = T, stringsAsFactors = F, dec =".")

  # Obtain name of this specific test to save outputs
  savename = gsub(".*results_MAGeCK_", "", filenames[i])
  savename = gsub(".gene.*", "", savename)

  # Take hits
  hits_m = mageck[mageck$neg.fdr < fdr | mageck$pos.fdr < fdr, 1]
  hits_p = pbnpa[pbnpa$neg.fdr < fdr | pbnpa$pos.fdr < fdr, 1]

  # Plot Venn diagram
  venndia(A = hits_m, B = hits_p, test = savename)

  # Mix MAGeCK and PBNPA results into one file
  combination = merge(
    mageck, pbnpa, by.x = "id", by.y = "Gene")[,c(1,4,17,10,15,8)]
  # Combine positive and negative p-values using Fisher's method
  neg.comb = apply(combination[,c(2,3)], 1, fishersMethod)
  pos.comb = apply(combination[,c(4,5)], 1, fishersMethod)
  # Compute FDRs of the combined p-values
  neg.fdr = p.adjust(neg.comb, method = "BH")
  pos.fdr = p.adjust(pos.comb, method = "BH")
  # Generate new dataframe with results
  comb.result = data.frame(Gene = combination$id,
    neg.pval = neg.comb, neg.fdr = neg.fdr,
    pos.pval = pos.comb, pos.fdr = pos.fdr,
    lfc = combination$neg.lfc)

  # Save useful information
  write.table(comb.result,
              file = paste(outputdir,
                           "/outputs/comb_result_", savename, ".txt", sep = ""),
              row.names = F, col.names = T, quote = F, sep = "\t")
}

message = dev.off()
