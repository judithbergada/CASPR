#!/usr/bin/env Rscript

# Relationship between input parameters and the ones used here
inputdata = commandArgs(trailingOnly=TRUE)
controls = inputdata[1]
fdr_required = inputdata[2]
countsfile = inputdata[3]
expdesign = inputdata[4]
idx_sample = inputdata[5]
outputdir = inputdata[6]
controlsfile = inputdata[7]

# Check if PBNPA package is properly installed and load it
if (!require("PBNPA",character.only = TRUE)) {
  install.packages("PBNPA", repos = "https://stat.ethz.ch/CRAN/")
}
library("PBNPA")

# Open pdf
pdf(paste(
  outputdir, "intermediate/Selected_PBNPA_", idx_sample, ".pdf", sep = ""))

############################
## Write useful functions ##
############################

# Function to create QQ-plots
qq_plot <- function(p_val, Experiment_ID, controls, ci = 0.95, nlabs = 4) {
  # Define needed parameters
  n = length(p_val)
  observed = -log10(sort(p_val))
  expected = -log10(ppoints(n))
  clower = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1))
  cupper = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  log10Pe = expression(paste("Expected -log"[10], plain((p-values))))
  log10Po = expression(paste("Observed -log"[10], plain((p-values))))
  # Generate qq-plot
  plot(expected, observed, pch = 1, bty = "l",
       xlab = log10Pe, ylab = log10Po,
       xlim = c(0, max(expected)),
       ylim = c(0, max(c(observed, expected))))
  title(paste("QQ-plot of p-values,", Experiment_ID), line = 2)
  if (nlabs > 0) {
    text(expected, observed,
         label = ifelse(
           observed %in% observed[1:nlabs],names(observed)[1:nlabs],""),
         pos = 1, srt = 90, col = "blue", cex = 0.4)
  }
  # Add confidence intervals
  abline(a = 0, b = 1, col = "red", lwd = 2, xlim = c(0, max(expected)))
  lines(expected, cupper, lty = 2, col = "red", xlim = c(0, max(expected)))
  lines(expected, clower, lty = 2, col = "red", xlim = c(0, max(expected)))

  # Add controls if input was given
  positive_ctrls = controls[controls$V2 == "Positive", 1]
  points(expected[names(observed) %in% positive_ctrls],
         observed[names(observed) %in% positive_ctrls],
         pch = 19, col = "green")
  negative_ctrls = controls[controls$V2 == "Negative", 1]
  points(expected[names(observed) %in% negative_ctrls],
         observed[names(observed) %in% negative_ctrls],
         pch = 19, col = "yellow")
  neutral_ctrls = controls[controls$V2 == "Neutral", 1]
  points(expected[names(observed) %in% neutral_ctrls],
         observed[names(observed) %in% neutral_ctrls],
         pch = 19, col = "red")
  legend(x = "topleft", bty = "n", pch = 19,
         legend = c("Positive controls",
                    "Negative controls",
                    "Neutral controls"),
         col = c("green", "yellow", "red"))
}

###########################################
## Analysis of the data using PBNPA in R ##
###########################################

# Read input files
dat = read.table(countsfile, header = T, stringsAsFactors = F)
design = read.table(expdesign, header = F, stringsAsFactors = F)
if (controlsfile != "") {
  ctrls_file = read.table(controlsfile)
} else {
  ctrls_file = data.frame(V1 = "NoGenesCtrls", V2 = "NoGenesCtrls")
}

# Put data in the required format for PBNPA
all_info = list()
counter = 1
for (k in unlist(strsplit(controls, split=","))) {
  sample_idx = design[design$V1 == k, 2]
  treatm_name = design[
    design$V2 == sample_idx & grepl(
      paste("treated", idx_sample, sep = ""), design$V3),1]
  ctr = dat[,colnames(dat) == k]
  exp = dat[,colnames(dat) == treatm_name]
  to_analyse = data.frame(
    sgRNA = dat[,1], Gene = dat[,2], contrl = ctr, treatm = exp)
  all_info[[counter]] = to_analyse
  counter = counter + 1
}

# Perform PBNPA analysis
result = PBNPA(all_info, fdr = fdr_required)

# If any of the p-values is 0, modify it to avoid future errors
result$final.result$pos.pvalue[result$final.result$pos.pvalue == 0] = 1E-10
result$final.result$neg.pvalue[result$final.result$neg.pvalue == 0] = 1E-10


######################
## Show the results ##
######################

# QQ-Plot of positively selected genes
pos_pvals = result$final.result$pos.pvalue
names(pos_pvals) = result$final.result$Gene
qq_plot(pos_pvals,
        Experiment_ID = paste("positive selection, test", idx_sample),
        controls = ctrls_file, nlabs = result$pos.no)

# QQ-Plot of negatively selected genes
neg_pvals = result$final.result$neg.pvalue
names(neg_pvals) = result$final.result$Gene
qq_plot(neg_pvals,
        Experiment_ID = paste("negative selection, test", idx_sample),
        controls = ctrls_file, nlabs = result$neg.no)

# Save information of all hits (positive and negative)
if (result$pos.no + result$neg.no > 0) {
  info_needed = rbind(
    result$final.result[
      result$final.result$Gene %in% result$pos.gene,],
    result$final.result[
      result$final.result$Gene %in% result$neg.gene,])
  # Save useful information
  write.table(info_needed,
              file = paste(outputdir,
                "intermediate/hits.all_PBNPA_", idx_sample, ".txt", sep = ""),
              row.names = F, col.names = T, quote = F, sep = "\t")
} else {
  # Create empty table to save
  empty_table = data.frame(
                  Gene = "", pos.pvalue = "",
                  pos.fdr = "", neg.pvalue = "", neg.fdr = "")
  write.table(empty_table,
              file = paste(outputdir,
                "intermediate/hits.all_PBNPA_", idx_sample, ".txt", sep = ""),
              row.names = F, col.names = T, quote = F, sep = "\t")
}

# Save all the results
write.table(result$final.result,
            file = paste(outputdir,
              "outputs/results_PBNPA_", idx_sample, "_summary.txt", sep = ""),
            row.names = F, col.names = T, quote = F, sep = "\t")

# Close pdf
message = dev.off()

###########
## DONE ##
##########
