#!/usr/bin/env Rscript

# Relationship between input parameters and the ones used here
inputargs = commandArgs(trailingOnly=TRUE)
fdr = inputargs[1]
outputdir = inputargs[2]
inputed_ctrls = inputargs[3]
filenames = inputargs[4:length(inputargs)]

# Open pdf
name_pdf = paste(outputdir, "outputs/Selected_genes_MAGeCK.pdf", sep = "")
pdf(name_pdf)


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


####################################
## Analysis of the MAGeCK outputs ##
####################################

# Import controls if provided
if (inputed_ctrls != "") {
  ctrls_file = read.table(inputed_ctrls)
} else {
  ctrls_file = data.frame(V1 = "NoGenesCtrls", V2 = "NoGenesCtrls")
}

for (i in seq(1,length(filenames), by = 2)) {

  # Import data
  genes_pos = read.table(filenames[i],
                          header = T, stringsAsFactors = F, dec =".")
  genes_neg = read.table(filenames[i+1],
                          header = T, stringsAsFactors = F, dec =".")

  # Obtain name of this specific experiment to save outputs
  savename = gsub(".*results_MAGeCK_", "", filenames[i])
  savename = gsub(".gene.*", "", savename)

  # Filter genes with lower than 2 statistical significant gRNAs
  genes_pos_filtered = genes_pos#[-which(genes_pos$goodsgrna<2),]
  genes_neg_filtered = genes_neg#[-which(genes_neg$goodsgrna<2),]

  if (nrow(genes_pos_filtered) == 0) {
    genes_pos_filtered = genes_pos
  }
  if (nrow(genes_neg_filtered) == 0) {
    genes_neg_filtered = genes_neg
  }

  # Re-compute FDR
  genes_pos_filtered$FDR = p.adjust(genes_pos_filtered$p, method = "BH")
  genes_neg_filtered$FDR = p.adjust(genes_neg_filtered$p, method = "BH")

  # If any of the p-values is 0, modify it to avoid future errors
  genes_pos_filtered$p[genes_pos_filtered$p == 0] = 1E-10
  genes_neg_filtered$p[genes_neg_filtered$p == 0] = 1E-10

  # Take top genes to show in the qq-plots
  pos_select = genes_pos_filtered[genes_pos_filtered$FDR < fdr[1],]
  neg_select = genes_neg_filtered[genes_neg_filtered$FDR < fdr[1],]

  # Join positive and negative selected genes for the QQ-plots
  total_genes = rbind(genes_pos_filtered, genes_neg_filtered)
  total_select = rbind(pos_select, neg_select)

  # Add column to hits tables showing % of significant gRNAs per gene
  total_select$Significant_gRNAs = round(
    total_select$goodsgrna/total_select$items_in_group * 100, 2)
  pos_select$Significant_gRNAs = round(
    pos_select$goodsgrna/pos_select$items_in_group * 100, 2)
  neg_select$Significant_gRNAs = round(
    neg_select$goodsgrna/neg_select$items_in_group * 100, 2)

  # Add column to general tables showing % of significant gRNAs per gene
  genes_pos$Significant_gRNAs = round(
    genes_pos$goodsgrna/genes_pos$items_in_group * 100, 2)
  genes_neg$Significant_gRNAs = round(
    genes_neg$goodsgrna/genes_neg$items_in_group * 100, 2)

  # QQ-Plot of positively selected genes
  pos_pvals = genes_pos_filtered$p
  names(pos_pvals) = genes_pos_filtered$group_id
  qq_plot(pos_pvals,
          Experiment_ID = paste("positive selection, test", savename),
          controls = ctrls_file, nlabs = nrow(pos_select))

  # QQ-Plot of negatively selected genes
  neg_pvals = genes_neg_filtered$p
  names(neg_pvals) = genes_neg_filtered$group_id
  qq_plot(neg_pvals,
          Experiment_ID = paste("negative selection, test", savename),
          controls = ctrls_file, nlabs = nrow(neg_select))

  # Boxplot with the percentage of significant gRNAs per gene
  boxplot(list(pos_select$Significant_gRNAs,
               neg_select$Significant_gRNAs,
               genes_pos$Significant_gRNAs,
               genes_neg$Significant_gRNAs),
          horizontal = T,
          outline = F,
          names = c("Pos hits",
                    "Neg hits",
                    "Pos screen",
                    "Neg screen"),
          xlab = "Percentage (%)",
          main = paste("Significant gRNAs per gene, test", savename))
  text(x=median(pos_select$Significant_gRNAs),
       labels = round(median(pos_select$Significant_gRNAs), 2),
       y = 1.5)
  text(x=median(neg_select$Significant_gRNAs),
       labels = round(median(neg_select$Significant_gRNAs), 2),
       y = 2.5)
  text(x=median(genes_pos$Significant_gRNAs),
       labels = round(median(genes_pos$Significant_gRNAs), 2),
       y = 3.5)
  text(x=median(genes_neg$Significant_gRNAs),
       labels = round(median(genes_neg$Significant_gRNAs), 2),
       y = 4.5)

  #########################################
  ## Arrange and save useful information ##
  #########################################

  # Save information of all hits (positive and negative)
  if (nrow(total_select) > 0) {
   # Modify the column with the percentage of significant pgRNAs per gene
   total_select$Significant_gRNAs = paste(
     total_select$Significant_gRNAs, "%", sep = "")

   # Remove unuseful data
   total_select = total_select[-c(2,3,6)]

   # Save useful information
   write.table(total_select,
               file = paste(outputdir,
                 "intermediate/hits.all_MAGeCK_", savename, ".txt", sep = ""),
               row.names = F, col.names = T, quote = F, sep = "\t")
  } else {
    # Create empty table to save
    empty_table = data.frame(
                  group_id = "", p = "", FDR = "", Significant_gRNAs = "")
    write.table(empty_table,
                file = paste(outputdir,
                  "intermediate/hits.all_MAGeCK_", savename, ".txt", sep = ""),
                row.names = F, col.names = T, quote = F, sep = "\t")
  }
}

message = dev.off()


##########
## DONE ##
##########
