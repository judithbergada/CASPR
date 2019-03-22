# CASPR Usage


## Input Parameters

For a complete analysis of CRISPR screens, CASPR requires at least three or
four parameters, depending on whether the library is composed of sgRNAs or pgRNAs.
These are as follows:

* `-f` or `--fastq-forward`: sequencing files that contain the forward reads,
for which fastq and fastq.gz formats are accepted. Files must be separated by
commas and written in quotes. Example:
```
-f "ctrl1_1.fastq, d1_forward.fq.gz"
```

* `-r` or `--fastq-reverse`: sequencing files that contain the reverse reads,
only needed to handle pgRNA libraries. If not provided, sgRNAs will be
considered. The accepted formats are fastq and fastq.gz. All files
must be separated by commas and written in quotes, and the order of
the samples must be the same as used in `-f`. Example:
```
-r "ctrl1_2.fastq, d1_reverse.fq.gz"
```

* `-l` or `--library`: text file with 3 or 4 columns, without header.
Column1 shows the IDs of the guide-RNAs. Column2 shows the
targeted genes. Column3 contains the sequence of the 1st sgRNA. Column4 contains
the sequence of the 2nd sgRNA, only if pgRNAs are used. Example:
```
-l pgrna_library.txt
```
```
The file pgrna_library.txt should have this format:

E2F2_1	     E2F2	   GCCGCGGGCCGTGTGAAAGGG	ATAAAGACTGACAGTCTGAA
E2F2_2	     E2F2	   GCAGACGGCGCCTCCCGCAGG	ATAAAGACTGACAGTCTGAA
E2F2_3	     E2F2	   GTTGTGCGATGCCTGCCGGCG	ATAAAGACTGACAGTCTGAA
E2F2_4	     E2F2	   GCCGCGGGCCGTGTGAAAGGG	ACGCTTGTAAGGATCTAGCA
E2F2_5	     E2F2	   GCCGCGGGCCGTGTGAAAGGG	TTAGATGGTTGAGGCCAAGG
ADAM17_1	 ADAM17	   GCTCCCGCCCCCCCATTCCGG	AGAGACTCGGTAAAGACCCA
ADAM17_2	 ADAM17	   GCTCCCGCCCCCCCATTCCGG	GCTCGTTTGAGGGAAGAGCA
ADAM17_3	 ADAM17	   GCTCCCGCCCCCCCATTCCGG	TATGCAGTAATCTCGTTGGG
ADAM17_4	 ADAM17	   GTTTTCGTGACGACAGACGGA	AGAGACTCGGTAAAGACCCA
ADAM17_5	 ADAM17	   GCTCGTTTGAGGGAAGAGCAG	TTTTCGTGACGACAGACGGA
```

* `-e` or `--exper-design`: text file showing the experimental design of the
CRISPR screen. It has three columns, without header. Column1 contains the
different samples in the analysis, named according to their forward fastq files.
Column2 shows the replicates, and each replicate is represented as an integer.
Column3 contains the condition of each sample, "control" or "treated", followed
by the number of test. Example:
```
-e expdesign.txt
```
```
The file expdesign.txt should have this format:

                    ctr1.fastq      1  control1
                    ctr2.fastq      2  control1
                    day30_1.fastq   1  treated1
                    day30_2.fastq   2  treated1
                    ctr3.fastq      3  control2
                    ctr4.fastq      4  control2
                    day30_3.fastq   3  treated2
                    day30_4.fastq   4  treated2

-> Test 1: ctr1 & ctr2 vs day30_1 & day30_2
-> Test 2: ctr3 & ctr4 vs day30_3 & day30_4
```

Additionally, CASPR offers a versatile control of the analysis through several
optional parameters. Some of them are set by default if not modified by the
users; others are ignored if they are not provided. The whole set of optional
arguments is the following:

* `-c` or `--controls`: text file with 2 columns, without header. Column1 shows
the name of the control genes. Column2 shows the type of control, which must be
"Positive", "Negative" or "Neutral". Positive and negative controls are only
used in the plots in order to show these genes with different colors. Neutral
controls are also used during the step of the analysis, and MAGeCK takes them to
predict the null distribution. If controls are not specified, they are not used.
Example:
```
-c controlfile.txt
```
```
The file controlfile.txt should have this format:

E2F2        Negative
ADAM17      Negative
BRCA2       Positive
FAM188A     Positive
AAVS1       Neutral
NonTarget   Neutral
```

* `-a` or `--adapter-f`: sequence of the adapter that must be considered to trim
the forward reads. Characters must be A, T, G or C, as needed to represent
nucleotides. Default: ACCG. Example:
```
-a AGCG
```

* `-A` or `--adapter-r`: sequence of the adapter that must be considered to trim
the reverse reads. Characters must be A, T, G or C, as needed to represent
nucleotides. Default: AAAC. Example:
```
-A ATAC
```

* `-m` or `--mismatches`: maximum number of mismatches to tolerate during the
mapping process. Default: 0. Example:
```
-m 2
```

* `-b` or `--bases-aligned`: minimum number of bases that must be aligned during
the mapping process. Default: 20 bp for sgRNA libraries and 30 bp for pgRNA
libraries. Example:
```
-b 35
```

* `-y` or `--fdr-threshold`: FDR threshold for the identification of hits. It
is used to display the hits in the plots, but it doesn't affect the analysis.
Default: 0.1. Example:
```
-y 0.25
```

* `-o` or `--orientation`: orientation of the paired guides in the plasmid
vector, only needed for pgRNA libraries. If guides are provided such that their
sequence is found in the plasmid vector as 5’-pgRNA-3’, this argument should
be 53. If users have 3’-pgRNA-5’, it should be 35. Other values are not
accepted. Default: 53.  Example:
```
-o 53
```

* `-t` or `--threads`: number of threads that will be used. Default: 1. Example:
```
-t 10
```

* `-q` or `--output-dir`: path to the directory where outputs will be saved. It
must exist and will not be created by CASPR. Default: current directory. Example:
```
-q ./Analysis/survival_screen
```

* `-s` or `--start`: step from which to start the analysis. This argument is
helpful if users want to modify any of the parameters of CASPR and recompute
only part of the analysis, because the previous steps will not be executed again.
Available options are: "qc", "trim", "map" or "test".
Default: start by indexing the genome. Example:
```
-s map
```

* `-p` or `--pause`: step at which to stop the analysis. This step will be
included in the execution. The argument is helpful if users are only interested
in performing a specific part of the analysis, avoiding CASPR to arrive at the
last step. Available options are: "indexing", "qc", "trim", "map".
Default: stop when test is finished. Example:
```
-p map
```

* `-k` or `--keep-tmp`: flag that, when used, forces CASPR to keep the intermediate
files after completing the analysis. This option is very helpful if users want to
repeat any step of the analysis, because all of the intermediate files are saved
and can be used again. Without this option, parameter `-s` doesn't make sense in
a posterior analysis, and CASPR needs to compute all the steps from the beginning.
Default: remove intermediate files.

* `-i` or `--info-alignment`: flag that, when used, forces CASPR to perform
multiple alignments. It is only implemented for the analysis of pgRNA libraries,
with the aim to provide additional information on the mapping process. It
consists on the following steps: first, it maps the reads to the pgRNAs using
parameters specified by the user; second, it tries to align the unmapped reads
from the previous step to only one sgRNA, without tolerating
any mismatches; third, it maps the reads to the pgRNAs considering parameters
`-m` and `-b` to be, respectively, 0 and the length of the pgRNA sequences;
fourth, it maps the reads to the pgRNAs using the same
conditions as in the third step, but changing `-m` to 3. <br />
CASPR provides, by default, some basic statistics and graphs about the trimming and mapping processes. Nonetheless, using this specific flag, `-i`, CASPR generates new plots in which the users can see whether unmapped reads have more mismatches than expected, whether they are too short, whether pgRNAs are composed of repeated sgRNAs, or whether pgRNAs are composed of sgRNAs targeting different genes. <br />
The main drawback of this option, however, is the long computational time required to perform four different alignments. Thus, it is recommended to use `-i` only in cases where mapping outputs are atypical or unclear.

Further information about the usage of CASPR can be obtained using the following flag, or in [CASPR Tutorial](tutorial.md):

* `-h` or `--help`: flag that, when used, shows a help message and exits.



## Outputs

After completing the analysis of a CRISPR screen, CASPR generates a folder
named `outputs`. All the results and final graphs can be found there.
They are as follows:

* **Alignment_statistics.pdf**: It contains a bar-plot that shows, for each sample,
the percentage of uniquely mapped reads, unmapped reads,
and reads that were mapped to more than one guide-RNA.
If flag `-i` was used, it also displays some pie charts with additional
information. Please, check documentation of input parameter `-i` for more details.

* **Trimming_statistics.pdf**: It contains multiple histograms and a bar-plot.
The histograms show, for each of the samples, the positions in which the
adapters were found within the reads. The bar-plot indicates the percentage of
reads that could be trimmed in each sample, and the total of number of
reads per sample.

* **results_MAGeCK_{{idx}}.gene_summary.txt**: It is a table with the p-values
and FDRs of all genes computed by MAGeCK. There is one file per test
that has been carried out.

* **results_PBNPA_{{idx}}_summary.txt**: It is a table with the p-values and
FDRs of all genes computed by PBNPA. There is one file per test
that has been carried out.

* **comb_result_{{idx}}.txt**: It is a table with the Fisher combined p-values of
MAGeCK and PBNPA, along with the corresponding FDRs. There is one file
per test that has been carried out.

* **Selected_genes_MAGeCK.pdf**: It contains two quantile-quantile plots
per test. QQ-plots are based on the p-values computed by MAGeCK for the
positive and negative selection.

* **Selected_genes_PBNPA.pdf** : It contains two quantile-quantile plots
per test. QQ-plots are based on the p-values computed by PBNPA for the
positive and negative selection.

* **Comp_MAGeCK_PBNPA.pdf**: It contains one Venn-diagram per test,
showing the overlap between hits detected by MAGeCK and PBNPA.

* **Log_counts.pdf**: It contains five plots per test: the first one shows
the log2-fold-change of all sgRNAs or pgRNAs; the second one compares the
number of read counts per guide of untreated versus treated samples;
the third one displays the cumulative counts of all samples; the fourth and
fifth graphs are volcano plots showing the log2-fold change of all genes versus
the Fisher combined p-values of MAGeCK and PBNPA. In the first volcano plot,
genes detected as hits either by MAGeCK or PBNPA are colored. In the second
volcano plot, all genes classified as hits according to the Fisher combined
statistics are colored.

* **table.counts.txt**: It is a table with the quantification of the sgRNAs or
pgRNAs. First and second columns correspond to the IDs of the guides and the
targeted genes, respectively. The rest of the columns show the number of
read counts per guide of the different samples
-there is one column per sample-.

* **inputs.txt**: It is a text file with information of the input parameters
that have been used by CASPR to get the outputs provided.

In [CASPR Tutorial](tutorial.md) users can
find examples of the output figures that should be expected.

Furthermore, CASPR generates three additional folders:

* `qualitycontrol`: It shows the quality control of the reads.

* `genome`: It contains the indexed library of guide-RNAs, which is needed
for the mapping process.

* `intermediate`: It contains intermediate files such as the trimmed reads,
the unmapped reads, the trimming and mapping statistics used to generate output
figures, and other files created by MAGeCK and PBNPA during the analysis. By
default, when flag `k` is not used, this folder is removed.
