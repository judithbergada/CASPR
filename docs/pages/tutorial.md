# Tutorial

Running CASPR is extremely easy and convenient to analyze CRIPR-Cas9 screens using pgRNAs. The test data folder contains two brief examples to go through all steps. Simply copy the scripts provided in each example to visualize the demos.

## Get started: download CASPR code and data

If it is your first time with CASPR, please make sure you have downoaded the repository. This can be done a follows:

```bash
git clone https://github.com/judithbergada/CASPR.git $HOME/CASPR
```

For a proper installation of the tool, it is also recommended to run:

```bash
# Add CASPR to PATH
echo "export PATH=$HOME/CASPR/source:$PATH" >> $HOME/.bashrc
source $HOME/.bashrc
```

You can now check that CASPR is available:

```bash
CASPR --help
```

## Example 1: analysis of a CRISPR screen using pgRNAs

This example shows how to analyse a CRISPR screen with the following features:

* `Library of pgRNAs`
* `Two time points: week zero and week four`
* `Only one replicate per time point`
* `Start analysis from the raw sequencing data`

To compute the analysis, it is essential to have:

* `Paired reads`
* `Text file with a library of pgRNAs`
* `Text file with the experimental design (control and treated samples)`

Additionally, you may provide a text file with positive, negative and neutral controls.

Let's check that the data is available in the expected formats:

```bash
# Check that the data has been sucessfully downoaded
ls $HOME/CASPR/testdata/example_pgrna/

# Display the content of the files to see the format
cat $HOME/CASPR/testdata/example_pgrna/expdesign.txt
head $HOME/CASPR/testdata/example_pgrna/library.txt
head $HOME/CASPR/testdata/example_pgrna/controlfile.txt
zcat < $HOME/CASPR/testdata/example_pgrna/week0.1_forward.fastq.gz | head
```

Now you are ready to run the test. If you are not sure about the commands, you can copy them immediately from here:

```bash
# Create output directory
mkdir firstexample

# Run the tool
CASPR \
-f "$HOME/CASPR/testdata/example_pgrna/week0.1_forward.fastq.gz, $HOME/CASPR/testdata/example_pgrna/week4.1_forward.fastq.gz" \
-r "$HOME/CASPR/testdata/example_pgrna/week0.1_reverse.fastq.gz, $HOME/CASPR/testdata/example_pgrna/week4.1_reverse.fastq.gz" \
--library $HOME/CASPR/testdata/example_pgrna/library.txt \
-y 0.25 -k \
--output-dir ./firstexample \
--controls $HOME/CASPR/testdata/example_pgrna/controlfile.txt \
--exper-design $HOME/CASPR/testdata/example_pgrna/expdesign.txt
```

Note that, if you copy the commands from above, you are changing two of the default parameters of CASPR: (1) the FDR threshold -y will be set to 0.25, and (2) the intermediate files will be kept. Please, make sure you keep the intermediate files to follow the second part of the example.

After a few minutes, the outputs should appear in your computer. At this point, you can open them and see if they are as expected. The expected outputs are found in:

```bash
ls $HOME/CASPR/testdata/example_pgrna/expected_outputs/
```

Moreover, you can benefit from VISPR to get interactive results on you data. Please try the following:

```bash
vispr server ./firstexample/config*
```
If you are working from a cluster, you will also need to run this command locally on your computer:

```bash
ssh -f {user}@binfservms01.unibe.ch -L 5000:localhost:5000 -N
```
Now, copy the webpage on Internet and enjoy visualizing the data.

If you are still curious about other options of CASPR, let's try to get more information on the unmapped reads. This can be done esily just by adding the tag -i to the previous command. Furthermore, if you kept the intermediate files with the -k argument, neither the previous nor the next steps will be necessary anymore.

You can find a solution to quickly finish the example here:

```bash
CASPR \
-f "$HOME/CASPR/testdata/example_pgrna/week0.1_forward.fastq.gz, $HOME/CASPR/testdata/example_pgrna/week4.1_forward.fastq.gz" \
-r "$HOME/CASPR/testdata/example_pgrna/week0.1_reverse.fastq.gz, $HOME/CASPR/testdata/example_pgrna/week4.1_reverse.fastq.gz" \
--library $HOME/CASPR/testdata/example_pgrna/library.txt \
-i --start map --pause map -y 0.25 -k \
--output-dir ./firstexample \
--controls $HOME/CASPR/testdata/example_pgrna/controlfile.txt \
--exper-design $HOME/CASPR/testdata/example_pgrna/expdesign.txt
```

As before, your outputs should look similar to the ones provided. Check that everything worked out.

If you arrived here, you finished the first example successfully!

## Example 2: analysis of a CRISPR screen using sgRNAs

This example shows how to analyse a CRISPR screen with the following features:

* `Library of sgRNAs`
* `Two time points: week zero and week four`
* `Two replicates per time point`
* `Start analysis from the raw sequencing data`

To compute the analysis, it is essential to have:

* `Single-read sequencing`
* `Text file with a library of sgRNAs`
* `Text file with the experimental design (control and treated samples)`

Moreover, you may provide a text file with positive, negative and neutral controls, if interested.

Before starting the test, check that the data is available in the expected formats:

```bash
cat $HOME/CASPR/testdata/example_sgrna/expdesign.txt
head $HOME/CASPR/testdata/example_sgrna/library.txt
head $HOME/CASPR/testdata/example_sgrna/controlfile.txt
zcat < $HOME/CASPR/testdata/example_sgrna/week0.1_forward.fastq.gz | head
```

At this point, you are probably ready to try the analysis yourself. To get nicer visualization of the outputs, it is recommended to use an FDR threshold of 0.25.

If you need further help, you can also use these commands:

```bash
# Create output directory
mkdir secondexample

# Run the tool
CASPR \
-f "$HOME/CASPR/testdata/example_sgrna/week0.1_forward.fastq.gz, $HOME/CASPR/testdata/example_sgrna/week0.2_forward.fastq.gz, $HOME/CASPR/testdata/example_sgrna/week4.1_forward.fastq.gz, $HOME/CASPR/testdata/example_sgrna/week4.2_forward.fastq.gz" \
--library $HOME/CASPR/testdata/example_sgrna/library.txt \
-y 0.25 \
--output-dir ./secondexample \
--controls $HOME/CASPR/testdata/example_sgrna/controlfile.txt \
--exper-design $HOME/CASPR/testdata/example_sgrna/expdesign.txt
```

Quickly, the outputs should appear in your computer. You can open them and see if they look as expected. The expected outputs are found in:

```bash
ls $HOME/CASPR/testdata/example_sgrna/expected_outputs/
```

Finally, let's perfom the test again using an FDR threshold of 0.1. For that, you do not need to execute all the previous steps, only the assessment of gene significance.

It can be done like this:

```bash
CASPR \
--start test -y 0.1 \
--output-dir ./secondexample \
--controls $HOME/CASPR/testdata/example_sgrna/controlfile.txt \
--exper-design $HOME/CASPR/testdata/example_sgrna/expdesign.txt
```

The results should be exactly the same as above, but fewer hits are displayed in the plots. This option ensures that you can perform analyses of CRISPR screens not only starting rom the raw sequencing data, but also beginning from a table of read counts.

If you arrived here, you are totally ready to use CASPR with your CRISPR data. Good luck!

## Ubelix cluster: useful information

CASPR employs several software packages to carry out its different functions, but working from the Ubelix cluster has some advantages: all of the required tools are already installed. Therefore, please make sure you load everything in advance:

```bash
# Import software: it only works for users of Ubelix cluster
module add UHTS/Analysis/fastx_toolkit/0.0.13.2
module add UHTS/Aligner/STAR/2.6.0c
module add UHTS/Analysis/samtools/1.4
module add R/3.5.1
module add UHTS/Quality_control/mageck-vispr/0.5.4
```
