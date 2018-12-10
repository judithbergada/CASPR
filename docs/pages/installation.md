# Installation Guide

To run the CASPRg pipeline, some tools are software are required to run the different steps. The list of software required is the following:

* `fastqc`
* `cutadapt`
* `fastx-toolkit`
* `STAR`
* `samtools`
* `mageck` and `mageck-vispr`
* `vispr`
* `R`

The best way to install this software is using the [Miniconda](https://conda.io/miniconda.html) which will allow us to create an isolated environment where to install all this software. To install miniconda, please follow the [instuctions](https://conda.io/docs/user-guide/install/index.html) given the operative system and the architecture of your computer.

First step is to clone the repository:

```bash
git clone https://github.com/judithbergada/CASPRg.git
```

## Conda Environment Installation

```bash
conda create -n casprg python=3.5
conda install -n casprg -c bioconda cutadapt
conda install -n casprg -c bioconda star
conda install -n casprg -c bioconda samtools
conda install -n casprg -c bioconda vispr
conda install -n casprg -c bioconda mageck
conda install -n casprg -c bioconda mageck-vispr
conda install -n casprg r
conda install -n casprg -c biobuilds fastx-toolkit
conda install -n casprg -c bioconda fastqc
```

To install all the packages from a single command, run:

```bash
conda env create -n casprg -f environment.yml
```

This will install all the required packages defined in the `environement.yml` file which ensures that the pipeline will run smoothly.

## R Packages

Problems with R. Tracking bug: https://github.com/conda/conda/issues/6183.

```
pbnpa
PBNPA
```
