# Installation Guide

CASPR is being developed and tested under Linux and macOS.

In order to use it, first you need to clone the repository:

```bash
git clone https://github.com/judithbergada/CASPRg.git $HOME/CASPR
```

For a proper installation, it is also recommended to run:

```bash
# Add CASPR to your PATH
echo "export PATH=$HOME/CASPR/source:$PATH" >> $HOME/.bashrc
source $HOME/.bashrc
# Please, note that the previous commands do not work if you use Zsh
```

At this point, the tool should be available in your computer.

## Conda Environment Installation

CASPR employs several software packages to carry out its different functions:

* `fastqc`
* `cutadapt`
* `fastx-toolkit`
* `STAR`
* `samtools`
* `mageck` and `mageck-vispr`
* `vispr`
* `R`

The easiest way to download them is to use [Miniconda](https://conda.io/miniconda.html), which allows to create an isolated environment. To install miniconda, please follow the [instuctions](https://conda.io/docs/user-guide/install/index.html).

Once Miniconda is installed, you can obtain all the packages from a single command like this:

```bash
conda env create -n caspr -f environment.yml
```

This option is recommended because it will install the software defined in the `environement.yml` file, ensuring that the pipeline works smoothly.

Alternatively, you can also run:

```bash
conda create -n caspr python=3.5
conda install -n caspr -c bioconda cutadapt
conda install -n caspr -c bioconda star
conda install -n caspr -c bioconda samtools
conda install -n caspr -c bioconda vispr
conda install -n caspr -c bioconda mageck
conda install -n caspr -c bioconda mageck-vispr
conda install -n caspr r
conda install -n caspr -c biobuilds fastx-toolkit
conda install -n caspr -c bioconda fastqc
```

Once the installation is successful, you can start using CASPR.
