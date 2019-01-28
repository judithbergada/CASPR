# Installation Guide


## Get started: Download CASPR Code and Data

CASPR is being developed and tested under Linux and macOS.

In order to use it, first you need to clone the repository:

```bash
git clone https://github.com/judithbergada/CASPR.git $HOME/CASPR
```

For a proper installation, it is also recommended to run:

```bash
# Add CASPR to your PATH
echo "export PATH=$HOME/CASPR/source:$PATH" >> $HOME/.bashrc
source $HOME/.bashrc
# Please, note that the previous commands do not work if you use Zsh
```

At this point, the tool should be available in your computer.
You can check it:

```bash
CASPR --help
```


## Conda Environment Installation

CASPR employs several software packages to carry out its different functions:

*   `fastqc`
*   `cutadapt`
*   `fastx-toolkit`
*   `STAR`
*   `samtools`
*   `mageck` and `mageck-vispr`
*   `vispr`
*   `R`

The easiest way to download them is to use
[Miniconda](https://conda.io/miniconda.html), which allows to create an
isolated environment. To install miniconda, please follow the
[instuctions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Once Miniconda is installed, you can obtain the required packages with a single
command like this:

```bash
conda env create -n caspr -f environment.yml
```

This option is recommended because it will install all the tools defined in the
`environement.yml` file, ensuring that CASPR works smoothly.

Alternatively, you can also run:

```bash
# Create environment
conda create -n caspr python=3.5
# Install required software
conda install -n caspr --yes r
conda install -n caspr_new -c bioconda --yes \
  cutadapt star samtools vispr mageck mageck-vispr fastqc
conda install -n caspr_new -c biobuilds --yes fastx-toolkit
```


!!! warning "Installation on macOS"

    Since macOS uses custom versions of `sed` and `grep`, it is also necessary to
    install the latest GNU version to run CASPR. For that, first you need to
    have [HomeBrew](https://brew.sh/).

    Once HomeBrew works at your device, you only need to copy the following commands
    on your terminal:

    ```bash
    # Install sed and grep
    brew install gnu-sed grep
    # Make sure the installed versions are the first ones available in your PATH
    echo "export PATH=\"/usr/local/opt/gnu-sed/libexec/gnubin:\$PATH\"" >> ~/.bashrc
    echo "export PATH=\"/usr/local/opt/grep/libexec/gnubin:\$PATH\"" >> ~/.bashrc
    source ~/.bashrc
    ```

Now, you can start using CASPR on Linux or macOS!

!!! info "Important"

    To use the tool, please do not forget to activate the Conda Environment!

    ```bash
    conda activate caspr
    ```


## Users of IBU cluster

Users of IBU cluster (Interfaculty Bioinformatics Unit, University of Bern)
do not need to follow the guide to install a Conda Environment.

Instead, you can load the required software packages,
which are all installed in the cluster. This can be done as follows:

```bash
# Import software: only for users of IBY cluster
module add UHTS/Analysis/fastx_toolkit/0.0.13.2
module add UHTS/Aligner/STAR/2.6.0c
module add UHTS/Analysis/samtools/1.4
module add R/3.5.1
module add UHTS/Quality_control/mageck-vispr/0.5.4
```
Neverthless, installing the Conda Environment would have an advantage:
VISPR would display the results without the need of additional commands
(see the user cases in the tutorial for more information).
