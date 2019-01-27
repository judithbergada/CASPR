# Installation Guide

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
[instuctions](https://conda.io/docs/user-guide/install/index.html).

Once Miniconda is installed, you can obtain all the packages from a single
command like this:

```bash
conda env create -n caspr -f environment.yml
```

This option is recommended because it will install the software defined in the
`environement.yml` file, ensuring that the pipeline works smoothly.

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

!!! info "Activate environment"

    Once you have installed all the required software don't forget to
    activate the environment to start using CASPR:

    ```bash
    conda activate caspr
    ```

Now, you can start using CASPR.

!!! warning "Installation on macOS"

    Because macOS uses custom versions of `sed` and `grep` it is necessary to
    install the latest GNU version to run CASPR. To install it, you need to have
    [HomeBrew](https://brew.sh/):

    ```bash
    # Install sed and grep
    brew install gnu-sed grep
    # Make sure the installed versions are the first ones available in your PATH
    echo "export PATH=\"/usr/local/opt/gnu-sed/libexec/gnubin:\$PATH\"" >> ~/.bashrc
    echo "export PATH=\"/usr/local/opt/grep/libexec/gnubin:\$PATH\"" >> ~/.bashrc
    source ~/.bashrc
    ```
