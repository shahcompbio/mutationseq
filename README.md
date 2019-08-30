# Mutationseq


## Install
1. Create a directory to store temporary files:
    ```
    mkdir -p $HOME/museq/
    cd $HOME/museq
    ```
2. Download and install Miniconda for python 2.7:

    Download miniconda
    ```
    wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
    ```
    Run
    ```
    sh Miniconda2-latest-Linux-x86_64.sh
    ```
    and follow the instructions. When you have finished following the instructions, you should have python installed:
    ```
    which python
    ~/miniconda2/bin/python
    ```
3. clone the mutationseq repo.
    ```
    git clone https://github.com/shahcompbio/mutationseq.git
    ```
4. create a conda environment:
    ```
    conda create --name mutationseq python=2.7 --file $HOME/museq/mutationseq/conda_packages.txt
    ```

3. activate the environment:
    ```
    conda activate mutationseq
    ```
4. Install mutationseq:

    download boost:
    ```
    wget https://sourceforge.net/projects/boost/files/boost/1.57.0/boost_1_57_0.tar.gz/download -O boost.tar.gz
    ```
    At the moment mutationseq only supports boost 1.57.0 or older. Newer versions will run into installation issues.
    extract boost:
    ```
    tar -xvf boost.tar.gz
    ```
    install mutationseq:
    ```
    cd $HOME/museq/mutationseq/
    python setup.py install --boost_source=$HOME/museq/boost_1_57_0
    ```

## Running Mutationseq:
To call variants using MutationSeq, we use the following command:
```
mkdir -p museq/results; \
museq \
  normal:bam/HCC1395_exome_normal.sort.markdup.17.7MB-8MB.bam \
  tumour:bam/HCC1395_exome_tumour.sort.markdup.17.7MB-8MB.bam \
  reference:refs/GRCh37-lite.fa \
  -o museq/results/HCC1395_exome_tumour_normal_17.vcf
```
