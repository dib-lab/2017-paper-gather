https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
https://conda.io/miniconda.html

conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda

conda create -n snakemake snakemake=4.0 python=3.6

Running the driver locally,
submit all other jobs:
```bash
. hpcc/submit.sh
```

Submit the driver too:
```bash
qsub hpcc/submit.qsub
```
