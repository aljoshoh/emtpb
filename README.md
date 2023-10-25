# The pharmacogenomic assessment of molecular epithelial-mesenchymal transition signatures reveals drug susceptibilities in cancer cell lines

This repository contains the analysis steps for assessing epithelial-mesenchymal transition scores as predictive biomarkers (emtpb) for drug response in cancer cell lines.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Code](#code)
- [Environment](#environment)
- [Usage](#usage)
- [Citation](#citation)

## Prerequisites

Most of the required software is available in a [docker](https://docs.docker.com/) container. Furthermore, [nextflow](https://www.nextflow.io/docs/latest/index.html) and [charliecloud](https://hpc.github.io/charliecloud/#) need to be installed for running the benchmarks.

## Code 

Pull the repository with
```bash
# Clone the repository
git clone https://github.com/mendenlab/emtpb.git
```

## Environment

The required software for running the analysis is available on a docker container. It can be either pulled from [dockerhub](https://hub.docker.com/) by typing `docker pull aljoshoh/emtpb` (recommended), or it can be rebuilt from scratch:

```bash
# Build and export container
cd emtpb/environment
docker build --no-cache -f Dockerfile -t "aljoshoh/emtpb" .
docker export $(docker create aljoshoh/emtpb) | gzip -c > metadata/emtpb.tar.gz
```

After having the container available, you can run the image by
```bash
docker run --rm -it -p 3838:3838 aljoshoh/emtpb /bin/bash
```

Within the image can run the analysis steps in `scripts/`. You can also run them within rstudio, which you can access by typing
```bash
/bin/start_rstudio 3838
```
while accessing the session in `http://localhost:3838` with `rstudio` as username and the generated password.

## Usage

The ordered analysis pipeline in `scripts/` can be run within the created environment. Make sure to adjust your system paths in `scripts/config.yaml`. While each of the R scripts only needs to be run once, the benchmarks in `scripts/02_benchmark_*` will need to be run in batch mode using nextflow (with adjusted paths). Then, un each of the two benchmarks by:

``` bash
# Run 1
cd nf
nextflow run train_chc.nf --file=scripts/02_benchmark_PANCAN.ipynb --amount=313600 --desc=exp3 -profile slurm -bg -with-trace -ansi-log false --Xms500M --Xmx2G > nf.log
```

``` bash
# Run 2
cd nf_2
nextflow run train_chc.nf --file=scripts/02_benchmark_CI.ipynb --amount=156800 --desc=exp6 -profile slurm -bg -with-trace -ansi-log false --Xms500M --Xmx2G > nf.log
```

The figures for the manuscript can be reproduced by running `paper/figures_02.Rmd`.

## Citation

This repository refers to our latest manuscript (citation pending). 