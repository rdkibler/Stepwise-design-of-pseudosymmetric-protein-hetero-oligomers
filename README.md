# Stepwise design of pseudosymmetric protein hetero-oligomers

This is the repo for the preprint https://www.biorxiv.org/content/10.1101/2023.04.07.535760v1

## System requirements
* RAM: 8GB
* OS: Ubuntu 22.04.2 LTS (Jammy Jellyfish) or similar version of linux

## Installation

1. Acquire a license to use [Rosetta](https://els2.comotion.uw.edu/product/rosetta) and [pyrosetta](https://els2.comotion.uw.edu/product/pyrosetta). 
2. Install [DeepAccNet](https://github.com/hiranumn/DeepAccNet)
3. Install [conda](https://docs.conda.io/en/latest/) 
4. Build the conda environment used by the pipeline by running the following command in the root directory of this repository: 
```conda env create -f pyro.yml```
This could take anywhere from a few minutes to a few hours depending on internet connection. 

Note that several scripts depend on `ss_grouped_vall_helix_shortLoop.h5`. This file is not included in this repository due to its large size and is available upon request. In many cases, the steps requiring this file can be replaced with more modern and objectively better methods, like ProteinMPNN and RFDiffusion.

## Execution

To run the pipeline, pick either BGL or RTR as the example path to follow. Under both paths, the pipeline is run by visiting, adjusting, and running the scripts in the numbered subdirs. In most cases, the `gentasks.py` or `gentasks.ipynb` files both provide examples of how to run the scripts as well as prepare a list of tasks to be run. Like any software downloaded from the internet, the scripts may need to be adjusted to run on your system. That may includes paths to software, paths to input files, and paths to output files.

Expected runtime for the complete pipeline is difficult to estimate, but expect to spend about a week of person-time and several years of CPU time. 
