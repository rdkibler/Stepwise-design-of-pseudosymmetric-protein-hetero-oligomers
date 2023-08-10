# Stepwise design of pseudosymmetric protein hetero-oligomers


You must first acquire a license to use [Rosetta](https://els2.comotion.uw.edu/product/rosetta) and [pyrosetta](https://els2.comotion.uw.edu/product/pyrosetta). 

The python environment used to run the scripts and ipython notebooks is specified in the file pyro.yml. To use this environment, you must have [conda](https://docs.conda.io/en/latest/) installed. Build the conda environment used by the pipeline by running the following command in the root directory of this repository:
```conda env create -f pyro.yml```

Install [DeepAccNet](https://github.com/hiranumn/DeepAccNet)


To run the pipeline, pick either BGL or RTR as the example path to follow. Under both paths, the pipeline is run by visiting, adjusting, and running the scripts in the numbered subdirs. In most cases, the `gentasks.py` or `gentasks.ipynb` files both provide examples of how to run the scripts as well as prepare a list of tasks to be run. Like any software downloaded from the internet, the scripts may need to be adjusted to run on your system. That may includes paths to software, paths to input files, and paths to output files.


Note that several scripts depend on `ss_grouped_vall_helix_shortLoop.h5`. This file is not included in this repository due to its large size and is available upon request. 



