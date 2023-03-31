# SPARKS

Snakemake workflow for SPARKS

## Dependencies

[Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) must be installed first. Then [./install.sh](install.sh) can install other dependencies using conda

## Install

[./install.sh](install.sh)

* Creates a conda environment with the required dependencies.
* Installs the SPARKS R package
* Sets some absolute file paths in [snakemake_config.yaml](snakemake_config.yaml).

## Run

* set input files in [snakemake_config.yaml](snakemake_config.yaml)
* [snakemake_profile/](snakemake_profile/) is used to submit jobs to the compute cluster
* [./run](./run)

### Output

* The main output file is `{cancer_type}.SPARKS.rds`
* Other output files are in:
  + bam/
  + CLIP/
  + expression/
  + rmats/
  + rMATS_matrix/
  + references/
* Log files are generated in the same directories as output files

## TODO

* Requires existing files:
  + `/mnt/isilon/xing_lab/aspera/yangt3/ENCODE/bed`
  + `/mnt/isilon/xing_lab/aspera/yangt3/ENCODE/SPARKS/ENCODE_library.KD_and_KO.count_20.rds`
