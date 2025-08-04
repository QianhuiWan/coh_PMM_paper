These files are used for environment set up for analyses in this paper:

# Environment Setup for Project

This folder includes Conda environment files used in the analyses for this paper.

## Available Environments

| Environment        | YAML File               | Purpose                                                  |
|--------------------|-------------------------|----------------------------------------------------------|
| `beers2_env`       | `beers2_env.yml`        | Used for RNAseq data simulation                          |
| `coh_PMM`          | `methylGrapher_env.yml` | Used for mapping WGBS to CHM13 graph genome              |
| `taco_env`         | `taco_env.yml`          | Used for TE identification with TEProf3 package          |
| `coh`              | `general_env.yml`       | Used for all other analyses in this paper                |


## Installation Instructions

> Requires [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io).

### Create Environment, for example:

```bash
conda env create -f beers2_env.yml
conda activate beers2_env
```

## Notes

- You can also use `mamba` to speed up environment creation, for example:

```bash
mamba env create -f beers2_env.yml
```

