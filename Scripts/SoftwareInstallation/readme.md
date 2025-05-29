# Software Installation Instructions

+ [NCBI-datasets](#ncbidatasets)
+ [RepeatMasker/RepeatModeller](#repeatmaskerrepeatmodeller)
+ [BRAKER](#braker)
+ [GENESPACE](#genespace)


## NCBI-datasets

```bash
#!/bin/bash

#SBATCH --job-name=installing_ncbi_datasets
#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=5g
#SBATCH --time=1:00:00
#SBATCH --output=/path/to/OnE/directory/%x.out
#SBATCH --error=/path/to/OnE/directory/%x.err

# create a new conda environment for samtools
conda create -y -n ncbi_datasets_env

#initialise conda
source ~/anaconda3/etc/profile.d/conda.sh

# activate environment with the conda packages
conda activate ncbi_datasets_env

# install samtools
conda install -y -c conda-forge ncbi-datasets-cli

# deactivate conda env
conda deactivate

# get job id
echo "The Job ID for this job is: $SLURM_JOB_ID"

```

## RepeatMasker/RepeatModeller

```bash
#!/bin/bash

#SBATCH --job-name=installing_repeat_masker_dependencies
#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=5g
#SBATCH --time=00:20:00
#SBATCH --output=/path/to/OnE/directory/%x.out
#SBATCH --error=/path/to/OnE/directory/%x.err

# create a new conda environment for quast
conda create -y -n RMaskerModeller_env python=3.9

#initialise conda
source ~/anaconda3/etc/profile.d/conda.sh

# activate environment with the conda packages
conda activate RMaskerModeller_env

# install h5py
pip install h5py

# deactivate conda env
conda deactivate

# get job id
echo "The Job ID for this job is: $SLURM_JOB_ID"
```

## BRAKER
Instructions to install braker

## GENESPACE

```bash
#!/bin/bash

#SBATCH --job-name=installing_R
#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=40g
#SBATCH --time=01:30:00
#SBATCH --output=/path/to/OnE/directory/%x.out
#SBATCH --error=/path/to/OnE/directory/%x.err

# initialise conda
source ~/anaconda3/etc/profile.d/conda.sh

# install and load software
source $HOME/.bash_profile

conda config --remove channels defaults

conda create --name genespace4 orthofinder=2.5.5 mcscanx r-base=4.4.1 r-devtools r-BiocManager bioconductor-biostrings -y
conda activate genespace4
R
devtools::install_github("jtlovell/GENESPACE")
BiocManager::install("rtracklayer")
library(GENESPACE)

# deactivate conda
conda deactivate

# get job id
echo "The Job ID for this job is: $SLURM_JOB_ID"
```
