# GENESPACE analysis

This github page will explain how to replicate the GENESPACE analysis conducted on the diploid Cardamine amara, and two other diploid Cardamine hirsuta assemblies.

+ [Prerequisites](#prerequisites)
  - [Tool Version and Links](#tool-version-and-links)
  - [Tool Installation](#tool-installation)
  - [Data Acquisition](#data-acquisition)
    	*[Cardamine hirsuta (Sanger)](#cardamine-hirsuta-sanger)
    
+ [Analysis](#the-analysis)
  


# Prerequisites

## Tool version and links

## Tool Installation
Yaml files

## Data Acquisition

To run genespace we need protein fastas for Cardamine amara (Haplome 1), Cardamine amara (Haplome 2), Cardamine hirsuta (Sanger), Cardamine hirsuta (Max Planck).

### Cardamine hirsuta (Sanger)
The Cardamine hirsuta assembly produced by Sanger is hosted on NCBI. To download the assembly run the code below:

```bash
#!/bin/bash

#SBATCH --job-name=downloading_cardamine_hirsuta
#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=24g
#SBATCH --time=01:00:00
#SBATCH --output=/path/to/Output/and/Error/dir/OnE/%x.out
#SBATCH --error=/path/to/Output/and/Error/dir/OnE/%x.err

#initialise conda
source ~/anaconda3/etc/profile.d/conda.sh

# activate conda env with ncbi datasets
conda activate ncbi_datasets_env

# define directory for output
OUTPUTDIR=~/Cardamine_Annotation_Haplomes/Shared_Input_Data/NCBI_Data/Sanger_Hirsuta

# create output directory if it does not exist already
mkdir -p $OUTPUTDIR

# move into output directory
cd $OUTPUTDIR

# download cardamine hirsuta dataset
datasets download genome accession GCA_964212585.1 \
	--include genome \
	--filename Sanger_Hirsuta_dataset.zip

# unzip
unzip Sanger_Hirsuta_dataset.zip

# deactivate conda env
conda deactivate

```

# The analysis


