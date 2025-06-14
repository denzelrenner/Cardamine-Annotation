# GENESPACE analysis

This github page will explain how to replicate the GENESPACE analysis conducted on the diploid Cardamine amara, and two other diploid Cardamine hirsuta assemblies.

+ [Prerequisites](#prerequisites)
  - [Tool Version and Links](#tool-version-and-links)
  - [Tool Installation](#tool-installation)
  - [Data Acquisition](#data-acquisition)
 
    - [Cardamine amara (Haplome 1)](#cardamine-amara-haplome-1)
    - [Cardamine amara (Haplome 2)](#cardamine-amara-haplome-2)
    - [Cardamine hirsuta (Sanger)](#cardamine-hirsuta-sanger)
    - [Cardamine hirsuta (Max Planck)](#cardamine-hirsuta-max-planck)
    
+ [Analysis](#the-analysis)
  - [RepeatMasker/RepeatModeller](#repeatmaskerrepeatmodeller)
  - [BRAKER](#braker)
  - [GENESPACE](#genespace)
  


# Prerequisites

## <ins>Tool version and links<ins>

## <ins>Tool Installation<ins>
Yaml files

## <ins>Data Acquisition<ins>

To run genespace we need protein fastas for Cardamine amara (Haplome 1), Cardamine amara (Haplome 2), Cardamine hirsuta (Sanger), Cardamine hirsuta (Max Planck).

### <ins>Cardamine amara (Haplome 1)<ins>
Haplome 1 for Cardamine amara was created by R.

```bash
mkdir -p ~/Cardamine_Annotation_Haplomes/Haplome1/Input_Seqs
cp /path/to/haplome1/assembly ~/Cardamine_Annotation_Haplomes/Haplome1/Input_Seqs/haplome1.fa
```

### <ins>Cardamine amara (Haplome 2)<ins>
Haplome 2 for Cardamine amara was created by R.

```bash
mkdir -p ~/Cardamine_Annotation_Haplomes/Haplome2/Input_Seqs
cp /path/to/haplome2/assembly ~/Cardamine_Annotation_Haplomes/Haplome2/Input_Seqs/haplome2.fa
```

### <ins>Cardamine hirsuta (Sanger)<ins>
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

### <ins>Cardamine hirsuta (Max Planck)<ins>
The Cardamine hirsuta assembly and annotations generated by MaxPlanck can be found [here](https://chi.mpipz.mpg.de/download/annotations)

```bash
# define directory for output
OUTPUTDIR=~/Cardamine_Annotation_Haplomes/Shared_Input_Data/MaxPlanck_Hirsuta

# create output directory if it does not exist already
mkdir -p $OUTPUTDIR

# move into output directory
cd $OUTPUTDIR

# download gff to current dir with
wget -O Chirsuta.gff https://chi.mpipz.mpg.de/download/annotations/carhr38.gff

# download peptides to current dir
wget -O Chirsuta.pep.fa https://chi.mpipz.mpg.de/download/annotations/carhr38.aa.fa

# download fasta
wget -O Chirsuta.nucl.fa https://chi.mpipz.mpg.de/download/sequences/chi_v1.fa
```
# The Analysis

## <ins>RepeatMasker/RepeatModeller<ins>

```bash
#!/bin/bash

#SBATCH --job-name=running_repeatmodellermasker_all_assemblies
#SBATCH --partition=hmemq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=100g
#SBATCH --time=60:00:00
#SBATCH --output=/path/to/output/and/error/directory/%x.out
#SBATCH --error=/path/to/output/and/error/directory/%x.err

#initialise conda
source ~/anaconda3/etc/profile.d/conda.sh

# activate Rmodel/mask env
conda activate RModeller-Masker

## ENTIRE HAPLOME 1
# define directory for output and where input sequences can be found
INPUTSEQ=~/Cardamine_Annotation_Haplomes/Haplome1/Input_Seqs/haplome1.fa
OUTPUTDIR=~/Cardamine_Annotation_Haplomes/Haplome1/Output/RMasker

# create output directory if it does not exist already
mkdir -p $OUTPUTDIR

# move into output directory
cd $OUTPUTDIR

# mask first haplome
BuildDatabase -name haplome1_db $INPUTSEQ
RepeatModeler -database haplome1_db -threads 48 -LTRStruct
RepeatMasker -pa 48 -dir $OUTPUTDIR -lib haplome1_db-families.fa -xsmall $INPUTSEQ

# chnage name to soft mask
cp haplome1.fa.masked haplome1.softmasked.fa

## ENTIRE HAPLOME 2
# define directory for output and where input sequences can be found
INPUTSEQ=~/Cardamine_Annotation_Haplomes/Haplome2/Input_Seqs/haplome2.fa
OUTPUTDIR=~/Cardamine_Annotation_Haplomes/Haplome2/Output/RMasker

# create output directory if it does not exist already
mkdir -p $OUTPUTDIR

# move into output directory
cd $OUTPUTDIR

# mask first haplome
BuildDatabase -name haplome2_db $INPUTSEQ
RepeatModeler -database haplome2_db -threads 48 -LTRStruct
RepeatMasker -pa 48 -dir $OUTPUTDIR -lib haplome2_db-families.fa -xsmall $INPUTSEQ

# chnage name to soft mask
cp haplome2.fa.masked haplome2.softmasked.fa

# sanger hirsuta
# define directory for output and where input sequences can be found
INPUTSEQ=~/Cardamine_Annotation_Haplomes/Shared_Input_Data/NCBI_Data/Sanger_Hirsuta/ncbi_dataset/data/GCA_964212585.1/GCA_964212585.1_ddCarHirs1.hap1.1_genomic.fna
OUTPUTDIR=~/Cardamine_Annotation_Haplomes/Shared_Output_Data/Sanger_Hirusta/RMasker

# create output directory if it does not exist already
mkdir -p $OUTPUTDIR

# move into output directory
cd $OUTPUTDIR

# mask first haplome
BuildDatabase -name sanger_haplome1_db $INPUTSEQ
RepeatModeler -database sanger_haplome1_db -threads 48 -LTRStruct
RepeatMasker -pa 48 -dir $OUTPUTDIR -lib sanger_haplome1_db-families.fa -xsmall $INPUTSEQ

# chnage name to soft mask
cp GCA_964212585.1_ddCarHirs1.hap1.1_genomic.fna.masked sanger_hirsuta_haplome1.softmasked.fa

# deactivate conda env
conda deactivate

# get job id
echo "The Job ID for this job is: $SLURM_JOB_ID"
```
## <ins>BRAKER<ins>

We are not using RNA-seq data so BRAKER will default to running BRAKER2

```bash
#!/bin/bash

#SBATCH --job-name=running_braker_entire_assemblies
#SBATCH --partition=hmemq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=165g
#SBATCH --time=48:00:00
#SBATCH --output=/path/to/output/and/error/directory/%x.out
#SBATCH --error=/path/to/output/and/error/directory/%x.err

# load singularity
module load singularity

# set important variables
BRAKER_INSTALL_DIR=~/BrakerInstall
SHARED_DATA_DIR=~/Cardamine_Annotation_Haplomes/Shared_Input_Data

## ENTIRE HAPLOME 1
# set shared data, braker, input and output dirs
GENOMESEQDIR=~/Cardamine_Annotation_Haplomes/Haplome1/Output/RMasker
OUTPUTDIR=~/Cardamine_Annotation_Haplomes/Haplome1/Output/Braker

# move into shared data dir
cd $SHARED_DATA_DIR

# get orthodb partition for plants
if [ ! -f Viridiplantae.fa.gz ]; then
    echo "File not found! Downloading..."
    wget https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Viridiplantae.fa.gz
    gunzip Viridiplantae.fa.gz
fi

# move into directory with barker.sif file
cd $BRAKER_INSTALL_DIR

# set new variable for BRAKER SIF file location
BRAKER_SIF=$PWD/braker3.sif

# run braker3, below is the original code from the test3.sh script
# Author: Katharina J. hoff
# Contact: katharina.hoff@uni-greifswald.de
# Date: Jan 12th 2023

# Copy this script into the folder where you want to execute it, e.g.:
# singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test3.sh .
# Then run "bash test3.sh".

# Check whether braker3.sif is available

if [[ -z "${BRAKER_SIF}" ]]; then
    echo ""
    echo "Variable BRAKER_SIF is undefined."
    echo "First, build the sif-file with \"singularity build braker3.sif docker://teambraker/braker3:latest\""
    echo ""
    echo "After building, export the BRAKER_SIF environment variable on the host as follows:"
    echo ""
    echo "export BRAKER_SIF=\$PWD/braker3.sif"
    echo ""
    echo "You will have to modify the export statement if braker3.sif does not reside in \$PWD."
    echo ""
    exit 1
fi

# Check whether singularity exists

if ! command -v singularity &> /dev/null
then
    echo "Singularity could not be found."
    echo "On some HPC systems you can load it with \"module load singularity\"."
    echo "If that fails, please install singularity."
    echo "Possibly you misunderstood how to run this script. Before running it, please copy it to the directory where you want to execute it by e.g.:"
    echo "singularity exec -B \$PWD:\$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test3.sh ."
    echo "Then execute on the host with \"bash test3.sh\"".
    exit 1
fi

# remove output directory if it already exists
#    - viridiplantae_odb12
#                - brassicales_odb12

# output directory set previously is the wording directory. If already exists it is replaced, otherwise it is ren
wd=$OUTPUTDIR

if [ -d $wd ]; then
    rm -r $wd
fi

singularity exec -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl --genome=$GENOMESEQDIR/haplome1.softmasked.fa \
		--prot_seq $SHARED_DATA_DIR/Viridiplantae.fa --workingdir=${wd} \
		--species=Cardamine_Haplome1_Ver1 \
		--threads 40 --busco_lineage brassicales_odb10 &> brakerrun.log





## ENTIRE HAPLOME 2
# reset input and output dirs
GENOMESEQDIR=~/Cardamine_Annotation_Haplomes/Haplome2/Output/RMasker
OUTPUTDIR=~/Cardamine_Annotation_Haplomes/Haplome2/Output/Braker

# run braker3, below is the original code from the test3.sh script
# Author: Katharina J. hoff
# Contact: katharina.hoff@uni-greifswald.de
# Date: Jan 12th 2023

# Copy this script into the folder where you want to execute it, e.g.:
# singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test3.sh .
# Then run "bash test3.sh".

# Check whether braker3.sif is available

if [[ -z "${BRAKER_SIF}" ]]; then
    echo ""
    echo "Variable BRAKER_SIF is undefined."
    echo "First, build the sif-file with \"singularity build braker3.sif docker://teambraker/braker3:latest\""
    echo ""
    echo "After building, export the BRAKER_SIF environment variable on the host as follows:"
    echo ""
    echo "export BRAKER_SIF=\$PWD/braker3.sif"
    echo ""
    echo "You will have to modify the export statement if braker3.sif does not reside in \$PWD."
    echo ""
    exit 1
fi

# Check whether singularity exists

if ! command -v singularity &> /dev/null
then
    echo "Singularity could not be found."
    echo "On some HPC systems you can load it with \"module load singularity\"."
    echo "If that fails, please install singularity."
    echo "Possibly you misunderstood how to run this script. Before running it, please copy it to the directory where you want to execute it by e.g.:"
    echo "singularity exec -B \$PWD:\$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test3.sh ."
    echo "Then execute on the host with \"bash test3.sh\"".
    exit 1
fi

# remove output directory if it already exists
#    - viridiplantae_odb12                - brassicales_odb12

# output directory set previously is the wording directory. If already exists it is replaced, otherwise it is ren
wd=$OUTPUTDIR

if [ -d $wd ]; then
    rm -r $wd
fi

singularity exec -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl --genome=$GENOMESEQDIR/haplome2.softmasked.fa \
                --prot_seq $SHARED_DATA_DIR/Viridiplantae.fa --workingdir=${wd} \
		--species=Cardamine_Haplome2_Ver1 \
		--threads 40 --busco_lineage brassicales_odb10 &> brakerrun.log



## Sanger Hirsuta
# reset input and output dirs
GENOMESEQDIR=~/Cardamine_Annotation_Haplomes/Shared_Output_Data/Sanger_Hirusta/RMasker
OUTPUTDIR=~/Cardamine_Annotation_Haplomes/Shared_Output_Data/Sanger_Hirusta/Braker

# run braker3, below is the original code from the test3.sh script
# Author: Katharina J. hoff
# Contact: katharina.hoff@uni-greifswald.de
# Date: Jan 12th 2023

# Copy this script into the folder where you want to execute it, e.g.:
# singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test3.sh .
# Then run "bash test3.sh".

# Check whether braker3.sif is available

if [[ -z "${BRAKER_SIF}" ]]; then
    echo ""
    echo "Variable BRAKER_SIF is undefined."
    echo "First, build the sif-file with \"singularity build braker3.sif docker://teambraker/braker3:latest\""
    echo ""
    echo "After building, export the BRAKER_SIF environment variable on the host as follows:"
    echo ""
    echo "export BRAKER_SIF=\$PWD/braker3.sif"
    echo ""
    echo "You will have to modify the export statement if braker3.sif does not reside in \$PWD."
    echo ""
    exit 1
fi

# Check whether singularity exists

if ! command -v singularity &> /dev/null
then
    echo "Singularity could not be found."
    echo "On some HPC systems you can load it with \"module load singularity\"."
    echo "If that fails, please install singularity."
    echo "Possibly you misunderstood how to run this script. Before running it, please copy it to the directory where you want to execute it by e.g.:"
    echo "singularity exec -B \$PWD:\$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test3.sh ."
    echo "Then execute on the host with \"bash test3.sh\"".
    exit 1
fi

# remove output directory if it already exists
#    - viridiplantae_odb12                - brassicales_odb12

# output directory set previously is the wording directory. If already exists it is replaced, otherwise it is ren
wd=$OUTPUTDIR

if [ -d $wd ]; then
    rm -r $wd
fi

singularity exec -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl --genome=$GENOMESEQDIR/sanger_hirsuta_haplome1.softmasked.fa \
		--prot_seq $SHARED_DATA_DIR/Viridiplantae.fa --workingdir=${wd} \
		--species=Sanger_Hirsuta_Ver1 \
		--threads 40 --busco_lineage brassicales_odb10 &> brakerrun.log

# get job id
echo "The Job ID for this job is: $SLURM_JOB_ID"
```

## <ins>GENESPACE<ins>

```bash
#!/bin/bash

#SBATCH --job-name=genespace
#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=50g
#SBATCH --time=12:00:00
#SBATCH --output=/path/to/output/and/error/directory/%x.out
#SBATCH --error=/path/to/output/and/error/directory/%x.err

# conda init, init?
source ~/anaconda3/etc/profile.d/conda.sh

# load conda env
conda activate genespace4

# make required directories for genespace
OutputDir=~/Cardamine_Annotation_Haplomes/Shared_Output_Data/Genespace_Results/Genespace_Amara_Vs_Hirsuta
ScriptsDir=~/Cardamine_Annotation_Haplomes/Shared_Output_Data/Genespace_Results/Genespace_Amara_Vs_Hirsuta/Scripts
Haplome1DataDir=~/Cardamine_Annotation_Haplomes/Haplome1/Output/Braker
Haplome2DataDir=~/Cardamine_Annotation_Haplomes/Haplome2/Output/Braker
PlanckHirsutaDataDir=~/Cardamine_Annotation_Haplomes/Shared_Input_Data/MaxPlanck_Hirsuta
SangerHirsutaDataDir=~/Cardamine_Annotation_Haplomes/Shared_Output_Data/Sanger_Hirsuta/Braker

# make output dir
mkdir -p $OutputDir

# move into main output dir
# make sub directories
cd $OutputDir
mkdir -p InputData
cd InputData

# set variable for genespace input dir, do not change this
GenespaceInputDataDir=$OutputDir/InputData

# copy protein fasta files to input data directory
cp $Haplome1DataDir"/braker.aa" haplome1_braker.fa
cp $Haplome2DataDir"/braker.aa" haplome2_braker.fa
cp $PlanckHirsutaDataDir"/Chirsuta.pep.fa" cardamine_hirsuta_planck.aa
cp $SangerHirsutaDataDir"/braker.aa" cardamine_hirsuta_sanger.aa

# copy gtfs to input data directory
cp $Haplome1DataDir"/braker.gtf" haplome1_braker.gtf
cp $Haplome2DataDir"/braker.gtf" haplome2_braker.gtf
cp $PlanckHirsutaDataDir"/Chirsuta.gff" cardamine_hirsuta_planck.gtf
cp $SangerHirsutaDataDir"/braker.gtf" cardamine_hirsuta_sanger.gtf

# remove debris scaffolds and proteins found in debris from gff and protein fasta for planck assembly of Cardamine hirsuta
python3 ~/Cardamine_Annotation_Haplomes/Scripts/Python_Scripts/remove_debris_scaffolds.py \
        --nucl ~/Cardamine_Annotation_Haplomes/Shared_Input_Data/MaxPlanck_Hirsuta/Chirsuta.nucl.fa \
        --gff cardamine_hirsuta_planck.gtf --chr 8 -o cardamine_hirsuta_planck \
        -d ./ --sep _ --chr_prefix Chr --debris RL_9 --rename --prt cardamine_hirsuta_planck.aa --feature mRNA

# remove debris scaffolds and proteins found in debris from gff and protein fasta for sanger assembly of Cardamine hirsuta
python3 ~/Cardamine_Annotation_Haplomes/Scripts/Python_Scripts/remove_debris_scaffolds.py \
	--nucl ~/Cardamine_Annotation_Haplomes/Shared_Input_Data/NCBI_Data/Sanger_Hirsuta/ncbi_dataset/data/GCA_964212585.1/GCA_964212585.1_ddCarHirs1.hap1.1_genomic.fna \
	--gff cardamine_hirsuta_sanger.gtf --chr 8 -o cardamine_hirsuta_sanger \
	-d ./ --sep _ --chr_prefix Chr --debris RL_9 --rename --prt cardamine_hirsuta_sanger.aa

# replace RL_ prefix with Chr
sed -i 's/RL_/Chr/g' haplome1_braker.gtf
sed -i 's/RL_/Chr/g' haplome2_braker.gtf

# convert to bed
grep -P '\ttranscript\t' haplome1_braker.gtf | awk '{print $1,$4,$5,$9}' > haplome1_braker.bed
grep -P '\ttranscript\t' haplome2_braker.gtf | awk '{print $1,$4,$5,$9}' > haplome2_braker.bed
grep -P '\ttranscript\t' cardamine_hirsuta_sanger.gff | awk '{print $1,$4,$5,$9}' > cardamine_hirsuta_sanger.bed

# filter gff for planck to only contain mRNA lines, then convert gff to a bed file
grep -P '\tmRNA\t' cardamine_hirsuta_planck.gff > filtered_cardamine_hirsuta_planck.gff
python3 ~/Cardamine_Annotation_Haplomes/Scripts/Python_Scripts/convert_gff_to_bed.py --gff filtered_cardamine_hirsuta_planck.gff \
       -o cardamine_hirsuta_planck.bed --gene_index 0 --feature mRNA

# remove any intermediate files
rm haplome1_braker.gtf
rm haplome2_braker.gtf
rm cardamine_hirsuta_sanger.gtf
rm cardamine_hirsuta_sanger.gff
rm cardamine_hirsuta_sanger.aa
rm cardamine_hirsuta_planck.gtf
rm cardamine_hirsuta_planck.gff
rm filtered_cardamine_hirsuta_planck.gff
rm cardamine_hirsuta_planck.aa

# run python script to sort all this input data and create genespace scripts
python3 ~/Cardamine_Annotation_Haplomes/Scripts/Python_Scripts/create_riparian_plots.py -od $OutputDir \
	-id $GenespaceInputDataDir \
	--reference haplome1_braker \
	--threads 48 \
	--bg_colour 'gray90'

# move into scripts dir
cd $ScriptsDir

# run script
for file in *; do
    echo $file
    if [ -f $file ]; then
        echo "Running ${file} .."
	R -f $file
    fi
done

# deactivate env
conda deactivate

# echo job id
echo "The Job ID for this job is: $SLURM_JOB_ID"
```


