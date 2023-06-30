#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=48:00:00

#mamba activate snps 
#mamba activate snps

i=$1  #SRA Run ID or similar
name=$2  #The desired name of the output directory
mkdir $name  

# Prefetch command from the SRA Toolkit fetches SRA, WGS, or dbGaP data. Here it's fetching a data file with the ID provided.
prefetch $i --max-size 1000gb  

# vdb-validate command from the SRA Toolkit validates a SRA data file.
vdb-validate $i/$i.sra  

# fastq-dump command from the SRA Toolkit converts a SRA data file to a FASTQ file. Here it's also applying several additional options like filtering and clipping.
fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $i/$i.sra  

# rename.sh from BBTools renames a pair of FASTQ files removing specific SRA-related naming issues.
rename.sh in=${i}_pass_1.fastq.gz in2=${i}_pass_2.fastq.gz fixsra=t out=${name}/${name}.fq.gz  

# bbduk.sh from BBTools trims adapters and does quality filtering. It's given a list of adapter sequences, a minimum length for output reads, and a kmer size to use for trimming. IT ALSO MAXES OUT AT 50GB SEQUENCE! MODIFY AS NECESSSARY
bbduk.sh t=10 -Xmx36g overwrite=true in=${name}/${name}.fq.gz maxbasesout=50000000000 out=${name}/reads.trimmed.fastq.gz ref=~/modules/bbmap/adapters.fa minlen=25 ktrim=r k=23  
