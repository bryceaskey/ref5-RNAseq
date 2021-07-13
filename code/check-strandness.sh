#!/bin/bash
#SBATCH --job-name=check-strandness                     # Job name
#SBATCH --mail-type=END,FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu                     # Where to send mail	
#SBATCH --account=jkim6                                 # Group providing CPU and memory resources
#SBATCH --qos=jkim6                                     # QOS to run job on (investment or burst)
#SBATCH --ntasks=1                                      # Number of CPU cores to use
#SBATCH --mem=4gb                                       # Job memory request
#SBATCH --time=24:00:00                                 # Time limit hrs:min:sec (max is 744:00:00)
#SBATCH --output=check-strandness_%j.log                # Standard output and error log

pwd; hostname; date

module load salmon/1.1

echo 'Checking strandness (e.g. unstranded, FR, RF) of PRJNA388948 RNAseq reads'

reads=/blue/jkim6/share/braskey/data/ref5/trimmed/
index=/blue/jkim6/share/braskey/data/TAIR10/salmon-index/
output=/blue/jkim6/share/braskey/data/ref5/salmon/H101/
mkdir -p ${output}

# Build index from TAIR10 genome and AtRTD2 transcriptome (https://ics.hutton.ac.uk/atRTD/RTD2/AtRTD2_19April2016.fa)
#grep "^>" < ${index}TAIR10.fa | cut -d " " -f 1 > ${index}decoys.txt
#sed -i.bak -e 's/>//g' ${index}decoys.txt
#cat ${index}AtRTD2.fa ${index}TAIR10.fa > ${index}AtRTD2_TAIR10_merged.fa
#salmon index -p 1 -t ${index}AtRTD2_TAIR10_merged.fa -d ${index}decoys.txt -i ${index}

# Quantify reads against index, and infer library type
salmon quant -p 1 -l A -i ${index} -o ${output} \
  -1 ${reads}H101_1_trimmed.fq -2 ${reads}H101_2_trimmed.fq