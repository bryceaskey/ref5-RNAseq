#!/bin/bash
#SBATCH --job-name=align-reads              # Job name
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu         # Where to send mail	
#SBATCH --account=jkim6                     # Group providing CPU and memory resources
#SBATCH --qos=jkim6                         # QOS to run job on (investment or burst)
#SBATCH --ntasks=1                          # Number of CPU cores to use
#SBATCH --mem=1gb                           # Job memory request
#SBATCH --time=24:00:00                     # Time limit hrs:min:sec (max is 744:00:00)
#SBATCH --output=align-reads_%j.log         # Standard output and error log

pwd; hostname; date

module load hisat2/2.2.1

echo 'Aligning reads to TAIR10 reference genome'

index=/blue/jkim6/share/braskey/data/TAIR10/HISAT2-index/TAIR10
reads=/blue/jkim6/share/braskey/data/ref5/trimmed/
aln=/blue/jkim6/share/braskey/data/ref5/HISAT2/
mkdir -p ${aln}

for id in H101 H102 H103 H551 H552 H553 R51 R52 R53 WT1 WT2 WT3 Y61 Y62 Y63
do
  # Apply HISAT2 to align trimmed and filtered reads to TAIR10 reference genome
  hisat2 -1 ${reads}${id}_1_trimmed.fq -2 ${reads}${id}_2_trimmed.fq \
    -x ${index} \
    -S ${aln}${id}.sam
done