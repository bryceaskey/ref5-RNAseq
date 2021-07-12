#!/bin/bash
#SBATCH --job-name=count-expression                 # Job name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=...@ufl.edu                     # Where to send mail	
#SBATCH --account=jkim6                             # Group providing CPU and memory resources
#SBATCH --qos=jkim6                                 # QOS to run job on (investment or burst)
#SBATCH --ntasks=1                                  # Number of CPU cores to use
#SBATCH --mem=1gb                                   # Job memory request
#SBATCH --time=24:00:00                             # Time limit hrs:min:sec (max is 744:00:00)
#SBATCH --output=count-expression_%j.log            # Standard output and error log

pwd; hostname; date

module load subread/2.0.0 samtools/1.10

echo 'Counting mapped reads'

aln=/ufrc/jkim6/...
index=/ufrc/jkim6/...
counts=/ufrc/jkim6/...
mkdir -p ${counts}
mkdir -p ${counts}just-counts

for id in ...
do
  # Use featureCounts to count reads which map to each gene
  featureCounts -T 1 -t exon -g gene_id -s 2 \
    -a ${index}TAIR10.gff \
    -o ${counts}${id}-counts.txt \
    ${aln}${id}.sam

  # Extract gene expression counts from output files
  cut -f 1,7 ${counts}${id}-counts.txt > ${counts}just-counts/${id}-counts.txt
  tail -n +2 ${counts}just-counts/${id}-counts.txt > ${counts}just-counts/${id}-counts.txt.tmp && mv ${counts}just-counts/${id}-counts.txt.tmp ${counts}just-counts/${id}-counts.txt

  # Extract gene lengths - needed to calculate RPKM/FPKM and TPM
  cut -f 1,6 ${counts}${id}-counts.txt > ${counts}gene-lengths.txt
  tail -n +2 ${counts}gene-lengths.txt > ${counts}gene-lengths.txt.tmp && mv ${counts}gene-lengths.txt.tmp ${counts}gene-lengths.txt
done