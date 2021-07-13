#!/bin/bash
#SBATCH --job-name=count-expression                 # Job name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu                 # Where to send mail	
#SBATCH --account=jkim6                             # Group providing CPU and memory resources
#SBATCH --qos=jkim6                                 # QOS to run job on (investment or burst)
#SBATCH --ntasks=1                                  # Number of CPU cores to use
#SBATCH --mem=1gb                                   # Job memory request
#SBATCH --time=24:00:00                             # Time limit hrs:min:sec (max is 744:00:00)
#SBATCH --output=count-expression_%j.log            # Standard output and error log

pwd; hostname; date

module load subread/2.0.0 samtools/1.12

echo 'Counting mapped reads'

aln=/blue/jkim6/share/braskey/data/ref5/HISAT2/
index=/blue/jkim6/share/braskey/data/TAIR10/
counts=/blue/jkim6/share/braskey/data/ref5/expr_counts/
mkdir -p ${counts}
mkdir -p ${counts}just-counts

for id in H101 H102 H103 H551 H552 H553 R51 R52 R53 WT1 WT2 WT3 Y61 Y62 Y63
do
  # Use featureCounts to count reads which map to each gene
  featureCounts -t exon -g gene_id -p -s 0 -M -O --fraction \
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