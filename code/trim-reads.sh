#!/bin/bash
#SBATCH --job-name=trim-reads               # Job name
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu         # Where to send mail	
#SBATCH --account=jkim6                     # Group providing CPU and memory resources
#SBATCH --qos=jkim6                         # QOS to run job on (investment or burst)
#SBATCH --ntasks=1                          # Number of CPU cores to use
#SBATCH --mem=1gb                           # Job memory request
#SBATCH --time=24:00:00                     # Time limit hrs:min:sec (max is 744:00:00)
#SBATCH --output=trim-reads_%j.log          # Standard output and error log

pwd; hostname; date

module load adapterremoval/2.2.2

echo 'Trimming and filtering reads'

reads=/blue/jkim6/share/ref5/raw_data/
trimmed=/blue/jkim6/share/braskey/data/ref5/trimmed/

mkdir -p ${trimmed}

for id in H101 H102 H103 H551 H552 H553 R51 R52 R53 WT1 WT2 WT3 Y61 Y62 Y63
do
  # Apply AdapterRemoval to remove adapter sequences and filter low quality reads
  AdapterRemoval --file1 ${reads}${id}/${id}_1.fq.gz --file2 ${reads}${id}/${id}_2.fq.gz \
    --basename ${trimmed}${id} --output1 ${trimmed}${id}_1_trimmed.fq --output2 ${trimmed}${id}_2_trimmed.fq \
    --trimns --trimqualities --minlength 120
done