#!/bin/bash
#SBATCH --job-name=call-variants        # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu         # Where to send mail	
#SBATCH --account=jkim6                 # Group providing CPU and memory resources
#SBATCH --qos=jkim6                     # QOS to run job on (investment or burst)
#SBATCH --ntasks=1                      # Number of CPU cores to use
#SBATCH --mem=1gb                       # Job memory request
#SBATCH --time=24:00:00                 # Time limit hrs:min:sec (max is 744:00:00)
#SBATCH --output=call-variants_%j.log   # Standard output and error log

pwd; hostname; date

module load bcftools/1.12 samtools/1.12

echo "Calling variants (counting SNPs)"

aln=/blue/jkim6/share/braskey/data/ref5/HISAT2/
index=/blue/jkim6/share/braskey/data/TAIR10/TAIR10.fa
variants=/blue/jkim6/share/braskey/data/ref5/variants/
mkdir -p ${variants}

samtools faidx ${index}

for grp in H10 H55 R5 WT Y6
do
    for rep in 1 2 3
    do
        samtools fixmate -O bam ${aln}${grp}${rep}.sam - | samtools sort - -O bam > ${aln}${grp}${rep}.bam
        samtools index ${aln}${grp}${rep}.bam
        #rm ${aln}${grp}${rep}.sam
    done
    bcftools mpileup -Ob -Q 20 -f ${index} ${grp}1.bam,${grp}2.bam,${grp}3.bam > ${variants}${grp}_raw.bcf
    bcftools call -m -Ov ${variants}${grp}_raw.bcf - | bcftools filter - -Ov -i '%QUAL>20' > ${variants}${grp}.vcf
    bcftools index ${variants}${grp}.vcf
    bcftools stats ${variants}${grp}.vcf
done