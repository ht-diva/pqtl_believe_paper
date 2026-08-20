#!/bin/bash

#SBATCH --job-name plink_ds
#SBATCH --output %j_indexing_b38.log
#SBATCH --partition cpuq
#SBATCH --cpus-per-task 1
#SBATCH --mem 25G
#SBATCH --time 30-00:00:00


source /exchange/healthds/singularity_functions

geno_path="/processing_data/shared_datasets/plasma_proteome/interval/genotypes/b38_TOPMED_imputed/filtered"

# Smallest chromosome with 22G
chrom=21
beg=35467738
end=36447522

#e.g.: 9_132740402_133566789
loc="${chrom}:${beg}-${end}"
VCF="${geno_path}/chr_${chrom}_interval.vcf.gz"
ofile="chr${chrom}_${beg}_${end}.vcf.gz"
prefix="chr${chrom}_${beg}_${end}"

# Create samples list for hg38
#bcftools query -l ${geno_path}/chr_20_interval.vcf.gz > interval_hg38.samples

# Inspect missingness in genotype
#bcftools stats ${geno_path}/chr_${chrom}_interval.vcf.gz > output/chr${chrom}.stats

# Report hard-call missingness,
#  - per-sample  -> --missing sample-only
#  - per-variant -> --missing variant-only

#Claudia: At the time with Giulia de Sanctis,
# we had asked the plink guy and he told us we needed 
# to add filters --vcols=+nmissdosage,+fmissdosage

# Report dosage missingness using -> --genotyping-rate ['dosage']
# plink2 \
#     --vcf ${geno_path}/chr_${chrom}_interval.vcf.gz \
#     --missing sample-only \
#     --vcols=+nmissdosage,+fmissdosage \
#     --out output/chr_${chrom}

# Fastest way to check missing dosage
# bcftools query -f '[%DS\t]\n' ${VCF} \
# | grep -wo '\.' \
# | wc -l

# awk '
# {
#     for(i=1;i<=NF;i++)
#         if($i=="."){print "Missing"; exit}
# }
# END{if(NR>0) print "No missing"}


#bcftools query -f '[%DS\n]\n' ${VCF} | head -50000 | uniq -c

date
# Creating index file
# echo "\nIndexing VCF..."
# #tabix -p vcf ${VCF}
# bcftools index ${VCF}

# date
# sleep 60
# echo "\nSubseting VCF..."

# # Extract region from the VCF and index VCF subset
# bcftools view ${VCF} -r ${loc} -Oz -o ${ofile}

# date
# sleep 60

# echo "\nIndexing VCF subset..."
# tabix -p vcf ${ofile}

# sleep 10

plink2 \
  --vcf ${VCF} dosage=DS \
  --out ${prefix} \
  --chr ${chrom} \
  --from-bp ${beg} \
  --to-bp ${end} \
  --memory 6000 \
  --threads 3 \
  --make-pgen


date
echo "\nDone."