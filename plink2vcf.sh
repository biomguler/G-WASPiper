#!/bin/bash
#############################################################################
#Title:    TopMed plink to vcf
#Function: Convert plink bfiles to vcf
#Author:   Murat Guler (murat.guler@dkfz.de, muratgmbg@gmail.com)
#Date:     Dec 29th 2023
#Note:     Please let me know if you have any trouble
#############################################################################

#Set arguments
cd QCsteps
# Download reference genome
wget -O hs37d5.fa.zst https://www.dropbox.com/s/d0grxbmwi0bdc2x/hs37d5.fa.zst?dl=1
zstd --decompress hs37d5.fa.zst

# Align reference allele
plink2 --bfile EUR_final_1k_merged --ref-from-fa --fa hs37d5.fa --make-bed --out EUR_final_1k_merged_aligned

#Create vcf files for uploading to imputation server for QC
#Note that the encoding for chromosome is e.g. chr22, not chr
for ((chr=1; chr<=22; chr++)); do
    plink --bfile EUR_final_1k_merged_aligned --keep-allele-order --chr $chr --recode vcf-iid --out tmp_chr${chr}
    vcf-sort tmp_chr${chr}.vcf | bgzip -c > chr${chr}_post_qc.vcf.gz
done

# Check Strand issues
for ((chr=1; chr<=22; chr++)); do
    bcftools +fixref chr${chr}_post_qc.vcf.gz -- -f hs37d5.fa > strand_issues${chr}.txt 2>&1
	grep -wE "SC|ST|NS" strand_issues*.txt > output_strand_issues.txt
done
#Cleanup
rm tmp_*
rm strand_issues*
