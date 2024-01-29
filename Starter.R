#############################################################################
#Title:    GWAS QC Main Script
#Function: Finalreport2plink, QC, ancestry, 1000G merge, imputation
#Author:   Murat Guler (murat.guler@dkfz.de, muratgmbg@gmail.com)
#Date:     Dec 29th 2023
#Note:     Please let me know if you have any trouble
#############################################################################
# Optinal to clean your memory/workspace
rm(list=ls())
gc()
#############################################################################

# Install necessary package(s)
# Vector of package names
packages <- c("tidyverse", "data.table", "dplyr", "ggplot2")

# Install and load packages
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

#############################################################################

### Step 0: Setup

# Download scripts that needed, please update file names or URL in the scripts
system ("git clone https://github.com/biomguler/G-WASPiper.git")
system ("unzip G-WASPiper.zip")

# Create folders for outputs
system("mkdir final2plink")
system("mkdir raw_plink")
system("mkdir QCsteps")
system("mkdir 1000G")
system("mkdir results")
system("mkdir final")
system("mkdir topmed")
system("mkdir PGS")

#############################################################################

### Step 1: Final report to plink files
# Usage: Rscript finalreport2plink.R --frn <final_report_name> --pfr <final_report_path> --sfn <strand_file_name> --psf <strand_file_path>
# Rscript --no-save finalreport2plink.R --frn FinalReport.txt --pfr /path/to/final/report/ --sfn strand_file.csv --psf /path/to/strand/file/

system ("cd G-WASPiper/scripts ;  Rscript --no-save finalreport2plink.R --frn FinalReport.txt --pfr /path/to/final/report/ --sfn strand_file.csv --psf /path/to/strand/file/")

#############################################################################

### Step 2: QC
system("cd raw_plink ; plink --bfile Pheno03 --missing --freq --hardy ; cd .. ; mv raw_plink/plink* results")
 
# Generate plots and tables to visualize the genotyping rate, maf, hwe
system ("cd scripts ; Rscript --no-save basic_stat.R")

# Remove SNPs (and individuals) with high levels of missingness, low MAF and hwe
system("cd raw_plink ; plink --bfile Pheno03 --geno 0.02 --mind 0.02 --maf 0.01 --hwe include-nonctrl 5e-08 --make-bed --out Pheno04 ; cd .. ; mv raw_plink/Pheno04* QCsteps")

# Check for sex discrepancy
system("cd QCsteps ; plink --bfile Pheno04 --check-sex ; cd .. ; cp QCsteps/plink* results") 

# Create list of individuals with sex discrepancy
system("cd QCsteps; grep PROBLEM plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt")

# This command removes the list of individuals with the status “PROBLEM”.
system("cd QCsteps ; plink --bfile Pheno04 --remove sex_discrepancy.txt --make-bed --out Pheno05")


# Select autosomal SNPs only (i.e., from chromosomes 1 to 22).
system("cd QCsteps ; plink --bfile Pheno04 --autosome --make-bed --out Pheno06")


# Generate a plot of the distribution of the heterozygosity rate of your subjects.
system("plink --bfile QCsteps/Pheno06 --exclude scripts/high-LD-regions-hg19-GRCh37.txt --range --indep-pairwise 50 5 0.3 --out indepSNP ; mv indepSNP* QCsteps/") #highld region for 37 and 38 provided

# Calculate heterozygosity and remove outliers
system("cd QCsteps ; plink --bfile Pheno06 --extract indepSNP.prune.in --het --out check")

# Plot of the heterozygosity rate distribution
system("Rscript --no-save heterozygosity.R")

# Remove heterozygosity rate outliers
system("plink --bfile QCsteps/Pheno06 --remove results/fail-het.txt --make-bed --out Pheno07 ; mv Pheno07* QCsteps")

#############################################################################

### Step 4: Relatedness, PCA, ancestry inference, 1000 Genome Merge
## a) Related Individual

# Check for relationships between individuals with a pihat > 0.2.
system("cd QCsteps ; plink --bfile Pheno06 --extract indepSNP.prune.in --genome --out pihat ; cd .. ; mv QCsteps/pihat* results")

# Plot of the relatedness rate distribution
system("Rscript --no-save Relatedness.R")

# To remove this individual (if there is any) run code below 
system("plink --bfile QCsteps/Pheno06 --remove results/relatedness_fail.txt --make-bed --out Pheno08 ; mv Pheno08* QCsteps")

## b) PCA, ancestry inference, 1000 Genome Merge

# Download and proccess 1000 genome data are avaible on the plink2 webpage
system("Rscript --no-save 1000G_download.R")

# Merge your study population with 1000G
system("Rscript --no-save 1000G_merge.R")

# Run fraposa 
system("bash scripts/fraposa.sh")

# Summary of fraposa
system("Rscript --no-save scripts/fraposa.R")

# Subset sample with an ancestry (Defult EUR)
system("Rscript --no-save scripts/ancestry.R")


# Pre Topmed/michigan imputation
system("bash scripts/plink2vcf.sh")

cd topmed
# Unzip and make copy
# TopMed store the imputation results for a week
# It is good idea to make copy of the this raw files before doing any modification/filter
# Optional make copy of raw files
mkdir topmed_raw
cp *zip topmed_raw


# After you do not have need this files you can safely delete them

cd topmed
mkdir topmed_raw
cp *zip topmed_raw

# Now we can continue with files in the topmed folder
# Unzip files
# Please change "your-password" argument with password that you get from topmed after your imputation done via email
system("bash scripts/unziptopmed.sh")

#!/bin/bash
for ((chr=1; chr<=22; chr++)); do
    unzip -P 'your-password' chr_${chr}.zip
done



# download reference genome and decompress
wget -O GRCh38_full_analysis_set_plus_decoy_hla.fa.zst https://www.dropbox.com/s/xyggouv3tnamh0j/GRCh38_full_analysis_set_plus_decoy_hla.fa.zst?dl=1
zstd --decompress GRCh38_full_analysis_set_plus_decoy_hla.fa.zst


#concat combines all the chromosomes into a single file
#view filters by info score
#norm normalises indels. Split multiallelic sites into biallelic records. SNPs and indels merged into a single record
#create final gzipped VCF file and annotate. Remove original SNP ID and assign new SNP ID as chrom:position:ref:alt

bcftools concat chr*.dose.vcf.gz -Ou | 
bcftools view -Ou -i 'R2>0.3' |
bcftools norm -Ou -m -any |
bcftools norm -Ou -f hs37d5.fa |
bcftools annotate -Oz -x ID -I +'%CHROM:%POS:%REF:%ALT' -o allchr.converted.R2_0p3.vcf.gz

#--double-id means that both family and within-family IDs to be set to the sample ID

plink --vcf allchr.converted.R2_0p3.vcf.gz \
--double-id \
--allow-extra-chr 0 \
--maf 0.0000001 \
--make-bed \
--out allchr.converted.R2_0p3.MAF_1e7

#You can also include these other flags - depends on what you want to do

# If you want to filter by minimum posterior probability
--vcf-min-GP 0.9 \ 

# If you have any spaces in your IDs, it converts to _ because plink does not allow spaces in IDs
--vcf-idspace-to _ \ 

# QC of imputed data

system("bash scripts/QC_impute.R")
# info score
# 
# PRS

# GWAS







# Sources
# 1) Marees, A. T., de Kluiver, H., Stringer, S., Vorspan, F., Curis, E., Marie‐Claire, C., & Derks, E. M. (2018). A tutorial on conducting genome‐wide association studies: Quality control and statistical analysis. International journal of methods in psychiatric research, 27(2), e1608.






















