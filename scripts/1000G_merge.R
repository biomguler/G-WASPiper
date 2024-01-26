#############################################################################
#Title:    relatedness
#Function: relatedness detection
#Author:   Murat Guler (murat.guler@dkfz.de, muratgmbg@gmail.com)
#Date:     Dec 29th 2023
#Note:     Please let me know if you have any trouble
#############################################################################
rm(list=ls())
gc()
#
# Test strand inconsistency and positions
# In the first step SNPs are in the same strand and position
# If you see "Warning: Multiple positions seen for variant 'rsXXXXX'"
# or Error: YY variants with 3+ alleles present.
# Here important part in how big YY, if it is thounds probably it is due to strand inconsistency, and we need to do --flip with trial1-merge.missnp.

# First test strand problems
system("plink --bfile QCsteps/Pheno06 --bmerge 1000G/1000g --make-bed --out trial1 ; mv trial1* QCsteps")

# Try to solve strand problem with flipping
system("cd QCsteps ; plink --bfile Pheno06 --flip trial1-merge.missnp --make-bed --out Pheno06_flip")

# second test for strand issues
system("plink --bfile QCsteps/Pheno06_flip --bmerge 1000G/1000g --make-bed --out trial2 ; mv trial2* QCsteps")

# Obtain list of variants with multiple positions/chromosomes
system("cd QCsteps ; grep -Eo 'rs[0-9]{1,}' trial2.log | sort | uniq > multpos_var.txt")

# Combine problematic SNPs
system("cd QCsteps ; cat trial2-merge.missnp multpos_var.txt > problematic_snps.txt")

# Exclude problematic SNPs 
system("cd QCsteps ; plink --bfile Pheno06_flip --exclude problematic_snps.txt --make-bed --out Pheno06_flip_tmp")

system("plink --bfile 1000G/1000g --exclude QCsteps/problematic_snps.txt --make-bed --out 1000g_tmp ; mv 1000g_tmp* QCsteps")
# Extracting data for variants common in both file sets

# Extract SNP IDs from both .bim files and sort them numerically.
system("cd QCsteps ; awk '{print $2}' 1000g_tmp.bim | sort > 1000g_tmp_sorted.txt")

system("cd QCsteps ; awk '{print $2}' Pheno06_flip_tmp.bim | sort > Pheno06_flip_tmp_sorted.txt")

# Find and extract the IDs common in both files
system("cd QCsteps ; comm -12 Pheno06_flip_tmp_sorted.txt 1000g_tmp_sorted.txt > intersecting_snps.txt")

# Check the number of overlapping SNPs
system("cd QCsteps ; wc -l intersecting_snps.txt > number_of_overlap.txt")

# Extract data using intersecting_snps.txt from origninal file sets
system("cd QCsteps ; plink --bfile Pheno06_flip_tmp --extract intersecting_snps.txt --make-bed --out common_Pheno06")
system("cd QCsteps ; plink --bfile 1000g_tmp --extract intersecting_snps.txt --make-bed --out common_1000g")

# Merge the file sets to obtain a list of variants with multiple positions/chromosomes and multiple alleles
system("cd QCsteps ; plink --bfile common_Pheno06 --bmerge common_1000g --make-bed --out final_1k_merged")


# Identify  high quality indepentent variants
system("plink --bfile QCsteps/final_1k_merged --exclude scripts/high-LD-regions-hg19-GRCh37.txt --range --indep-pairwise 50 5 0.2 --out merged_indepSNP ; mv merged_indepSNP* QCsteps/") 

#Extract indepent high QC variants

system("cd QCsteps ; plink --bfile common_Pheno06 --extract merged_indepSNP.prune.in --make-bed --out in_common_Pheno06")
system("cd QCsteps ; plink --bfile common_1000g --extract merged_indepSNP.prune.in --make-bed --out in_common_1000g")

rm(list=ls())
gc()
#End
