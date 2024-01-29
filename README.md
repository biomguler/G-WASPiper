<p align="center">
  <img width="600" src="Doc/oo_70943.jpg" alt="Alternative Text">
</p>

<p align="rigth">
  <small>Part of: Paukkunen J, Berg A, Soon V, Ødegaard F, Rosa P (2015) An illustrated key to the cuckoo wasps (Hymenoptera, Chrysididae) of the Nordic and Baltic countries, with the description of a new species. ZooKeys 548: 1-116. <a href="https://doi.org/10.3897/zookeys.548.6164">



# Short description and main usage

## What is G-WASPiper?
+ This repository is my collection of pipelines written in R to simplify, please cite original tools/approaches referred in the pipelines. The main idea here is putting every step in a pipe to ensure repredociblity and simplyfy all process. 
+ The pipeline designed to start from genotyping results to create ancestery / GWAS / TWAS / PRS / MR / xQTL / GWIS / FINE mapping (will be updated) analysis.
+ The main idea here creating standart pipeline for all process.
+ It is my personal reposotry to keep track all pipeline that I am using.
+ Of course anyone can used as it or with modifaction.
+ I will acknowledge any resourse/pipeline/code inclueded this codes.

 ## What is NOT G-WASPiper?
 + It is not automatic pipeline or click and run pipeline!
 + You need to modify some arguments (MAF, INFO, HWE, p value etc.) in the codes, so please be carefull before running any pipeline.
 + Some steps and pipelines are need strong computational resources and running this process with out any paralellization/optimization running pipelines as it will waste your time.
 + If you have access  to any HPC, please run this analysis in side to HPC. The pipelines not optimazated for parallel work.
 + It is not a novel package/software, published work. I will try to answer/fix any question/bug, but it would be regulary basic.
   
# Step-by-step G-WASPiper


## Explanation of the main script 
# This script create a final_PhenoXX plink files and ancestry inference results.
# After getting QCed and genetic ethnicity your target and 1000G data will be 
# merged for the assocation study. The 1000G samples (same genetic ancestry) 
# will be used as control to increase statictical power.

# The script use R, plink and python and bash commands.
# You can see these four language sign at the beging of the each code.
# Of course you need this softwares/packages installed before. 
# Because we want to keep track each step your output name will be changed in 
# the some steps like final_Pheno03,final_Pheno04 ... here Pheno is your phenotype
# or trait.
# If your plink files (bed, bim, fam) different then this format, I strongly 
# recommend to use this format.
# Note if you are running this script in LSF cluster use the bsub

# Setup your working directory
setwd("your/finalreport/path")
# Please create a new empyt folder and mv your data to this folder
# For example  "mkdir data" will create a data folder and your path will be "HOME/{user}/data}"
# Please set your working directory
setwd("HOME/{user}/{path/folder containing your final report}")
### Step 1: Final report to plink files
# Optinal- If your data in the plink format skip this step.
# If you have Finalreport file and want to convert it to plink (bed, bim,fam)
# Illumina final report to Igen files
# Note please change input name and path in the R script
# Please add phenotype and sex information by modifing R script to .fam file if it is availble
# Before running this script you need to find strand file from service provider
# Please change strand file link before running!!! 

# outputs: 
# folder name: final2plink/ flip.missnp  Pheno01.fam  Pheno01.lgen  Pheno01.map  Pheno02.log  Pheno02.map  Pheno02.nosex  Pheno02.ped
# folder name: raw_plink/ Pheno03.bed  Pheno03.bim  Pheno03.fam  Pheno03.log  Pheno03.nosex
# folder name: your main folder/ InfiniumOmniExpress-24v1-4_A1.csv  InfiniumOmniExpress-24v1-4_A1_csv.zip
### Step 2: QC
# Now, you should have Pheno03 .fam .bim .bed in the raw_plink folder
# If your data alread in the plink format you can start from this step
# Please add phenotype and sex information if you didnt put it in the first step
# The before removing any SNPs or individuals, we need to understand the data.
# The pipeline created for large scale studies, so sample size should be minimum 1000s the conduct this QC steps
# If you sample size below to this numbers, please dont use very strick tresholds to filter like MAF
# Also please change this thresholds accourding to your aim and consider to merge your data with public data if possible

# The very first step of QC is make a decision for the MAF, mind, geno, hwe

system("cd raw_plink ; plink --bfile Pheno03 --missing --freq --hardy ; cd .. ; mv raw_plink/plink* results")

# output: plink.frq  plink.hh  plink.hwe  plink.imiss  plink.lmiss  plink.log  plink.nosex
 
# Generate plots and tables to visualize the genotyping rate, maf, hwe
system ("cd scripts ; Rscript --no-save basic_stat.R")

# output(s): outputs in the results folder
# genothreshold.txt: shows how many SNP will remain after setting --geno filter with certain threshold
# mindthreshold.txt: shows how individual will remain after setting --mind filter with certain threshold
# missing.pdf: Histogram of genotyping rate both individual and SNP level
# hwethreshold.txt: Show how many SNP will remain after setting --hwe filter. all, aff(cases only),unaff(controls)
# hwe.pdf: Two pages. First page: histogram of x=P values y=SNP count, second page x=P values y=SNP percentage
# mafthreshold.txt: shows how many SNP will remain after setting --maf filter with certain threshold
# mafsummary.txt: show number of SNP in the MAF categories below:
# MAF >= 0 & MAF < 0.001 ~ "very_rare", MAF >= 0.001 & MAF < 0.005 ~ "rare", MAF >= 0.005 & MAF < 0.05 ~ "low_freq", MAF >= 0.05 & MAF <= 0.5 ~ "common  
# maf.pdf: Two pages. First page: histogram of x=MAF y=SNP count, second page x=MAF y=SNP by= categories
  

# Remove SNPs (and individuals) with high levels of missingness, low MAF and hwe
# Warning! if your sample size less thound or only a few hundered  be carefull to apply maf filter may results loss of too many variant
# By default the --hwe option in plink only filters for controls. Please check hwethreshold.
# To hwe all samples --hwe include-nonctrl can be used.
# Alternatively two step hwe filter can be conduct first for control and than all 

system("cd raw_plink ; plink --bfile Pheno03 --geno 0.02 --mind 0.02 --maf 0.01 --hwe include-nonctrl 5e-08 --make-bed --out Pheno04 ; cd .. ; mv raw_plink/Pheno04* QCsteps")

#outputs: Folder: QCsteps/ Pheno04.bed  Pheno04.bim  Pheno04.fam  Pheno04.hh  Pheno04.log

# Check for sex discrepancy.
# Subjects who were a priori determined as females must have a F value of <0.2, and subjects who were a priori determined as males must have a F value >0.8. This F value is based on the X chromosome inbreeding (homozygosity) estimate.
# Subjects who do not fulfill these requirements are flagged "PROBLEM" by PLINK.
# Be carefull and check also ycount and y-only argument of plink.
# Assing chr25 to PAR regions if not assinged correctly

system("cd QCsteps ; plink --bfile Pheno04 --check-sex ; cd .. ; cp QCsteps/plink* results") 

#outputs: Folder: results/ plink.sexcheck

# Remove individuals with sex discrepancy.
# Before removing any individual please check ycount and do LD-pruning
# If you can not solve sex-mismatch problem then remove these samples
# This command generates a list of individuals with the status “PROBLEM”.
system("cd QCsteps; grep PROBLEM plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt")

# This command removes the list of individuals with the status “PROBLEM”.
# Optional(if needed): system("cd QCsteps ; plink --bfile Pheno04 --remove sex_discrepancy.txt --make-bed --out Pheno05")


# Select autosomal SNPs only (i.e., from chromosomes 1 to 22).

system("cd QCsteps ; plink --bfile Pheno04 --autosome --make-bed --out Pheno06")


# Generate a plot of the distribution of the heterozygosity rate of your subjects.
# And remove individuals with a heterozygosity rate deviating more than 3 sd from the mean.

# Checks for heterozygosity are performed on a set of SNPs which are not highly correlated.
# Therefore, to generate a list of non-(highly)correlated SNPs, we exclude high-LD-regions.
# The parameters ‘50 5 0.3’ stand respectively for: the window size, the number of SNPs to shift 
#the window at each step, and the multiple correlation coefficient for a SNP being regressed on all
#other SNPs simultaneously.
# !!! IF your sample size not bigger than a few hundered please change parameters to less strong
# In the case of small sample size please use --read-freq to use external freq
# In the case of small sample size please do not remove any individual before merging your data with public data like 1000G; see 1000G steps

system("plink --bfile QCsteps/Pheno06 --exclude scripts/high-LD-regions-hg19-GRCh37.txt --range --indep-pairwise 50 5 0.3 --out indepSNP ; mv indepSNP* QCsteps/") #highld region for 37 and 38 provided

# Calculate heterozygosity and remove outliers
system("cd QCsteps ; plink --bfile Pheno06 --extract indepSNP.prune.in --het --out check")

#Output:QCsteps/check.het

# Plot of the heterozygosity rate distribution
system("Rscript --no-save heterozygosity.R")

# Output: results/heterozygosity.txt heterozygosity.pdf fail-het.txt
# add more detail

# Remove heterozygosity rate outliers (optinal,small sample size can skip)
# Optional(if needed): system("plink --bfile QCsteps/Pheno06 --remove results/fail-het.txt --make-bed --out Pheno07 ; mv Pheno07* QCsteps")


#############################################################################
#############################################################################
#############################################################################
### Step 4: Relatedness, PCA, ancestry inference, 1000 Genome Merge
## a) Related Individual
# It is essential to check datasets you analyse for cryptic relatedness.
# Assuming a random population sample we are going to exclude all individuals above the pihat threshold of 0.2 in this tutorial.

# Check for relationships between individuals with a pihat > 0.2.
# Recent	 shared	 ancestry	 for	 a	pair	 of	individuals	 (identity	 by	 descent,	IBD) can	 be	 estimated	 using	 genome-wide	 IBS	
# data	 using	 Plink.	 (IBD=pi_hat	 in	plink)	
# The	expectation	is	that	:	
# IBD	 =	 1	 for	 duplicates	 or	 monozygotic twins	
# IBD	=	0.5	for	first-degree	relatives,	
# IBD	=	0.25	for	second-degree	relatives		
# IBD	=	0.125	for	third-degree	relatives	
# Genotyping	 error,	 LD	 and	 population	structure	 cause	 variation	 around	 these	theoretical	values	and	it	is	 typical	 to	remove	
# one	 individual	 from	 each	 pair	 with	 an	 IBD	 >	0.1875	 (halfway	 between	 the	 expected	 IBD	for	third-	and	second-degree	relaIves).		
# For	 same	 reasons	 an	 IBD	 >	 0.98	 identifies	duplicates.

system("cd QCsteps ; plink --bfile Pheno06 --extract indepSNP.prune.in --genome --out pihat ; cd .. ; mv QCsteps/pihat* results")

# Output : results/ pihat.genome  pihat.log


# Plot of the relatedness rate distribution
# For the practical reason IBD based relatedness threshold can be vary but so thresholds below used to identify related and relatedness_degree
# PI_HAT >= 0.1875 ~ "releted",
# PI_HAT < 0.1875 ~ "unrelated"))
# relatedness_degree 
# PI_HAT >= 0.95 ~ "duplicates_monozygotic twins" 
# PI_HAT >= 0.45 & PI_HAT <= 0.55 ~ "first-degree"
# PI_HAT >= 0.20 & PI_HAT <= 0.25 ~ "second-degree"
# PI_HAT >= 0.120 & PI_HAT <= 0.1875 ~ "third-degree"

system("Rscript --no-save Relatedness.R")

# Output : results/ relatedness.txt relatedness_degree.txt relatedness_missing.txt  relatedness_fail.txt relatedness.pdf

# After checking results related individual need to be removed
# The R script will create relatedness_fail.txt which select related individuals with lowest genotyping rate
# To remove this individual (if there is any) run code below 

# Optional(if needed): system("plink --bfile QCsteps/Pheno06 --remove results/relatedness_fail.txt --make-bed --out Pheno08 ; mv Pheno08* QCsteps")


## b) PCA, ancestry inference, 1000 Genome Merge

# Download and proccess 1000 genome data are avaible on the plink2 webpage
# This script will download and convert raw files to plink
# Before runun script check filter parameters like maf (line 23)

system("Rscript --no-save 1000G_download.R")

# Outputs: Folder 1000G/ 1000g.bed  1000g.fam  1000g.bim 1000g.log 
# You can delete all other raw files (all_phase3.pgen all_phase3.pvar clean.bed clean.fam     genetic_map_b37 phase3_corrected.psam raw.bim raw.log
# all_phase3.pgen.zst all_phase3.pvar.zst clean.bim exclude.snps genetic_map_b37.zip raw.bed raw.fam)


# Merge your study population with 1000G
# Please open R script before running it and check what script doing.
# Script checking strand inconsistency and positions
# Trying to solve strand problem by flipping them
# Finally removing SNPs that are unsolved

system("Rscript --no-save 1000G_merge.R")

# Output(s) : Folder: QCsteps/ final_1k_merged.bed/bim/fam/log common_1000g.bed/bim/fam/log    common_Pheno06.bed/bim/fam/log


# Run fraposa 

system("bash scripts/fraposa.sh")

# Outputs: Folder fraposa/Pheno06_stupref.png Pheno06_stupref1.png Pheno06_stupref.popu Pheno06_stupref1.popu

# Summary of fraposa
system("Rscript --no-save scripts/fraposa.R")
#Outputs !!!!
# Subset sample with an ancestry (Defult EUR)
system("Rscript --no-save scripts/ancestry.R")


# Pre Topmed imputation
# This script checks ref allele aligment to reference (defualt hg37)
# Fix strand issues and report 
# Convert to vcf by chr
system("bash scripts/plink2vcf.sh")

#Outputs: Folder: QCsteps/vcf files for each chr to updated TOPMed imputation server
# Please download you imputation results to contineu next steps
# please download your topmed files in the topmed folder
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




























