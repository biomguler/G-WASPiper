#############################################################################
#Title:    relatedness
#Function: relatedness detection
#Author:   Murat Guler (murat.guler@dkfz.de, muratgmbg@gmail.com)
#Date:     Dec 29th 2023
#Note:     Please let me know if you have any trouble
#############################################################################
rm(list=ls())
gc()
# Download files
system("cd 1000G ; wget https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst")
system("cd 1000G ; wget https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst")
system("cd 1000G ; wget https://www.dropbox.com/s/6ppo144ikdzery5/phase3_corrected.psam")

# Use PLINK2 to decompress the pgen and pvar files

system("cd 1000G ; plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen")
system("cd 1000G ; plink2 --zst-decompress all_phase3.pvar.zst > all_phase3.pvar")

# Use PLINK2 to convert to binary PLINK format, restricting to autsomal SNPs with MAF>0.01 (and excluding duplicates and SNPs with name ".")

system("cd 1000G ; echo '.' > exclude.snps")
system("cd 1000G ; plink2 --make-bed --out raw --pgen all_phase3.pgen --pvar all_phase3.pvar --psam phase3_corrected.psam --maf 0.01 --autosome --snps-only just-acgt --max-alleles 2 --rm-dup exclude-all --exclude exclude.snps")


# Download genetic distances, then insert these using PLINK1.9

system("cd 1000G ; wget https://www.dropbox.com/s/slchsd0uyd4hii8/genetic_map_b37.zip")
system("cd 1000G ; unzip genetic_map_b37.zip")
system("cd 1000G ; plink --bfile raw --cm-map genetic_map_b37/genetic_map_chr@_combined_b37.txt --make-bed --out 1000g")

rm(list=ls())
gc()
#End
