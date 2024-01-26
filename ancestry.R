#############################################################################
#Title:    Ancestry
#Function: Remove samples from other ancestry
#Author:   Murat Guler (murat.guler@dkfz.de, muratgmbg@gmail.com)
#Date:     Dec 29th 2023
#Note:     Please let me know if you have any trouble
#############################################################################
rm(list=ls())
gc()

# Remove all non_EUR or other ancestry using output of previous step
system("plink --bfile QCsteps/Pheno06 --remove results/non_EUR.txt --make-bed --out Pheno07 ; mv Pheno07* QCsteps")

# Select and combine non_EUR samples from 1000G
system("cat scripts/1000g_notEUR.txt results/non_EUR.txt > no_EUR.txt ; mv no_EUR.txt results")
# Remove no_EUR samples and create plink files
system("plink --bfile QCsteps/final_1k_merged --remove results/no_EUR.txt --make-bed --out EUR_final_1k_merged ; mv EUR_final_1k_merged* QCsteps")

# Assing control to 1000G EUR samples
system("cd QCsteps ; awk '{$6 = ($6 == -9) ? 1 : $6; print}' EUR_final_1k_merged.fam > modified_EUR_final_1k_merged.fam")
system("cd QCsteps ; mv modified_EUR_final_1k_merged.fam EUR_final_1k_merged.fam")

rm(list=ls())
gc()
#End