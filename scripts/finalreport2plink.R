#############################################################################
#Title:    finalreport2plink
#Function: Convert illumina finalreport2plink
#Author:   Murat Guler (murat.guler@dkfz.de, muratgmbg@gmail.com)
#Date:     Dec 29th 2023
#Note:     Please let me know if you have any trouble
#############################################################################

# Read in final report file
finalReport <- fread("FinalReport.txt", sep = "\t", skip = 9)

# Modify sample IDs contains not an alphanumeric character with "x"
finalReport$`Sample ID` <- gsub("[^[:alnum:]]", "x", finalReport$`Sample ID`)

# Filter GC Score. Please change value 0.15 if required. The value is recommeded by illumina
finalReport <- finalReport %>% filter (!is.na(finalReport$`GC Score`) & finalReport$`GC Score` > 0.15)

# Download strand report file from the illumina webpage
system("wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/humanomniexpress-24/v1-4/InfiniumOmniExpress-24v1-4_A1_csv.zip")
system("unzip InfiniumOmniExpress-24v1-4_A1_csv.zip")
rsmap <- fread("InfiniumOmniExpress-24v1-4_A1.csv", 
               sep = ",", skip = 5)

# Merge finalreport and rsmap
merged_df <- merge(finalReport, rsmap, by.x = "SNP Name", by.y = "Name", all = TRUE)

# Remove unnecessery objects
rm(finalReport, rsmap)

# Remove Chr = 0 and NA and Sample ID NA
merged_df <- filter(merged_df, Chr != 0 & !is.na(Chr) & !is.na(`Sample ID`))

# Select SNP need to be flipped to + strand
# This code will create a file "flip.missnp" and it will used in the next step
merged_df %>%
  distinct(`SNP Name`, .keep_all = TRUE) %>%
  mutate(Flip = ifelse((IlmnStrand == "BOT" & RefStrand == "+") | (IlmnStrand == "TOP" & RefStrand == "-"), "TopM", "TopP")) %>%
  filter(Flip == "TopM") %>% dplyr::select(`SNP Name`) %>%
  fwrite("final2plink/flip.missnp", col.names = FALSE, sep = " ")

# Please check samples ids: unique(merged_df$`Sample ID`)
# you can change sample ids, add sex and phenotype in this step or later
# here I will put pseudo info (0, -9) to go next step

# Main script to create FAM, Lgen and MAP

# Fam file
merged_df %>%
  distinct(`Sample ID`) %>%
  mutate(FID = 0, IDF = 0, IDM = 0, sex = 0, phenotype = -9) %>%
  relocate(`Sample ID`, .after = FID) %>%
  fwrite("final2plink/Pheno01.fam", col.names = FALSE, sep = " ")

# Lgen file
merged_df %>%
  mutate(FID = 0) %>%
  select(FID, `Sample ID`, `SNP Name`, `Allele1 - Top`, `Allele2 - Top`) %>%
  write_delim("final2plink/Pheno01.lgen", col_names = F)

# Map file
merged_df %>%
  distinct(`SNP Name`, .keep_all = TRUE) %>%
  mutate(`deCODE(cM)` = 0) %>%
  dplyr::select(Chr, `SNP Name`, `deCODE(cM)`, MapInfo) %>%
  fwrite("final2plink/Pheno01.map", col.names = FALSE, sep = " ")

# Remove unnecessery objects
rm(merged_df)
gc()

# lfile to plink (plink1.9)

# lgen to ped file with PLINK 
system("cd final2plink ; plink --nonfounders --allow-no-sex --lfile Pheno01 --missing-genotype -  --keep-allele-order --output-missing-genotype 0 --recode --out Pheno02")
# output: ped, map

#convert ped and map to bed/bim/fam
system("cd final2plink ; plink --file Pheno02 --keep-allele-order --flip flip.missnp  --make-bed --out Pheno03")
system("mv final2plink/Pheno03* raw_plink")
# output: bed, bim, fam
