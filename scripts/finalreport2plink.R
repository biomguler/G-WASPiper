#############################################################################
#Title:    finalreport2plink
#Function: Convert Illumina finalreport2plink
#Author:   Murat Guler (murat.guler@dkfz.de, muratgmbg@gmail.com)
#Date:     Dec 29th 2023
#Note:     Please let me know if you have any trouble
#############################################################################
# Load necessary package(s)
library(tidyverse)
library(data.table)
library(dplyr)
library(ggplot2)

# Function to get command-line argument value by option name
get_arg_value <- function(option, default = "") {
  option_index <- match(option, commandArgs(trailingOnly = TRUE))
  if (!is.na(option_index) && option_index + 1 <= length(commandArgs(trailingOnly = TRUE))) {
    return(commandArgs(trailingOnly = TRUE)[option_index + 1])
  } else {
    return(default)
  }
}

# Retrieve command-line arguments
final_report_name <- get_arg_value("--frn", "FinalReport.txt")
final_report_path <- get_arg_value("--pfr", "data")
strand_file_name <- get_arg_value("--sfn", "strand.csv")
strand_file_path <- get_arg_value("--psf", "data")


# Read the final report file
finalReport <- fread(file.path(final_report_path, final_report_name), sep = "\t", skip = 9)

# Check if finalReport does not exist
if (!exists("finalReport")) {
  stop("Error: Final report is not found! Please check the name of the final report and path.")
}

# Check if finalReport is empty or not standart
if (ncol(finalReport) <= 4 || nrow(finalReport) <= 100) {
  stop("Error: Final report is too small! Please check the name of the final report and path.")
}

# Modify sample IDs that contain not an alphanumeric character with "x"
finalReport$`Sample ID` <- gsub("[^[:alnum:]]", "x", finalReport$`Sample ID`)

# Filter GC Score. Please change value 0.15 if required. The value is recommeded by illumina
finalReport <- finalReport %>% filter (!is.na(finalReport$`GC Score`) & finalReport$`GC Score` > 0.15)

# Check if any rows are left after filtering
if (nrow(finalReport) == 0) {
  stop("Error: No rows left after filtering by GC Score.")
}


# Read the strand file
rsmap <- fread(file.path(strand_file_path, strand_file_name), 
               sep = ",", skip = 5)

# Check if strand does not exist
if (!exists("rsmap")) {
  stop("Error: Strand file is not found! Please check the name of the file and path.")
}

# Check if strand file is empty or not standart
if (ncol(rsmap) <= 4 || nrow(rsmap) <= 100) {
  stop("Error: Strand file is too small! Please check the name of the strand file and path.")
}

# Merge finalreport and rsmap
merged_df <- merge(finalReport, rsmap, by.x = "SNP Name", by.y = "Name", all = TRUE)

# Remove unnecessary objects and free the memory
rm(finalReport, rsmap)
gc()
# Remove Chr = 0 and NA and Sample ID NA
merged_df <- filter(merged_df, Chr != 0 & !is.na(Chr) & !is.na(`Sample ID`))

# Select SNP needs to be flipped to + strand
# This code will create a file "flip.missnp" and it will used in the next step
merged_df %>%
  distinct(`SNP Name`, .keep_all = TRUE) %>%
  mutate(Flip = ifelse((IlmnStrand == "BOT" & RefStrand == "+") | (IlmnStrand == "TOP" & RefStrand == "-"), "TopM", "TopP")) %>%
  filter(Flip == "TopM") %>% dplyr::select(`SNP Name`) %>%
  fwrite("final2plink/flip.missnp", col.names = FALSE, sep = " ")

# Please check samples ids: unique(merged_df$`Sample ID`)
# You can change sample IDs, add sex and phenotype in this step or later
# Here I will put pseudo info (0, -9) to go to the next step

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

# Remove unnecessary objects
rm(merged_df)
gc()

# lfile to plink (plink1.9)

# lgen to ped file with PLINK 
system("cd final2plink ; plink --nonfounders --allow-no-sex --lfile Pheno01 --missing-genotype -  --keep-allele-order --output-missing-genotype 0 --recode --out Pheno02")


#convert ped and map to bed/bim/fam
system("cd final2plink ; plink --file Pheno02 --keep-allele-order --flip flip.missnp  --make-bed --out Pheno03")
system("mv final2plink/Pheno03* raw_plink")

