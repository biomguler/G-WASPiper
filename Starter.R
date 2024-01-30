#############################################################################
#Title:    G-WASPiper-Starter.R
#Function: Install packages and create folders
#Author:   Murat Guler (murat.guler@dkfz.de, muratgmbg@gmail.com)
#Date:     Dec 29th 2023
#Note:     Please let me know if you have any trouble
#############################################################################
# Optinal to clean your memory/workspace
rm(list=ls())
gc()
#############################################################################
# Install and check necessary package(s)
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
# Download scripts that needed, please update file names or URL in the scripts
system ("git clone https://github.com/biomguler/G-WASPiper.git") #update if needed
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
# Optinal to clean your memory/workspace
rm(list=ls())
gc()
#End
