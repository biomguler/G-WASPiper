#############################################################################
#Title:    heterozygosity
#Function: illustrate and calculate --mind --geno, hwe and maf
#Author:   Murat Guler (murat.guler@dkfz.de, muratgmbg@gmail.com)
#Date:     Dec 29th 2023
#Note:     Please let me know if you have any trouble
#############################################################################

# Get heterozygosity output
het <- read.table("QCsteps/check.het", head=TRUE)

# Calculate heterozygosity rate and identify individual HET_RATE > +-3 SD
het <- het %>% 
  mutate(HET_RATE = (N.NM. - O.HOM.) / N.NM.) %>% 
  mutate(HET_FAIL = ifelse(HET_RATE < mean(HET_RATE) - 3 * sd(HET_RATE) | HET_RATE > mean(HET_RATE) + 3 * sd(HET_RATE), "FAIL", "NORMAL"))

# Save Results
fwrite(het, "results/heterozygosity.txt", sep = "\t")

# Create list of FAILED HET_QC individual
system("cd results; awk 'NR > 1' heterozygosity.txt | grep FAIL | awk '{print $1, $2}'> fail-het.txt")

# Plot heterozygosity rate

pdf("results/heterozygosity.pdf")
hist(het$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate", col="blue")
dev.off()

rm(list=ls())
gc()
#End