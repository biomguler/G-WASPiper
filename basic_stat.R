#############################################################################
#Title:    basic statistics
#Function: illustrate and calculate --mind --geno, hwe and maf
#Author:   Murat Guler (murat.guler@dkfz.de, muratgmbg@gmail.com)
#Date:     Dec 29th 2023
#Note:     Please let me know if you have any trouble
#############################################################################
# Get individual and SNP level missingness data
indmiss<-read.table(file="results/plink.imiss", header=TRUE)
snpmiss<-read.table(file="results/plink.lmiss", header=TRUE)

# Specify the threshold values can be used plink --geno and --mind arguments
thresholds <- c(0, 0.0001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.5, 0.75, 1)

# Create output dataframe
snpmisstable <- data.frame(thresholds = thresholds,
                        Count = sapply(thresholds, function(threshold) sum(snpmiss$F_MISS < threshold)))
indmisstable <- data.frame(thresholds = thresholds,
                        Count = sapply(thresholds, function(threshold) sum(indmiss$F_MISS < threshold)))
# Save SNP and individual level missingness results
fwrite(snpmisstable, "results/genothreshold.txt")
fwrite(indmisstable, "results/mindthreshold.txt")
# Plot genotyping rate
pdf("results/missing.pdf")
par(mfrow=c(1,2))
hist(1-indmiss$F_MISS, breaks="sturges",main="Individuals",col="blue", 
     xlab="Genotyping Rate", ylab="Number of Individuals")
hist(1-snpmiss$F_MISS, breaks="sturges", main="SNPs", col="blue", 
     xlab="Genotyping Rate", ylab="Number of SNPs")
dev.off()

rm(list=ls())
gc()
## HWE
# Get hwe results 
hwe<-read.table(file="results/plink.hwe", header=TRUE)
# Thresholds
thresholds <- c(1e-50, 1e-30, 1e-20, 1e-15, 1e-10, 5e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2)

# HWE summary table: please check plink hardy and hwe arguments
hwe_all <- hwe %>% filter(TEST == "ALL")
hwe_aff <- hwe %>% filter(TEST == "AFF")
hwe_unaff <- hwe %>% filter(TEST == "UNAFF")

hwe_all_table <- data.frame(thresholds = thresholds,
                        Count_hwe_all = sapply(thresholds, function(threshold) sum(hwe_all$P > threshold)))
hwe_aff_table <- data.frame(thresholds = thresholds,
                        Count_hwe_aff = sapply(thresholds, function(threshold) sum(hwe_aff$P > threshold)))
hwe_unaff_table <- data.frame(thresholds = thresholds,
                        Count_hwe_unaff = sapply(thresholds, function(threshold) sum(hwe_unaff$P > threshold)))
						
hwe_table <- Reduce(function(x, y) merge(x, y, by = "thresholds", all = TRUE), 
                    list(hwe_all_table, hwe_aff_table, hwe_unaff_table))
					
fwrite(hwe_table, "results/hwethreshold.txt")

# Plot HWE p values
custom_breaks <- c(5e-2, 5e-1, 1)

p1 <- ggplot(hwe, aes(x = hwe[, 9], fill = hwe[, 3])) +
  theme_classic() +
  geom_histogram(alpha = 1, color = 'gray80',
                 position = "identity", bins = 6) +
  labs(x = "P values", y = "Frequency") +
  scale_x_continuous(breaks = custom_breaks) +  # Set custom breaks on the x-axis
  scale_fill_manual(values = c("ALL" = "red", "AFF" = "blue", "UNAFF" = "green")) +
  guides(fill = guide_legend(title = "Test Methods")) + facet_grid(vars(hwe[, 3])) + coord_cartesian(expand = FALSE)

p2 <- ggplot(hwe, aes(x = hwe[, 9], fill = hwe[, 3])) +
  theme_classic() +
  geom_histogram(aes(y = after_stat(count / sum(count))), alpha = 1, color = 'gray80',
                 position = "identity", bins = 6) +
  labs(x = "P values", y = "Frequency") +
  scale_x_continuous(breaks = custom_breaks) +  # Set custom breaks on the x-axis
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("ALL" = "red", "AFF" = "blue", "UNAFF" = "green")) +
  guides(fill = guide_legend(title = "Test Methods")) + facet_grid(vars(hwe[, 3])) + coord_cartesian(expand = FALSE)

pdf("results/hwe.pdf")
par(mfrow=c(1,2))
p1
p2
dev.off()

rm(list=ls())
gc()

## MAF
maf<-read.table(file="results/plink.frq", header=TRUE)

# Define groups based on MAF
maf <- maf %>% filter(!is.na(MAF)) %>%
  mutate(type = case_when(
    MAF >= 0 & MAF < 0.001 ~ "very_rare",
    MAF >= 0.001 & MAF < 0.005 ~ "rare",
    MAF >= 0.005 & MAF < 0.05 ~ "low_freq",
    MAF >= 0.05 & MAF <= 0.5 ~ "common"))

# Summary of the variant frequency
mafthreshold <- data.frame(
  type = unique(maf$type),
  Count = sapply(unique(maf$type), function(type) sum(maf$type == type)))
mafthreshold <- mafthreshold %>% mutate(range = case_when(
type == "very_rare" ~ "0-0.002",
type == "rare" ~ "0.001-0.005",
type == "low_freq" ~ "0.005-0.05",
type == "common" ~ "0.05-0.005"))
rownames(mafthreshold) <- NULL


# Specify the threshold values can be used plink --maf arguments
thresholds <- c(0, 0.0001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.5)

# Create output dataframe
maftable <- data.frame(thresholds = thresholds,
                        Count = sapply(thresholds, function(threshold) sum(maf$MAF > threshold)))


# Save MAF results
fwrite(mafthreshold, "results/mafsummary.txt")
fwrite(maftable, "results/mafthreshold.txt")	
# Plot MAF distribition
custom_breaks <- c(1e-2, 5e-2, 4e-2, 3e-2, 2e-2, 1e-1, 2e-1, 3e-1, 4e-1, 5e-1)
# Histogram of all MAF
p1 <- ggplot(maf, aes(x = maf[, 5])) + 
  geom_histogram(alpha = 1, fill = 'blue', color = 'blue', position = "identity", binwidth = 0.01) +
  labs(x = "MAF", y = "Number of SNPs") + 
  scale_x_continuous(breaks= custom_breaks) +
  coord_cartesian(expand = FALSE) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# Histogram of MAF based on MAF category
p2 <- ggplot(maf, aes(x = maf[, 5])) +
  geom_histogram(alpha = 1, fill = 'blue', color = 'blue', position = "identity", binwidth = 0.01) +
  labs(x = "MAF", y = "Number of SNPs") + 
  scale_x_continuous(breaks= custom_breaks) +
  coord_cartesian(expand = FALSE) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(vars(maf[, 7]))
# Save plots  
pdf("results/maf.pdf")
par(mfrow=c(2,1))
p1
p2
dev.off()

rm(list=ls())
gc()
# End