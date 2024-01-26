#############################################################################
#Title:    ancestry inference with fraposa_summary
#Function: summary of fraposa
#Author:   Murat Guler (murat.guler@dkfz.de, muratgmbg@gmail.com)
#Date:     Dec 29th 2023
#Note:     Please let me know if you have any trouble
#############################################################################
rm(list=ls())
gc()
# Read results of fraposa
super_pop <-read.table(file="fraposa/Pheno06_stupref.popu", header=FALSE)
sub_pop <-read.table(file="fraposa/Pheno06_stupref1.popu", header=FALSE)

# Fix column names and remove extra
colnames(super_pop) <- c("FID", "IID", "Main_pop", "Probality", "Distance", super_pop[1,11], super_pop[1,12], super_pop[1,13], super_pop[1,14], super_pop[1,15])
super_pop <- super_pop[,-11:-15]

colnames(sub_pop) <- c("FID", "IID", "Main_pop", "Probality", "Distance", c(sub_pop[1,32:57]))
sub_pop <- sub_pop[,-32:-57]

# These populations have been divided into 5 super populations
# 
# Non_eur
non_EUR <- super_pop %>% filter(Main_pop != "EUR")

# Non_AFR
non_AFR <- super_pop %>% filter(Main_pop != "AFR")


# Non_AMR
non_AMR <- super_pop %>% filter(Main_pop != "AMR")

# Non_EAS
non_EAS <- super_pop %>% filter(Main_pop != "EAS")

# Non_SAS
non_SAS <- super_pop %>% filter(Main_pop != "SAS")

# Summary of results

counts1 <- data.frame(table(super_pop$Main_pop))
counts2 <- data.frame(table(sub_pop$Main_pop))

# Create plots
# Super population bar chart
p1 <- ggplot(super_pop, aes(x = Main_pop)) +
  geom_bar(fill = 'blue', color = 'blue', position = "identity") +
  geom_text(data = counts1, aes(x = Var1, y = Freq, label = Freq),
            vjust = -0.5, color = "black", size = 3) +  # Add labels at the top of each bar
  labs(title = "Infered Ancestry-Super Populations",
       x = "Super Populations",
       y = "Number of Individual") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5))

# Sub population bar chart
p2 <- ggplot(sub_pop, aes(x = Main_pop)) +
  geom_bar(fill = 'blue', color = 'blue', position = "identity") +
  geom_text(data = counts2, aes(x = Var1, y = Freq, label = Freq),
            vjust = -0.5, color = "black", size = 3) +  # Add labels at the top of each bar
  labs(title = "Infered Ancestry-Sub Populations",
       x = "Sub Populations",
       y = "Number of Individual") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5))

# Save plots
pdf("results/fraposa.pdf")
par(mfrow=c(2,1))
p1
p2
dev.off()

# Save results
# Main results editted
fwrite(super_pop, "results/super_pop.txt", sep = "\t")
fwrite(sub_pop, "results/sub_pop.txt", sep = "\t")
# Summary results
colnames(counts1) <- c("Super_populations", "Count")
colnames(counts2) <- c("Sub_populations", "Count")
fwrite(counts1, "results/summary_super_pop.txt", sep = "\t")
fwrite(counts2, "results/summary_sub_pop.txt", sep = "\t")
# plink related
# Write non_EUR data frame to a tab-separated text file
fwrite(non_EUR[, 1:2], "results/non_EUR.txt", sep = "\t", col.names = FALSE)

# Write non_AFR data frame to a tab-separated text file
fwrite(non_AFR[, 1:2], "results/non_AFR.txt", sep = "\t", col.names = FALSE)

# Write non_EAS data frame to a tab-separated text file
fwrite(non_EAS[, 1:2], "results/non_EAS.txt", sep = "\t", col.names = FALSE)

# Write non_AMR data frame to a tab-separated text file
fwrite(non_AMR[, 1:2], "results/non_AMR.txt", sep = "\t", col.names = FALSE)

# Write non_SAS data frame to a tab-separated text file
fwrite(non_SAS[, 1:2], "results/non_SAS.txt", sep = "\t", col.names = FALSE)
rm(list=ls())
gc()
#End

