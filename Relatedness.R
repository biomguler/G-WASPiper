#############################################################################
#Title:    relatedness
#Function: relatedness detection
#Author:   Murat Guler (murat.guler@dkfz.de, muratgmbg@gmail.com)
#Date:     Dec 29th 2023
#Note:     Please let me know if you have any trouble
#############################################################################
rm(list=ls())
gc()
# Get IBD data and individual level missingness

pihat <-read.table(file="results/pihat.genome", header=TRUE)
imiss <-read.table(file="results/plink.imiss", header=TRUE)
#Identify related individuals
# IBD	 =	 1	 for	 duplicates	 or	 monozygotic twins	
# IBD	=	0.5	for	first-degree	relatives,	
# IBD	=	0.25	for	second-degree	relatives		
# IBD	=	0.125	for	third-degree	relatives	
# Genotyping	 error,	 LD	 and	 population	structure	 cause	 variation	 around	 these	theoretical	values	and	it	is	 typical	 to	remove	
# one	 individual	 from	 each	 pair	 with	 an	 IBD	 >	0.1875

# Identify related individuals
pihat <- pihat %>% mutate(relatedness = case_when(
  PI_HAT >= 0.1875 ~ "releted",
  PI_HAT < 0.1875 ~ "unrelated"))
# Detech relatedness_degree
related <- pihat %>% filter (relatedness == "related") %>% 
  mutate(relatedness_degree = case_when(
    PI_HAT >= 0.95 ~ "duplicates_monozygotic twins", 
    PI_HAT >= 0.45 & PI_HAT <= 0.55 ~ "first-degree",
    PI_HAT >= 0.20 & PI_HAT <= 0.25 ~ "second-degree", 
    PI_HAT >= 0.120 & PI_HAT <= 0.1875 ~ "third-degree"))

# combine relatedness and missingness
imiss1 <- imiss %>% mutate (F_MISS1 = F_MISS) %>% select(FID, IID, F_MISS1)
imiss2 <- imiss %>% mutate (F_MISS2 = F_MISS) %>% select(FID, IID, F_MISS2)
merged_df <- merge(pihat, imiss2, by.x = c("FID2", "IID2"), by.y = c("FID", "IID"), all = FALSE)
merged_df <- merge(merged_df, imiss1, by.x = c("FID1", "IID1"), by.y = c("FID", "IID"), all = FALSE)
related_missing <- merged_df %>% filter(relatedness == "related") %>% select(FID1, IID1, F_MISS1, FID2, IID2, F_MISS2, PI_HAT)

# Identify related and highest missingness, create plink file to remove these individuals
IBD_fail <- related_missing %>%
  mutate(DELTA = F_MISS1 - F_MISS2) %>%
  mutate(FAILED = case_when(
    DELTA == 0 ~ "First",
    DELTA < 0  ~ "Second",
    DELTA > 0  ~ "First"
  ))

# Filter rows based on "FAILED" column
IBD_fail_first <- IBD_fail %>% filter(FAILED == "First") %>% select(FID1, IID1)
IBD_fail_second <- IBD_fail %>% filter(FAILED == "Second") %>% select(FID2, IID2)

# Rename columns in both data frames
colnames(IBD_fail_first) <- c("FID", "IID")
colnames(IBD_fail_second) <- c("FID", "IID")

# Combine data frames
combined_IBD_fail <- bind_rows(IBD_fail_first, IBD_fail_second)

# Save results
fwrite(pihat, "results/relatedness.txt", sep = "\t")
fwrite(related, "results/relatedness_degree.txt", sep = "\t")
fwrite(related_missing, "results/relatedness_missing.txt", sep = "\t")
fwrite(combined_IBD_fail, "results/relatedness_fail.txt", sep = "\t", col.names=FALSE)

# Plot IBD distribution
pdf("results/relatedness.pdf")
hist(pihat$PI_HAT, breaks="sturges", main="IBD", col="blue", xlab="IBD", ylab="Frequency")
dev.off() 

rm(list=ls())
gc()
#End
