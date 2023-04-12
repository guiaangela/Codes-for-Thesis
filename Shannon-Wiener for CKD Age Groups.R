# Shannon-Wiener Diversity Index (from https://rpubs.com/lconteville/713954 compiled by Liliane Conteville and https://github.com/jerichocarena/Thesis-Codes/blob/main/Alpha%20diversity%20(+boxplot%20generation compiled by Jericho Carena)
# This code compiles code for Shannon-Wiener Diversity Index dotplot and boxplots
# It also compiles statistical analysis tools for normal distributed data like ANOVA and its post hoc test, Tukey's Honestly Significant Difference
# It also compiles statistical analysis tools for non-normal distributed data like Kruskal-Wallis Test and its post hoc test, Wilcoxon Test

# Prepare three .csv files
# First file should contain the total OTUs on the Y-axis and samples on the X-axis with the raw counts of OTUs
# Second file should contain the samples on the first row and age group on the second row
# Make sure that the sample names used are in the same order as in file 1
# Third file should contain the total OTUs that is in the same order as in file 1
# Load the required packages below

library("phyloseq")
library("ggplot2")
library("dplyr")
library("ggpubr")

# Load the files: raw counts table (file 1), metadata (file 2), taxonomy table (file 3)
# Settings for file 1: encoding - automatic, heading - yes, row names - use first column, separator - comma, decimal - period, quote - double, comment - none, na.strings - NA
# Settings for file 2: encoding - automatic, heading - yes, row names - use first column, separator - comma, decimal - period, quote - double, comment - none, na.strings - NA
# Settings for file 3: encoding - automatic, heading - yes, row names - automatic, separator - comma, decimal - period, quote - double, comment - none, na.strings - NA

# In this example, the RawCounts_AD is the raw counts table, the OTUs_AD is the taxonomy table, and the MetadataAge_AD is the metadata file.

head(RawCounts_AD)
head(OTUs_AD)
head(MetadataAge_AD)

setdiff(rownames(RawCounts_AD),OTUs_AD$Species)
# Species is the column name of OTUs_AD
## character(0)

rownames(OTUs_AD) <- OTUs_AD$Species
OTUs_AD =  as.matrix(OTUs_AD)


OTU = otu_table(RawCounts_AD, taxa_are_rows = TRUE)
TAX = tax_table(OTUs_AD)
sampledata = sample_data(MetadataAge_AD)
physeq1 = phyloseq(OTU,TAX,sampledata)

physeq1

sample_data(physeq1)$Age <- factor((sample_data(physeq1)$Age), levels=c("49 below","50-60","61-70","71-80","81 above","N/A"))
# The levels are the age groups of MetadataAge_AD

# dotplot generation
richness <- estimate_richness(physeq1)
# This will result to error, therefore copy the code in line 48 to the console and change it to richness <- estimate_richness(physeq1,measures="Shannon")
head(richness)
plot_richness(physeq1)
# This will also result to error, therefore copy the code in line 51 and change it to plot_richness(physeq1,measures="Shannon")

# boxplot generation
plot_richness(physeq1, x="Age", measures="Shannon", color = "Age")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

# A histogram can be created to check if the Shannon values estimated for the metagenomes are normally distributed:
hist(richness$Shannon, main="Shannon index", xlab="")

# If it is normally distributed, we can run anova tests and check if the variable Age impacts the Shannon diversity:
anova.sh = aov(richness$Shannon ~ sample_data(physeq1)$Age)
summary(anova.sh)
TukeyHSD(anova.sh)

# For post hoc test, generate the TukeyHSD in a graph:
tukey.plot.test <- TukeyHSD(anova.sh)
plot(tukey.plot.test, las = 1)

# For Non-normally distributed data, we can use Kruskal-Wallis Rank Sum Test:
kruskal.test(richness$Shannon ~ sample_data(physeq1)$Age)
# We can also get a list with the p-values resulted of the Wilcoxon Tests considering each pair of groups:
pairwise.wilcox.test(richness$Shannon, sample_data(physeq1)$Age, p.adj = "bonf")

# Modified and annotated by Guia Angela J. Castillo