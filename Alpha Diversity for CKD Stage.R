# Shannon-Wiener Diversity Index (from https://rpubs.com/lconteville/713954 compiled by Liliane Conteville and https://github.com/jerichocarena/Thesis-Codes/blob/main/Alpha%20diversity%20(+boxplot%20generation compiled by Jericho Carena)
# This code compiles code for Shannon-Wiener and Simpson Diversity and ACE and Observed indices' dotplots and boxplots
# It also compiles statistical analysis tools for normal distributed data like ANOVA and its post hoc test, Tukey's Honestly Significant Difference
# It also compiles statistical analysis tools for non-normal distributed data like Kruskal-Wallis Test and its post hoc test, Wilcoxon Test

# Prepare three .csv files
# First file should contain the total OTUs on the Y-axis and samples on the X-axis with the raw counts of OTUs
# Second file should contain the samples on the first row and clinical status on the second row
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

# In this example, the RawCounts_AD is the raw counts table, the OTUs_AD is the taxonomy table, and the MetadataStage_AD is the metadata file.

head(RawCounts_AD)
head(OTUs_AD)
head(MetadataStage_AD)

setdiff(rownames(RawCounts_AD),OTUs_AD$Species)
# Species is the column name of OTUs_AD
## character(0)

rownames(OTUs_AD) <- OTUs_AD$Species
OTUs_AD =  as.matrix(OTUs_AD)


OTU = otu_table(RawCounts_AD, taxa_are_rows = TRUE)
TAX = tax_table(OTUs_AD)
sampledata = sample_data(MetadataStage_AD)
physeq1 = phyloseq(OTU,TAX,sampledata)

physeq1

sample_data(physeq1)$Stage <- factor((sample_data(physeq1)$Stage), levels=c("Healthy","Stage 3","Stage 4","Stage 5","Non-ESRD","Generally Afflicted"))
# The levels are the clinical status of MetadataStage_AD

# dotplot generation of Shannon Diversity Index
richness1 <- estimate_richness(physeq1,measures="Shannon")
head(richness1)
plot_richness(physeq1,measures="Shannon")

# boxplot generation of Shannon Diversity Index
plot_richness(physeq1, x="Stage", measures="Shannon", color = "Stage")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

# A histogram can be created to check if the Shannon values estimated for the metagenomes are normally distributed:
hist(richness1$Shannon, main="Shannon index", xlab="")

# If it is normally distributed, we can run anova tests and check if the variable Stage impacts the Shannon diversity:
anova.sh = aov(richness1$Shannon ~ sample_data(physeq1)$Stage)
summary(anova.sh)
TukeyHSD(anova.sh)

# For post hoc test, generate the TukeyHSD in a graph:
tukey.plot.test <- TukeyHSD(anova.sh)
plot(tukey.plot.test, las = 1)

# For Non-normally distributed data, we can use Kruskal-Wallis Rank Sum Test:
kruskal.test(richness1$Shannon ~ sample_data(physeq1)$Stage)
# We can also get a list with the p-values resulted of the Wilcoxon Tests considering each pair of groups:
pairwise.wilcox.test(richness1$Shannon, sample_data(physeq1)$Stage, p.adj = "bonf")

# dotplot generation of Simpson Diversity Index
richness2 <- estimate_richness(physeq1,measures="Simpson")
head(richness2)
plot_richness(physeq1,measures="Simpson")

# boxplot generation of Simpson Diversity Index
plot_richness(physeq1, x="Stage", measures="Simpson", color = "Stage")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

# A histogram can be created to check if the Simpson values estimated for the metagenomes are normally distributed:
hist(richness2$Simpson, main="Simpson index", xlab="")

# If it is normally distributed, we can run anova tests and check if the variable Stage impacts the Simpson diversity:
anova.sh = aov(richness2$Simpson ~ sample_data(physeq1)$Stage)
summary(anova.sh)
TukeyHSD(anova.sh)

# For post hoc test, generate the TukeyHSD in a graph:
tukey.plot.test <- TukeyHSD(anova.sh)
plot(tukey.plot.test, las = 1)

# For Non-normally distributed data, we can use Kruskal-Wallis Rank Sum Test:
kruskal.test(richness2$Simpson ~ sample_data(physeq1)$Stage)
# We can also get a list with the p-values resulted of the Wilcoxon Tests considering each pair of groups:
pairwise.wilcox.test(richness2$Simpson, sample_data(physeq1)$Stage, p.adj = "bonf")

# dotplot generation of ACE Index
richness3 <- estimate_richness(physeq1,measures="ACE")
head(richness3)
plot_richness(physeq1,measures="ACE")

# boxplot generation of ACE Index
plot_richness(physeq1, x="Stage", measures="ACE", color = "Stage")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

# A histogram can be created to check if the Simpson values estimated for the metagenomes are normally distributed:
hist(richness3$ACE, main="ACE index", xlab="")

# If it is normally distributed, we can run anova tests and check if the variable Stage impacts the ACE index:
anova.sh = aov(richness3$ACE ~ sample_data(physeq1)$Stage)
summary(anova.sh)
TukeyHSD(anova.sh)

# For post hoc test, generate the TukeyHSD in a graph:
tukey.plot.test <- TukeyHSD(anova.sh)
plot(tukey.plot.test, las = 1)

# For Non-normally distributed data, we can use Kruskal-Wallis Rank Sum Test:
kruskal.test(richness3$ACE ~ sample_data(physeq1)$Stage)
# We can also get a list with the p-values resulted of the Wilcoxon Tests considering each pair of groups:
pairwise.wilcox.test(richness3$ACE, sample_data(physeq1)$Stage, p.adj = "bonf")

# dotplot generation of Observed Index
richness4 <- estimate_richness(physeq1,measures="Observed")
head(richness4)
plot_richness(physeq1,measures="Observed")

# boxplot generation of Observed Index
plot_richness(physeq1, x="Stage", measures="Observed", color = "Stage")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

# A histogram can be created to check if the Simpson values estimated for the metagenomes are normally distributed:
hist(richness4$Observed, main="Observed index", xlab="")

# If it is normally distributed, we can run anova tests and check if the variable Stage impacts the Observed index:
anova.sh = aov(richness4$Observed ~ sample_data(physeq1)$Stage)
summary(anova.sh)
TukeyHSD(anova.sh)

# For post hoc test, generate the TukeyHSD in a graph:
tukey.plot.test <- TukeyHSD(anova.sh)
plot(tukey.plot.test, las = 1)

# For Non-normally distributed data, we can use Kruskal-Wallis Rank Sum Test:
kruskal.test(richness4$Observed ~ sample_data(physeq1)$Stage)
# We can also get a list with the p-values resulted of the Wilcoxon Tests considering each pair of groups:
pairwise.wilcox.test(richness4$Observed, sample_data(physeq1)$Stage, p.adj = "bonf")

# Modified and annotated by Guia Angela J. Castillo