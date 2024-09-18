### QIIME2R ANALYSIS SCRIPT ###

# Setting the working directory
path <- setwd("C:/Users/Cobb2/Desktop/18.09.20")


# Listing files in the directory
list.files(path)


# Loading the required packages
library(qiime2R)
library(ggplot2)
library(tidyverse)


# Reading feature table file into R
SVs <- read_qza("filtered-table-AGE-Healthy-only.qza") # ASV table

names(SVs) # Getting extra .qza-related data on the ASV table

SVs$data[1:5, 1:5] # Accessing the actual data stored in the .qza file

SVs$type # Getting artifact type for file

SVs$uuid # Getting unique identifier for file

SVs$contents # Getting complete list of files within the .qza file and their sizes


# Reading in the metadata file
metadata <- read_q2metadata("sample-metadata.tsv")
View(metadata)

head(metadata, 10) # Getting top ten lines from the metadata



# Reading in the taxonomy data
taxonomy <- read_qza("taxonomy.qza")

head(taxonomy$data, 10) # Viewing the top 10 lines in the taxonomy data

taxonomy <- parse_taxonomy(taxonomy$data) # Breaking up taxonomy strings

head(taxonomy, 10)


# Creating Phyloseq Object
physeq <- qza_to_phyloseq(
  features = "filtered-table-AGE-Healthy-only.qza",
  tree = "tree.qza",
  "taxonomy.qza",
  metadata = "sample-metadata.tsv"
)

physeq # Viewing the phyloseq objects

##############################################################################################
# Reading in the shannon vector for analysis
shannon <- read_qza("shannon_vector.qza")

shannon <- shannon$data %>% 
  rownames_to_column("SampleID") # Moving sample names to a new column to allow for merging later


gplots::venn(list(metadata = metadata$SampleID, shannon = shannon$SampleID)) # Finding samples with shannon diversity value(s) 
# 10 samples don't have shannon values because they were filtered out (pos and neg controls) 


# Adding shannon values as a column to the metadata
metadata <- 
  metadata %>% 
  left_join(shannon)

head(metadata)  

library(dplyr)
group_by(metadata, Status) %>% # getting summary statistics for shannon diversity
  filter(!is.na(shannon_entropy)) %>%
  summarise(
    count= n(), 
    mean = mean(shannon_entropy), 
    sd = sd(shannon_entropy), 
    median = median(shannon_entropy), 
    IQR = IQR(shannon_entropy),
    min = min(shannon_entropy),
    max = max(shannon_entropy)
  )



# Making a boxplot of shannon data and other variable(s)
D <- metadata %>% 
  filter(!is.na(shannon_entropy)) %>% 
  ggplot(aes(x = Status, y = shannon_entropy, color = Status)) +
  geom_boxplot() + 
  geom_point(position = "jitter", size = 5) + 
  xlab("Status") + 
  ylab("Shannon Diversity") + 
  theme_bw() + 
  scale_fill_manual("indianred", "cornflowerblue") + 
  theme(legend.position = "none") + 
  theme(text = element_text(size = 14, family = "TT Times New Roman")) +
  ggtitle("D")


# Reading in the evenness vector for analysis
evenness <- read_qza("evenness_vector.qza")

evenness <- evenness$data %>% 
  rownames_to_column("SampleID") # Moving sample names to a new column to allow for merging later


gplots::venn(list(metadata = metadata$SampleID, evenness = evenness$SampleID)) # Finding samples with shannon diversity value(s) 
# 10 samples don't have evenness values because they were filtered out 


# Adding evenness values as a column to the metadata
metadata <- 
  metadata %>% 
  left_join(evenness)

head(metadata)  


# Making a boxplot of evenness data and other variable(s)
B <- metadata %>% 
  filter(!is.na(pielou_evenness)) %>% 
  ggplot(aes(x = Status, y = pielou_evenness, color = Status)) +
  geom_boxplot() + 
  geom_point(position = "jitter", size = 5) + 
  xlab("Status") + 
  ylab("Pielou Evenness") + 
  theme_bw() + 
  scale_fill_manual("indianred", "cornflowerblue") + 
  theme(legend.position = "none")+
  theme(text = element_text(size = 14, family = "TT Times New Roman")) +
  ggtitle("B")



# Reading in the faith pd vector for analysis
faithpd <- read_qza("faith_pd_vector.qza")

faithpd <- faithpd$data %>% 
  rownames_to_column("SampleID") # Moving sample names to a new column to allow for merging later


gplots::venn(list(metadata = metadata$SampleID, faithpd = faithpd$SampleID)) # Finding samples with shannon diversity value(s) 
# 10 samples don't have faith pd values because they were filtered out 


# Adding faith pd values as a column to the metadata
metadata <- 
  metadata %>% 
  left_join(faithpd)

head(metadata)  


# Making a boxplot of faith pd data and other variable(s)
C <- metadata %>% 
  filter(!is.na(faith_pd)) %>% 
  ggplot(aes(x = Status, y = faith_pd, color = Status)) +
  geom_boxplot() + 
  geom_point(position = "jitter", size = 5) + 
  xlab("Status") + 
  ylab("Faith's Phylogenetic Diversity") + 
  theme_bw() + 
  scale_fill_manual("indianred", "cornflowerblue") + 
  theme(legend.position = "none") +
  theme(text = element_text(size = 14, family = "TT Times New Roman")) +
  ggtitle("C")
  


# Reading in the observed features vector for analysis
observed_feat <- read_qza("observed_features_vector.qza")

observed_feat <- observed_feat$data %>% 
  rownames_to_column("SampleID") # Moving sample names to a new column to allow for merging later


gplots::venn(list(metadata = metadata$SampleID, observed_feat = observed_feat$SampleID)) # Finding samples with shannon diversity value(s) 
# 10 samples don't have observed features values because they were filtered out 


# Adding observed features values as a column to the metadata
metadata <- 
  metadata %>% 
  left_join(observed_feat)

head(metadata)  


# Making a boxplot of observed features data and other variable(s)
A <- metadata %>% 
  filter(!is.na(observed_features)) %>% 
  ggplot(aes(x = Status, y = observed_features, color = Status)) +
  geom_boxplot() + 
  geom_point(position = "jitter", size = 5) + 
  xlab("Status") + 
  ylab("Observed Features") + 
  theme_bw() + 
  scale_fill_manual("indianred", "cornflowerblue") + 
  theme(legend.position = "none") +
  theme(text = element_text(size = 14, family = "TT Times New Roman")) +
  ggtitle("A")


# Making a grid of alpha diversity metrics
library(gridExtra)

grid.arrange(A, B, C, D) # Making a 2*2 grid for the four alpha diversity metrics
ggsave("Alpha diversity ", height = 10, width = 15, device = "pdf", useDingbats = F, dpi = 600)

###############################################################################################
### TRYING TO DO THE ALPHA CORRELATION WITH REFORMATTED METADATA FILE ###


# Reading in the metadata file and excluding empty cells and NCs and PCs from Status column
metadata_2 <- read.csv("sample-metadata.csv", header = T, na.strings = c("N/A", "", " ", "NC", "PC1", "PC2", "PC3", "PC4"))
View(metadata_2)




# Reading in the shannon vector for analysis
shannon <- read_qza("shannon_vector.qza")

shannon <- shannon$data %>% 
  rownames_to_column("SampleID") # Moving sample names to a new column to allow for merging later


gplots::venn(list(metadata = metadata_2$SampleID, shannon = shannon$SampleID)) # Finding samples with shannon diversity value(s) 
# 10 samples don't have shannon values because they were filtered out (pos and neg controls) 


# Adding shannon values as a column to the metadata
metadata_2 <- 
  metadata_2 %>% 
  left_join(shannon)

head(metadata_2)  
View(metadata_2)

# Shannon diversity correlation with age
cor.test(metadata_2$Age_In_Months, metadata_2$shannon_entropy, formula = "spearman")

# Visualising correlation in a plot
theme_set(
  theme_bw() # setting theme
)

ggplot(metadata_2, aes(
  x = Age_In_Months, y = shannon_entropy)) + 
  geom_point()

# Linear modelling of shannon diversity
lm1 <- lm(metadata_2$shannon_entropy ~ metadata_2$Status : metadata_2$Age_In_Months) # Interacting

lm2 <- lm(metadata_2$shannon_entropy ~ metadata_2$Status + metadata_2$Age_In_Months)

lm1
lm2




####################################################################################################
# Importing beta diversity metrics pcoa results from qiime2 for plotting in R
# Required packages (qiime2r and tidyverse) loaded from the top of the page

# Reading in required files
metadata <- read_q2metadata("sample-metadata.tsv")

braycurtis <- read_qza("bray_curtis_pcoa_results.qza") 

jaccard <- read_qza("jaccard_pcoa_results.qza")

uwunifrac <- read_qza("unweighted_unifrac_pcoa_results.qza")

wunifrac <- read_qza("weighted_unifrac_pcoa_results.qza")


# Bray-Curtis PCoA plot
braycurtis$data$Vectors %>% 
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata) %>% 
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(alpha = 0.5, size = 5, aes(color = Status)) + # modifying point appearance
  stat_ellipse(aes(x = PC1, y = PC2, color = Status), type = "norm") + # adding ellipse
  theme_bw() +
  theme(text = element_text(size = 14, family = "TT Times New Roman")) +
  ggtitle("A")
ggsave("BC.pdf", height = 3, width = 3, device = "pdf")

# Jaccard PCoA plot
jaccard$data$Vectors %>% 
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>% 
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(alpha = 0.5, size = 5, aes(color = Status)) + # modifying point appearance
  stat_ellipse(aes(x = PC1, y = PC2, color = Status), type = "norm") + # adding ellipse
  theme_bw() + 
  theme(text = element_text(size = 14, family = "TT Times New Roman")) +
  ggtitle("B")

# Unweighted UniFrac PCoA plot
uwunifrac$data$Vectors %>% 
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata) %>% 
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(alpha = 0.5, size = 5, aes(color = Status)) + # modifying point appearance
  stat_ellipse(aes(x = PC1, y = PC2, color = Status), type = "norm") + # adding ellipse
  theme_bw() +
  theme(text = element_text(size = 14, family = "TT Times New Roman")) +
  ggtitle("C")

# Weighted UniFrac PCoA plot
wunifrac$data$Vectors %>% 
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata) %>% 
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(alpha = 0.5, size = 5, aes(color = Status)) + # modifying point appearance
  stat_ellipse(aes(x = PC1, y = PC2, color = Status), type = "norm") + # adding ellipse
  theme_bw() +
  theme(text = element_text(size = 14, family = "TT Times New Roman")) +
  ggtitle("D")

# gridding the plots
library(gridExtra)
grid.arrange(bc, j, uw, wu)


# Changing font style (using wunifrac as e.g.)
wunifrac$data$Vectors %>% 
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata) %>% 
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(alpha = 0.5, size = 5, aes(color = Status)) + # modifying point appearance
  stat_ellipse(aes(x = PC1, y = PC2, color = Status), type = "norm") + # adding ellipse
  theme_bw() +
  theme(text = element_text(size = 16, family = "TT Times New Roman")) + # Times New Roman Font
  ggtitle("Trial")

# using the gridExtra package to get multipanel plots
library(gridExtra)
grid.arrange(bc, j, uw, wu) # Ignore


#################################################################################################
# Plotting heatmap


metadata <- read_q2metadata("sample-metadata.tsv")

SVs <- read_qza("filtered-table-AGE-Healthy-only.qza")$data

taxonomy <- read_qza("taxonomy.qza")$data

SVs <- apply(SVs, 2, function(x) x/sum(x)*100) # converting to percentage relative abundance


SVsToPlot <- 
  data.frame(MeanAbundance = rowMeans(SVs)) %>% # Finding the average abundance of ASVs
  rownames_to_column("Feature.ID") %>% 
  arrange(desc(MeanAbundance)) %>% 
  top_n(30, MeanAbundance) %>% # change to any number of choice
  pull(Feature.ID)


SVs %>% 
  as.data.frame() %>% 
  rownames_to_column("Feature.ID") %>% 
  gather(-Feature.ID, key = "SampleID", value = "Abundance") %>% 
  mutate(Feature.ID = if_else(Feature.ID %in% SVsToPlot, Feature.ID, "Remainder")) %>% # Features to be collapsed are flagged here
  group_by(SampleID, Feature.ID) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  left_join(metadata) %>% 
  mutate(NormAbundance = log10(Abundance + 0.01)) %>% # Doing a log10 transformation after adding a 0.01% pseudocount
  left_join(taxonomy) %>% 
  mutate(Feature = paste(Feature.ID, Taxon)) %>% 
  mutate(Feature = gsub("[kpcofgs]_", "", Feature)) %>% # Trimming out leading text (prefix) from taxonomy
  ggplot(aes(x = SampleID, y = Feature, fill = NormAbundance)) + 
  geom_tile() + 
  facet_grid(~'Status', scales = "free_x") + 
  theme_q2r() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 11, family = "TT Times New Roman")) + # Times New Roman font style
  scale_fill_viridis_c(name = "log10(% Abundance)")
ggsave("Heatmap-top-30-taxa.pdf", height = 5, width = 12, device = "pdf")  




##############################################################################################
## Differential abundance analysis

# Reset the directory for new files to appear on the path list of files
# Read in packages and required files

library(ggrepel)
library(ggtree)
library(ape)
library(ggplot2)
library(tidyverse)
library(qiime2R)

metadata <- read_q2metadata("sample-metadata.tsv")

SVs <- read_qza("filtered-table-AGE-Healthy-only.qza")$data

taxonomy <- read_qza("taxonomy.qza")$data

results <- read_qza("differentials.qza")$data

tree <- read_qza("tree.qza")$data

View(results)

# Generating a volcano plot

results %>%
  left_join(taxonomy) %>% 
  mutate(Significant = if_else(we.eBH < 0.1, TRUE, FALSE)) %>% 
  mutate(Taxon = as.character(Taxon)) %>% 
  mutate(TaxonToPrint = if_else(we.eBH < 0.1, Taxon, "")) %>% # Labeling only significant results
  ggplot(aes(x = diff.btw, y = -log10(we.ep), color = Significant, label = TaxonToPrint)) + 
  geom_text_repel(size = 1, nudge_y = 0.05) + # text labels on volcano plot
  geom_point(alpha = 0.6, shape = 16) + 
  theme_q2r() + 
  xlab("log2 (fold change)") + 
  ylab("-log10 (P-value)") + 
  theme(legend.position = "none") + 
  scale_color_manual(values = c("black", "red"))
ggsave("Volcano-plot.pdf", height = 3, width = 3, device = "pdf")


# Changing parameters from preceding volcano plot
results %>%
  left_join(taxonomy) %>% 
  mutate(Significant = if_else(we.eBH < 0.1, TRUE, FALSE)) %>% # Changed the BH corrected p-values
  mutate(Taxon = as.character(Taxon)) %>% 
  mutate(TaxonToPrint = if_else(we.eBH < 0.1, Taxon, "")) %>% # Labeling only significant results
  ggplot(aes(x = diff.btw, y = -log10(we.ep), color = Significant)) + 
  geom_point(alpha = 0.6, shape = 16, size = 3) + 
  theme_bw() + 
  xlab("log2 (fold change)") + 
  ylab("-log10 (P-value)") + 
  theme(legend.position = "none") + 
  scale_color_manual(values = c("black", "red")) +
  theme(text = element_text(size = 12, family = "TT Times New Roman"))
ggsave("Volcano-plot-no-labels.pdf", height = 3, width = 3, device = "pdf")


## Aldex2 results from qiime2 feature extraction ##
library(qiime2R)
sig_status <-read_qza("sig_status.qza")

View(sig_status$data) ## shows the significant differentially abundant feature list

# Importing and viewing same differentially abundant features exported from QIIME2
aldex_sig <- read.delim("aldex-sig-differentials.tsv", header = T)

View(aldex_sig)

write.table(aldex_sig, "aldex-sig-resaved.txt")

library(tidyverse)
library(writexl)

sel_aldex <- aldex_sig %>% # selecting only columns of interest
  select(featureid, Taxonomy, rab.win.AGE, rab.win.Healthy, diff.btw, we.eBH)

View(sel_aldex)

write_xlsx(sel_aldex, "selected-aldex-results.xlsx") # saving as .xlsx format

###############################################################################################



## Per-feature plotting ##

# first feature trial
clr1 <- apply(log2(SVs + 0.5), 2, function(x) x-mean(x))

clr1 %>% 
  as.data.frame() %>% 
  rownames_to_column("Feature.ID") %>% 
  gather(-Feature.ID, key = SampleID, value = CLR) %>% 
  filter(Feature.ID == "35ffcc3b809d667286737d79670b8de5") %>% 
  left_join(metadata) %>% 
  ggplot(aes(x = Status, y = CLR, fill = Status)) + 
  stat_summary(geom = "bar", color = "black") + 
  geom_jitter(width = 0.2, height = 0, shape = 21, size = 5) + 
  theme_bw() + 
  theme(legend.position = "none")

# trying the violin plot
clr1 %>% 
  as.data.frame() %>% 
  rownames_to_column("Feature.ID") %>% 
  gather(-Feature.ID, key = SampleID, value = CLR) %>% 
  filter(Feature.ID == "35ffcc3b809d667286737d79670b8de5") %>% 
  left_join(metadata) %>% 
  ggplot(aes(x = Status, y = CLR, fill = Status)) + 
  geom_violin() + 
  geom_boxplot(width = 0.1, color = "grey", alpha = 0.2) +
  theme_bw() + 
  theme(legend.position = "none")


## Plotting a phylogenetic tree

results <- 
  results %>% 
  mutate(Significant = if_else(we.eBH < 0.1, "*", ""))


tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% results$Feature.ID]) # Removing all features from the tree we don't data for


ggtree(tree, layout = "circular") %<+% results + 
  geom_tippoint(aes(fill = diff.btw), shape = 21, color = "grey50") + 
  geom_tiplab2(aes(label = Significant), size = 10) + 
  scale_fill_gradient2(low = "darkblue", high = "darkred", midpoint = 0, mid = "white", name = "log2 (fold-change)") + 
  theme(legend.position = "right")
ggsave("tree_2.pdf", height = 15, width = 10, device = "pdf", useDingbats = F, dpi = 600)


##############################################################################################
# saving qiime2 .tsv metadata as .csv
meta <- read.delim("sample-metadata.tsv")
View(meta)

write.csv(meta, 'sample-metadata.csv')

##############################################################################################






























































