## Decontam package to remove contaminants ##

# Setting working directory
path <- setwd("~/Desktop/Seqs/bacteria/R_analyses")

# Listing files in directory
list.files(path)


# Loading required packages
library(phyloseq)
library(qiime2R)
library(decontam)
library(ggplot2)
library(microbiome)

library(showtext) # to import font style from pc into R

# importing Helvetica
# font_add(family = "Arial", regular = "Arial.ttf")

font_add(family = "Helvetica", regular = "Helvetica.ttc") # use this

showtext_auto()



# Creating phyloseq object (original code used)

insertion-tree <- qza_to_phyloseq("insertion-tree-bac.qza") # code didn't work 
# Error in FUN(X[[i]], ...) : 
# numbers of left and right parentheses in Newick string not equal)..reported to qiime2R GitHub  

ps <- qza_to_phyloseq(
  features = "filtered-table-bac.qza", # phyloseq object without tree data
  taxonomy = "taxonomy-bac.qza", # taxonomy-bac-silva.qza for SILVA 138 , taxonomy-bac.qza for Greengenes
  metadata = "metadata.tsv"
)


# ps_silva <- qza_to_phyloseq(
  features = "filtered-table-bac.qza", # phyloseq object without tree data
  taxonomy = "taxonomy-bac-silva.qza", # taxonomy-bac-silva.qza for SILVA 138 , taxonomy-bac.qza for Greengenes
  metadata = "metadata.tsv"
)

ps # 3037 taxa/ASVs

View(data.frame(otu_table(ps)))

plot_bar(ps, fill = "Phylum")

# Viewing the metadata details 
head(sample_data(ps))
tail(sample_data(ps))
View(sample_data(ps))


## DECONTAM step to remove contaminants ##

# Inspecting library sizes and plotting
df <- as.data.frame(sample_data(ps)) # data frame

df$LibrarySize <- sample_sums(ps)

df <- df[order(df$LibrarySize),]

df$Index <- seq(nrow(df)) # [library size/sequence depth against samples]

ggplot(data = df,
       aes(x = Index, y = LibrarySize, color = Status)) +
  geom_point() +
  theme_bw()

# Identifying contaminants using Prevalence
sample_data(ps)$is.neg <- sample_data(ps)$Status == "NC"

contamdf.prev <- isContaminant(ps, 
                               method = "prevalence", 
                               neg = "is.neg")

table(contamdf.prev$contaminant)
# number of contaminants was 17 

head(which(contamdf.prev$contaminant))

# Using a stringent threshold of 0.5
contamdf.prev05 <- isContaminant(ps,
                                 method = "prevalence",
                                 neg = "is.neg",
                                 threshold = 0.5)

table(contamdf.prev05$contaminant)
# number of contaminants was 51

head(which(contamdf.prev05$contaminant))
which(contamdf.prev05$contaminant)

## Removing contaminants from 0.5 prevalence object ##
# Have a look at phyloseq obj before removal
ps

ps.noncontam <- prune_taxa(!contamdf.prev05$contaminant, ps)

# "noncontaminant" phyloseq object
ps.noncontam # 2986


plot_bar(ps.noncontam, fill = "Phylum")


# Removing negative controls (NCs), named as B1-B6 ("blanks")
ps.noncontam_noNC <- prune_samples(
  sample_names(ps.noncontam) !="B1", ps.noncontam)

ps.noncontam_noNC <- prune_samples(
  sample_names(ps.noncontam_noNC) !="B2", ps.noncontam_noNC)

ps.noncontam_noNC <- prune_samples(
  sample_names(ps.noncontam_noNC) !="B3", ps.noncontam_noNC)

ps.noncontam_noNC <- prune_samples(
  sample_names(ps.noncontam_noNC) !="B4", ps.noncontam_noNC)

ps.noncontam_noNC <- prune_samples(
  sample_names(ps.noncontam_noNC) !="B5", ps.noncontam_noNC)

ps.noncontam_noNC <- prune_samples(
  sample_names(ps.noncontam_noNC) !="B6", ps.noncontam_noNC)


plot_bar(ps.noncontam_noNC, fill = "Phylum")


# Removing positive controls (PCs), named as PC1-PC4
ps.noncontam_noPC <- prune_samples(
  sample_names(ps.noncontam_noNC) !="PC1", ps.noncontam_noNC)

ps.noncontam_noPC <- prune_samples(
  sample_names(ps.noncontam_noPC) !="PC2", ps.noncontam_noPC)

ps.noncontam_noPC <- prune_samples(
  sample_names(ps.noncontam_noPC) !="PC3", ps.noncontam_noPC)

ps.noncontam_noPC <- prune_samples(
  sample_names(ps.noncontam_noPC) !="PC4", ps.noncontam_noPC)


plot_bar(ps.noncontam_noPC, fill = "Phylum")

# Renaming phyloseq object
ps.noncontam_noctrls <- ps.noncontam_noPC

ps.noncontam_noctrls


# Confirming with taxa bar plot
plot_bar(ps.noncontam_noctrls, fill = "Phylum") +
  facet_wrap(~Status, scales = "free_x")





#### Different analysis set ####

# https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html
# https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/


# Trying something out
rank_names(ps.noncontam_noctrls)

table(tax_table(ps.noncontam_noctrls)[, "Phylum"],
      exclude = NULL)


# Removing Phylum with NA
ps.noncontam_noctrls <- subset_taxa(
  ps.noncontam_noctrls,
  !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

table(tax_table(ps.noncontam_noctrls)[, "Phylum"],
      exclude = NULL)

ps.noncontam_noctrls # 2897 taxa by 107 samples



######################### UNNECESSARY (EXPLORATORY) ############################
# Agglomerating to Phylum and genus levels
ps.phylum <- phyloseq::tax_glom(ps.noncontam_noctrls,
                                                  "Phylum")

phyloseq::taxa_names(ps.phylum) <- phyloseq::tax_table(ps.phylum)[, "Phylum"]

phyloseq::otu_table(ps.phylum)[1:5, 1:5] # viewing first five rows and columns



library(tidyverse)

# Melt and plot
phyloseq::psmelt(ps.phylum) %>%
  ggplot(data = .,
         aes(x = Status, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  facet_wrap(~ Phylum, scales = "free", nrow = 3) + # 3 rows
  scale_color_viridis_d() +
  theme_bw() +
  theme(text = element_text(size = 19)) + # font size of 19
  theme(legend.position = "none") # removing legend 
ggsave("Phylum_boxplot.pdf", width = 18, height = 14)


# Plot with relative abundance data 
phyloseq::psmelt(ps.rel) %>%
  ggplot(data = .,
         aes(x = Status, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  facet_wrap(~ Phylum, scales = "free", nrow = 3) +
  scale_color_viridis_d() +
  theme_bw() +
  theme(text = element_text(size = 18)) + # font size of 19
  theme(legend.position = "none")
################################################################################




################################################################################
######## Testing for overall phylum group differences using HMP package ########

# subsetting groups
HCs <- phyloseq::subset_samples(ps.phylum,
                                Status == "Healthy")

AGE <- phyloseq::subset_samples(ps.phylum,
                                Status == "AGE")


# subsetting ASV/OTU table
HC_otu <- data.frame(phyloseq::otu_table(HCs))

AGE_otu <- data.frame(phyloseq::otu_table(AGE))

# grouping rare phyla *those with ASV counts < 10*
HC_otu <- HC_otu %>%
  t(.) %>%
  as.data.frame(.) %>%
  mutate(Other = 
           Elusimicrobia +
           Euryarchaeota +
           Lentisphaerae +
           OP3 +
           Planctomycetes +
           Spirochaetes +
           SR1 +
           Synergistetes +
           TM6 +
           TM7 +
           Verrucomicrobia) %>%
  dplyr::select(-Elusimicrobia,
      -Euryarchaeota,
      -Lentisphaerae,
      -OP3,
      -Planctomycetes,
      -Spirochaetes,
      -SR1,
      -Synergistetes,
      -TM6,
      -TM7,
      -Verrucomicrobia)



AGE_otu <- AGE_otu %>%
  t(.) %>%
  as.data.frame(.) %>%
  mutate(Other = 
           Elusimicrobia +
           Euryarchaeota +
           Lentisphaerae +
           OP3 +
           Planctomycetes +
           Spirochaetes +
           SR1 +
           Synergistetes +
           TM6 +
           TM7 +
           Verrucomicrobia) %>%
  dplyr::select(-Elusimicrobia,
                -Euryarchaeota,
                -Lentisphaerae,
                -OP3,
                -Planctomycetes,
                -Spirochaetes,
                -SR1,
                -Synergistetes,
                -TM6,
                -TM7,
                -Verrucomicrobia)



# Testing for group differences
library(HMP)

group_data <- list(HC_otu, AGE_otu)

xdc <- HMP::Xdc.sevsample(group_data)

xdc # there are differences in phyla 

pchisq(85.11739, 8)
################################################################################



########################### Top 10 Taxa-level barplots #########################

library(ggplot2)
library(fantaxtic)


# Phylum
ps.phylum <- phyloseq::tax_glom(ps.noncontam_noctrls, # agglomerating to phylum level
                                "Phylum")


top10.phy <- get_top_taxa(physeq_obj = ps.phylum, 
                          n = 10, 
                          relative = T, 
                          discard_other = F, 
                          other_label = "Other")


phylum_bar <- fantaxtic_bar(top10.phy, 
              color_by = "Phylum", 
              other_label = "Other",
              order_alg = "alph", 
              facet_by = "Status",
              facet_cols = 2) + 
  ylab("Relative Abundance") +
  theme(text = element_text(size = 12, family = "Helvetica"),
        axis.text.x = element_blank(), # removing x-axis labels
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 12, family = "Helvetica"),# increasing facet label size
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 10, family = "Helvetica"),
        legend.title = element_text(size = 12, family = "Helvetica")) 

phylum_bar

# Class
ps.class <- phyloseq::tax_glom(ps.noncontam_noctrls,
                                "Class")


top10.class <- get_top_taxa(physeq_obj = ps.class, 
                          n = 10, 
                          relative = T, 
                          discard_other = F, 
                          other_label = "Other")


class_bar <- fantaxtic_bar(top10.class, 
              color_by = "Class", 
              other_label = "Other", 
              order_alg = "alph", 
              facet_by = "Status",
              facet_cols = 2) + 
  ylab("Relative Abundance") +
  theme(text = element_text(size = 14, family = "Helvetica"),
        axis.text.x = element_blank(), # removing x-axis labels
        axis.title.x = element_blank(),
        strip.text.x = element_blank(), # removing facet label
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 12, family = "Helvetica"),
        legend.title = element_text(size = 14, family = "Helvetica")) 



# Order
ps.order <- phyloseq::tax_glom(ps.noncontam_noctrls,
                               "Order")


top10.order <- get_top_taxa(physeq_obj = ps.order, 
                            n = 10, 
                            relative = T, 
                            discard_other = F, 
                            other_label = "Other")


fantaxtic_bar(top10.order, 
              color_by = "Order", 
              other_label = "Other", 
              order_alg = "alph", 
              facet_by = "Status") + 
  ylab("Relative Abundance") +
  theme(text = element_text(size = 14, family = "Helvetica"),
        axis.text.x = element_blank(), # removing x-axis labels
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 14, family = "Helvetica"),# increasing facet label size
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 12, family = "Helvetica"),
        legend.title = element_text(size = 14, family = "Helvetica")) 



# Family
ps.family <- phyloseq::tax_glom(ps.noncontam_noctrls,
                               "Family")


top10.family <- get_top_taxa(physeq_obj = ps.family, 
                            n = 10, 
                            relative = T, 
                            discard_other = F, 
                            other_label = "Other")


family_bar <- fantaxtic_bar(top10.family, 
              color_by = "Family", 
              other_label = "Other", 
              order_alg = "alph", 
              facet_by = "Status",
              facet_cols = 2) + 
  ylab("Relative Abundance") +
  theme(text = element_text(size = 12, family = "Helvetica"),
        axis.text.x = element_blank(), # removing x-axis labels
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),# removing facet labels cos it's already in the phylum stacked barplot
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 10, family = "Helvetica"),
        legend.title = element_text(size = 12, family = "Helvetica")) 


# Genus
ps.genera <- phyloseq::tax_glom(ps.noncontam_noctrls,
                                                  "Genus")

ps.genera

top10.genera <- get_top_taxa(physeq_obj = ps.genera, 
                          n = 10, 
                          relative = T, 
                          discard_other = F, 
                          other_label = "Other")


genus_bar <- fantaxtic_bar(top10.genera, 
              color_by = "Genus", 
              other_label = "Other", 
              order_alg = "alph", 
              facet_by = "Status",
              facet_cols = 2) + 
  ylab("Relative Abundance") +
  theme(text = element_text(size = 12, family = "Helvetica"),
        axis.text.x = element_blank(), # removing x-axis labels
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 10, family = "Helvetica"),
        legend.title = element_text(size = 12, family = "Helvetica")) 


library(ggpubr) # grid


ggarrange(phylum_bar, family_bar, genus_bar, # phylum, family, and genus
          ncol = 1,
          nrow = 3,
          align = "hv",
          labels = c("(A)", "(B)", "(C)")) # exported manually as .pdf file using the export feature



########################## Testing for taxa differences ########################
####### Wilcoxon rank sum test (not recommended for compositional data) ########

# Agglomerating to taxa level of choice

# Different names for agglomerated data than earlier---cos I did this on Windows pc before transfering the code onto the Mac

ps.noncontam.noctrls.phylum <- phyloseq::tax_glom(ps.noncontam_noctrls,
                                                  "Phylum")

ps.noncontam.noctrls.phylum


ps.noncontam.noctrls.class <- phyloseq::tax_glom(ps.noncontam_noctrls,
                                                 "Class")

ps.noncontam.noctrls.class


ps.noncontam.noctrls.order <- phyloseq::tax_glom(ps.noncontam_noctrls,
                                                 "Order")

ps.noncontam.noctrls.order


ps.noncontam.noctrls.family <- phyloseq::tax_glom(ps.noncontam_noctrls,
                                                  "Family")

ps.noncontam.noctrls.family


ps.noncontam.noctrls.genus <- phyloseq::tax_glom(ps.noncontam_noctrls,
                                                 "Genus")


ps.noncontam.noctrls.genus



ps.noncontam.noctrls.species <- phyloseq::tax_glom(ps.noncontam_noctrls,
                                                   "Species")


ps.noncontam.noctrls.species


# loading packages for data manipulation
library(magrittr)
library(broom)
library(Hmisc)
library(dplyr)


# transform_sample_counts(ps, function(x) x/sum(x)*100)
# PHYLUM
ps.noncontam.noctrls.phylum.rel <- transform_sample_counts(
  ps.noncontam.noctrls.phylum, function(x) x/sum(x)*100)


## Converting to data frame
ps.phy.rel <- psmelt(ps.noncontam.noctrls.phylum.rel)

View(ps.phy.rel)

## testing relative abundance differences
ps.phy.rel %>%
  group_by_("Phylum") %>%
  do(tidy(wilcox.test(Abundance ~ Status, 
                      data = .,
                      exact = F))) %>% 
  ungroup() %>% 
  mutate(p.adjust = p.adjust(p.value, method = "BH")) -> mw.phylum.results


View(mw.phylum.results)
# test summary

mw.phylum.results %>%
  subset(p.adjust < .05) %>%
  knitr::kable()



View(group_by(ps.phy.rel, Phylum) %>% # getting summary statistics for overall phyla abundances
       summarise(
         count= n(), 
         mean = mean(Abundance), 
         sd = sd(Abundance), 
         median = median(Abundance), 
         IQR = IQR(Abundance),
         min = min(Abundance),
         max = max(Abundance)
       ))



View(group_by(ps.phy.rel, Phylum, Status) %>% # getting summary statistics for phyla grouped by status
       summarise(
         count= n(), 
         mean = mean(Abundance), 
         sd = sd(Abundance), 
         median = median(Abundance), 
         IQR = IQR(Abundance),
         min = min(Abundance),
         max = max(Abundance)
       ))


#### CLASS

ps.noncontam.noctrls.class.rel <- transform_sample_counts(
  ps.noncontam.noctrls.class, function(x) x/sum(x)*100)


## Converting to data frame
ps.class.rel <- psmelt(ps.noncontam.noctrls.class.rel)

View(ps.class.rel)

## testing relative abundance differences
ps.class.rel %>%
  group_by_("Class") %>%
  do(tidy(wilcox.test(Abundance ~ Status, 
                      data = .,
                      exact = F))) %>% 
  ungroup() %>% 
  mutate(p.adjust = p.adjust(p.value, method = "BH")) -> mw.class.results

# test summary

mw.class.results %>%
  subset(p.adjust < .05) %>%
  knitr::kable()



View(group_by(ps.class.rel, Class) %>% # getting summary statistics for class abundances
       summarise(
         count= n(), 
         mean = mean(Abundance), 
         sd = sd(Abundance), 
         median = median(Abundance), 
         IQR = IQR(Abundance),
         min = min(Abundance),
         max = max(Abundance)
       ))


View(group_by(ps.class.rel, Class, Status) %>% # getting summary statistics for phyla grouped with status
       summarise(
         count= n(), 
         mean = mean(Abundance), 
         sd = sd(Abundance), 
         median = median(Abundance), 
         IQR = IQR(Abundance),
         min = min(Abundance),
         max = max(Abundance)
       ))


# FAMILY
ps.noncontam.noctrls.family.rel <- transform_sample_counts(
  ps.noncontam.noctrls.family, function(x) x/sum(x)*100) 



## Converting to data frame
ps.family.rel <- psmelt(ps.noncontam.noctrls.family.rel)

View(ps.family.rel)

## testing relative abundance differences
ps.family.rel %>%
  group_by_("Family") %>%
  do(tidy(wilcox.test(Abundance ~ Status, 
                      data = .,
                      exact = F))) %>% 
  ungroup() %>% 
  mutate(p.adjust = p.adjust(p.value, method = "BH")) -> mw.family.results

View(mw.family.results)

# test summary

mw.family.results %>%
  subset(p.adjust < .05) %>%
  knitr::kable()



View(group_by(ps.family.rel, Family) %>% # getting summary statistics for family abundances
       summarise(
         count= n(), 
         mean = mean(Abundance), 
         sd = sd(Abundance), 
         median = median(Abundance), 
         IQR = IQR(Abundance),
         min = min(Abundance),
         max = max(Abundance)
       ))


View(group_by(ps.family.rel, Family, Status) %>% # getting summary statistics for family grouped with status
       summarise(
         count= n(), 
         mean = mean(Abundance), 
         sd = sd(Abundance), 
         median = median(Abundance), 
         IQR = IQR(Abundance),
         min = min(Abundance),
         max = max(Abundance)
       ))




# GENUS
ps.noncontam.noctrls.genus.rel <- transform_sample_counts(
  ps.noncontam.noctrls.genus, function(x) x/sum(x)*100)

ps.noncontam.noctrls.genus.rel

## Converting to data frame
ps.genus.rel <- psmelt(ps.noncontam.noctrls.genus.rel)

View(ps.genus.rel)

## testing relative abundance differences
ps.genus.rel %>%
  group_by_("Genus") %>%
  do(tidy(wilcox.test(Abundance ~ Status, 
                      data = .,
                      exact = F))) %>% 
  ungroup() %>% 
  mutate(p.adjust = p.adjust(p.value, method = "BH")) -> mw.genus.results

# test summary
View(mw.genus.results)

mw.genus.results %>%
  subset(p.adjust < .05) %>%
  knitr::kable()



View(group_by(ps.genus.rel, Genus) %>% # getting summary statistics for class abundances
       summarise(
         count= n(), 
         mean = mean(Abundance), 
         sd = sd(Abundance), 
         median = median(Abundance), 
         IQR = IQR(Abundance),
         min = min(Abundance),
         max = max(Abundance)
       ))


View(group_by(ps.genus.rel, Genus, Status) %>% # getting summary statistics for phyla grouped with status
       summarise(
         count= n(), 
         mean = mean(Abundance), 
         sd = sd(Abundance), 
         median = median(Abundance), 
         IQR = IQR(Abundance),
         min = min(Abundance),
         max = max(Abundance)
       ))
################################################################################




################################################################################
########################## CORE MICROBIOTA FOR ALL #############################
library(microbiome)

ps.genera # 216 taxa by 107 samples

# keeping only taxa with positive sums
ps.genera.taxsums <- prune_taxa(
  taxa_sums(ps.genera) > 0,
  ps.genera
)

ps.genera.taxsums # 212

# compositional version
ps.genera.taxsums.rel <- microbiome::transform(
  ps.genera.taxsums,
  "compositional"
)

# core taxa members
taxa_names(ps.genera.taxsums.rel)[1:5]

core.mem.genera.50 <- core_members(
  ps.genera.taxsums.rel,
  detection = 0.0001, # 0.01%
  prevalence = 0.5 # 50% 
)

core.mem.genera.50

# core taxa
ps.genera.core.50 <- core(
  ps.genera.taxsums.rel,
  detection = 0.0001, # 0.01%
  prevalence = 0.5 # 50% 
)

core.genera <- taxa(ps.genera.core.50)

class(core.genera) # character

tax.mat.core.genera <- tax_table(ps.genera.core.50)

tax.df.core.genera <- as.data.frame(tax.mat.core.genera)

View(tax.df.core.genera) # 21 core genera


# adding compositionals
library(viridis)

p1 <- plot_core(ps.genera.taxsums.rel,
  plot.type = "heatmap",
  colours = gray(seq(0, 1, length = 5)), # previous gray colour
  prevalences = seq(0.05, 1, 0.05),
  detections = round(10^seq(log10(1e-4),
                      log10(0.2),
                      length = 10), 4),
  min.prevalence = 0.5) +
  scale_fill_viridis(option = "A") + # option A for magma palette in viridis package
  theme_bw()

# reprocessing to get taxa names on plot
library(knitr)

# data for plotting 
df <- p1$data

# getting ASV list
list <- df$Taxa

list # checking OTU ids

# getting taxonomy data
tax <- as.data.frame(tax_table(ps.genera.taxsums.rel))

# adding ASVs to last column
tax$ASV <- rownames(tax)

# selecting only those ASVs for plots
tax2 <- dplyr::filter(tax,
                      rownames(tax) %in% list)

head(tax2)

# merging select columns into one
tax.unit <- tidyr::unite(tax2, 
                         Taxa_level, 
                         c("ASV", "Genus"),
                         sep = ":",
                         remove = T)

tax.unit$Taxa_level <- gsub(pattern = "[a-z]_",
                            replacement = "",
                            tax.unit$Taxa_level)

# adding new info to plot data
df$Taxa <- tax.unit$Taxa_level

# checking taxa info
knitr::kable(head(df))
knitr::kable(df)[1:24] # seeing all 21 core members at 1e-4 cutoff

# replacing data in plot object
p1$data <- df

library(scales)

core_heatmap <- plot(p1 + 
       theme(axis.text.y = element_text(face = "italic")) +
       ylab("Genus") + xlab("Detection Threshold (Relative Abundance)") +
       theme(text = element_text(size = 12))) +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x.bottom = element_text())



# core microbiota venn diagram for per group analyses

# https://microbiome.github.io/tutorials/core_venn.html
library(eulerr)
library(microbiome)
library(microbiomeutilities)


ps.genera

ps.genera.taxsums

ps.genera.tax.rel <- microbiome::transform(ps.genera.taxsums,
                                       "compositional")


# number of samples per group
table(meta(ps.genera.taxsums)$Status, useNA = "always") # 57 AGE, 50 HC, No NA samples

# list of status variables
status_states <- unique(as.character(meta(ps.genera.taxsums)$Status))

status_states # AGE and Healthy

# for loop to identify core taxa in each sample
list_core <- c()

for (n in status_states){
  
  ps.genera.sub <- subset_samples(ps.genera.tax.rel, Status == n)
  
  core.genera_m <- core_members(ps.genera.sub,
                                detection = 0.01/100,
                                prevalence = 0.50)
  
  print(paste0("No. of core taxa in ", n, " : ", length(core.genera_m)))
  
  list_core[[n]] <- core.genera_m 
}

print(list_core)


mycols <- c(AGE = "#E69F00", "#56B4E9") # change to match group colours i.e  use scale_color_manual(values = c("#E69F00", "#56B4E9"))

# plot the venn diagram
plot(venn(list_core), # I prefer this to the venn plot below
     fills = mycols)
# 10 for AGE, 30 for HC, and 8 overlapping

taxa_names(ps.genera.tax.rel)[1:5] # top 5 IDs


# getting taxa names
ps.genera.rel.f <- format_to_besthit(ps.genera.tax.rel) # micribomeutilities package function

taxa_names(ps.genera.rel.f)[1:5] # has g__ prefix


# rerun for loop to incorporate new name labels

list_core <- c()

for (n in status_states){
  
  ps.genera.sub <- subset_samples(ps.genera.rel.f, Status == n)
  
  core.genera_m <- core_members(ps.genera.sub,
                                detection = 0.01/100,
                                prevalence = 0.50)
  
  print(paste0("No. of core taxa in ", n, " : ", length(core.genera_m)))
  
  list_core[[n]] <- core.genera_m
}


print(list_core)

# Rothia and Atopobium are unique core to AGE

# recoding venn diagram
mycols <- c(AGE = "#E69F00", "#56B4E9") # change to match group colours

# plot the venn diagram
core_venn <- plot(venn(list_core), 
     fills = mycols,
     element_text(family = "Helvetica"),
     element_text(size = 12))

core_venn


library(ggpubr)

core_heatmap

core_venn


ggarrange(core_heatmap, core_venn, 
          widths = c(1.7, 0.8), # heatmap and venn diagram
          ncol = 2,
          nrow = 1,
          labels = c("(A)", "(B)"))







################################################################################

##################### DIFFERENTIAL ABUNDANCE TESTING ###########################


# NB: remove less abundant otus, agglomerate to genus, and use an alpha cutoff of .05

ps.genera # 216

filter <- phyloseq::genefilter_sample(ps.genera,
                                      filterfun_sample(function(x) x >= 1), # at least 1 read
                                      A = 0.1*nsamples(ps)) # 10 percent of samples

ps.genera_filtered <- prune_taxa(filter, ps.genera)

ps.genera_filtered # 64 genera



#################### DESEQ2 differential abundance testing #####################

# adding pseudocount of 1
otu_table(ps.genera_filtered) <- otu_table(ps.genera_filtered) + 1

# converting to deseq format
library(DESeq2)

dds.data <- phyloseq_to_deseq2(ps.genera_filtered, ~ Status)

dds <- DESeq(dds.data)

res <- results(dds)

res <- res[order(res$padj, na.last = NA), ]

sigtab <- res[res$padj < 0.05, ] # alpha < 0.05

sigtab <- cbind(as(sigtab,
                   "data.frame"),
                as(tax_table(ps.genera_filtered)[rownames(sigtab), ],
                   "matrix"))

sigtab


# bayesian estimation of dispersion
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

plotDispEsts(dds)

sigtabgen <- subset(sigtab, !is.na(Family)) # should change to Phylum but I removed all those without labels already...so this doesn't hurt

# phylum order
x <- tapply(sigtabgen$log2FoldChange,
            sigtabgen$Phylum,
            function(x) max(x))

x <- sort(x, T)

sigtabgen$Phylum <- factor(as.character(sigtabgen$Phylum),
                           levels = names(x))

# family order
x <- tapply(sigtabgen$log2FoldChange, # Family taxa not used. I prefer to have "genus under phylum" rather 
            sigtabgen$Family,
            function(x) max(x))

x <- sort(x, T)

sigtabgen$Family <- factor(as.character(sigtabgen$Family),
                           levels = names(x))

# genus order
x <- tapply(sigtabgen$log2FoldChange, # use this, as it is more informative
            sigtabgen$Genus,
            function(x) max(x))

x <- sort(x, T)

sigtabgen$Genus <- factor(as.character(sigtabgen$Genus),
                           levels = names(x))

View(sigtabgen)
dim(sigtabgen)

fdr_deseq <- sigtabgen %>%
  dplyr::filter(padj < 0.05)

ggplot(sigtabgen,
       aes(y = Genus, x = log2FoldChange, color = Phylum)) +
  geom_vline(xintercept = 0.0, color = "black", size = 0.5, linetype = "dotted") +
  geom_point(size = 5) +
  labs(x = "\nLog2 Fold-Change (AGE Cases vs. Healthy Controls)", y = "Genus-level Features") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) +
  theme(text = element_text(family = "Helvetica")) +
  theme(text = element_text(size = 14))


# another plot for deseq2
write.csv(fdr_deseq, "fdr_dseq.csv") # necessary to manually correct the alloiococcus misclassification  formatted

fdr_deseq_manual_edit <- read.csv("fdr_dseq.csv", header = T)


deseq_diff_abun_plot <- ggplot(fdr_deseq_manual_edit,
                             aes(x = reorder(Genus, -log2FoldChange),
                                 y = log2FoldChange, 
                                 color = ifelse(log2FoldChange < 0, "#E69F00", "#56B4E9"))) + # colour per group
  geom_point(size = 3) +
  coord_flip() +
  ggtitle("DESeq2 Differential Abundance Test") +
  labs(x = "Genus-level Features", y = "Log2 Fold-Change (AGE Cases vs. Healthy Controls)") +
  geom_hline(yintercept = 0, color = "black", size = 0.5, linetype = "dotted") +
  theme_bw() +
  theme(text = element_text(family = "Helvetica", size = 12)) +
  scale_color_identity()

deseq_diff_abun_plot


##### Trying to have a specific order for genera names for easy comparisons ####







##################### ANCOMBC DIFFERENTIAL ABUNDANCE TEST ######################
library(ANCOMBC)

library(microbiomeutilities)

ps.genera_filtered_name <- ps.genera_filtered # renaming physeq obj, so I can use the new obj for other diff abun tests

ps.genera_filtered_name <- format_to_besthit(ps.genera_filtered_name) # to add genera labels to feature IDs



ancom_genera <- ancombc(phyloseq = ps.genera_filtered_name, 
                        formula = "Status",
                        p_adj_method = "BH",
                        zero_cut = 0.90, 
                        lib_cut = 1000, 
                        group = "Status", 
                        struc_zero = T, 
                        neg_lb = T, 
                        tol = 1e-5,
                        max_iter = 1000, # 1000 iterations
                        conserve = T,
                        alpha = 0.05, # alpha cut-off value of 5%, also known as the significant value
                        global = T)

# reformating ancom results
ancom_genera

ancom_genera_df <- data.frame(
  Genus = row.names(ancom_genera$res$beta),
  beta = unlist(ancom_genera$res$beta),
  se = unlist(ancom_genera$res$se),
  W = unlist(ancom_genera$res$W),
  p_val = unlist(ancom_genera$res$p_val),
  q_val = unlist(ancom_genera$res$q_val),
  diff_abn = unlist(ancom_genera$res$diff_abn)
)

# focusing on features with adjusted p vals < 0.05
#ancom_genera_df <- ancom_genera_df[
  (ancom_genera_df$q_val < 0.05), ]

#ancom_genera_df <- ancom_genera_df[
  order(ancom_genera_df$W, na.last = NA),]


fdr_ancom <- ancom_genera_df %>%
  dplyr::filter(q_val < 0.05)

dim(fdr_ancom) # 48 taxa significant

head(fdr_ancom)
tail(fdr_ancom)

View(fdr_ancom)


# plots with feature ID + Genus name (not appealing to the eyes)
#ggplot(fdr_ancom,


# could not figure out how to remove feature IDs using sub/gsub functions
# exporting, manually creating a new column for feature IDs, importing, and replotting

# write.csv(fdr_ancom, "fdr_ancom.csv") # to save ancom results as .csv for manual editting outside R

fdr_ancom_manual_edit <- read.csv("fdr_ancom.csv", header = T)

fdr_ancom_manual_edit



ancom_diff_abun_plot <- ggplot(fdr_ancom_manual_edit,
       aes(x = reorder(Genus, -W),
           y = W, 
           color = ifelse(W < 0, "#E69F00", "#56B4E9"))) + # colour per group
  geom_point(size = 3) +
  coord_flip() +
  ggtitle("ANCOM-BC Differential Abundance Test") +
  labs(x = "Genus-level Features", y = "W (AGE Cases vs. Healthy Controls)") +
  geom_hline(yintercept = 0, color = "black", size = 0.5, linetype = "dotted") +
  theme_bw() +
  theme(text = element_text(family = "Helvetica", size = 12)) +
  scale_color_identity()

ancom_diff_abun_plot


# making a grid plot for the differential abundance for deseq2 and ancom-bc 
library(ggpubr)

ggarrange(deseq_diff_abun_plot, ancom_diff_abun_plot, 
          ncol = 2,
          nrow = 1,
          labels = c("(A)", "(B)"))


################################################################################
############# manually combined the .csv files for deseq2 and ancom-bc results
############# ...so I could have one plot with shared genus axis and "change" for effect sizes
####### can be done using the R console (I couldn't do that here, unfortunately :( ...)

fdr_deseq2_ancom_combined <- read.csv("fdr_deseq2_ancom_combined.csv", 
                                      header = T,
                                      na.strings = c("N/A", "", " "))

View(fdr_deseq2_ancom_combined)

library(dplyr)

fdr_deseq2_ancom_combined %>%
ggplot(aes(x = Genus, 
           y = Change, 
           color = ifelse(Change < 0, "#E69F00", "#56B4E9"))) + # colour per group
  geom_point(size = 3) +
  coord_flip() +
  geom_hline(yintercept = 0, color = "black", size = 0.5, linetype = "dotted") +
  labs(x = "Genus-level Features", y = "AGE Cases vs. Healthy Controls") +
  theme_bw() +
  theme(text = element_text(family = "Helvetica", size = 13),
        strip.text.x = element_text(size = 15)) + # increasing font size for facet header
  scale_color_identity() +
  facet_wrap(~ factor(Tool, levels = c("DESeq2", "ANCOM-BC")), # positioning deseq2 first before ancom-bc
             scales = "free_y") # scales = "free" will allow for diff range of axis values
# note that facet_grid was less intuitive and cluttered



####################################################################################################

##### ANCOM-BC with age as covariate
ancom_genera_age <- ancombc(phyloseq = ps.genera_filtered_name, 
                            formula = "Status + Age_In_Months",
                            p_adj_method = "BH",
                            zero_cut = 0.90, 
                            lib_cut = 1000, 
                            group = "Status", 
                            struc_zero = T, 
                            neg_lb = T, 
                            tol = 1e-5,
                            max_iter = 1000, # 1000 iterations
                            conserve = T,
                            alpha = .05,
                            global = T)

ancom_genera_age

ancom_genera_age_df <- data.frame(
  Genus = row.names(ancom_genera_age$res$beta),
  beta = unlist(ancom_genera_age$res$beta),
  se = unlist(ancom_genera_age$res$se),
  W = unlist(ancom_genera_age$res$W),
  p_val = unlist(ancom_genera_age$res$p_val),
  q_val = unlist(ancom_genera_age$res$q_val),
  diff_abn = unlist(ancom_genera_age$res$diff_abn)
)

ancom_genera_age_df

View(ancom_genera_age_df)

fdr_ancom_age <- ancom_genera_age_df %>%
  dplyr::filter(q_val < 0.05)

dim(fdr_ancom_age)


################################################################################################
######################### PURELY EXPLORATORY DIFF ABUN TESTS ###################################
################## MAASLIN2 DIFFERENTIAL ABUNDANCE TESTING #####################
library(Maaslin2)

ps.genera_filtered_name


maaslin_genera <- Maaslin2(
  input_data = data.frame(otu_table(ps.genera_filtered_name)),
  input_metadata = data.frame(sample_data(ps.genera_filtered_name)),
  output = "~/Desktop/Seqs/bacteria/R_analyses/Maaslin2-LM",
  min_abundance = 0,
  min_prevalence = 0,
  normalization = "TSS",
  transform = "LOG",# NONE can be used too
  analysis_method = "LM",
  max_significance = 0.05,
  fixed_effects = "Status",
  correction = "BH",
  standardize = F,
  cores = 2
)


mas_res_df <- maaslin_genera$results # focusing on the results 

mas_res_df

fdr_mas <- mas_res_df %>% # filtering to include only significant features
  dplyr::filter(qval < 0.05)

dim(fdr_mas) # 47 significant features

# write.csv(fdr_mas, "fdr_mas.csv") # exporting, manually editing, and re-importing for further analysis


fdr_mas_manual_edit <- read.csv("fdr_mas.csv", header = T)

fdr_mas_manual_edit



mas_diff_abun_plot <- ggplot(fdr_mas_manual_edit,
                               aes(x = reorder(Genus, -coef),
                                   y = coef, 
                                   color = ifelse(coef < 0, "#E69F00", "#56B4E9"))) + # colour per group
  geom_point(size = 3) +
  coord_flip() +
  labs(x = "Genus-level Features", y = "Coefficient (AGE Cases vs. Healthy Controls)") +
  geom_hline(yintercept = 0, color = "black", size = 0.5, linetype = "dotted") +
  theme_bw() +
  theme(text = element_text(family = "Helvetica", size = 12)) +
  scale_color_identity()

mas_diff_abun_plot



# adding age of participants as a covariate
maaslin_genera_age <- Maaslin2(
  input_data = data.frame(otu_table(ps.genera_filtered_name)),
  input_metadata = data.frame(sample_data(ps.genera_filtered_name)),
  output = "~/Desktop/Seqs/bacteria/R_analyses/Maaslin2-LM-age-covar",
  min_abundance = 0,
  min_prevalence = 0,
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.05,
  fixed_effects = c("Status", "Age_In_Months"),
  correction = "BH",
  standardize = F,
  cores = 2
)


# making a grid plot for the differential abundance 

ggarrange(ancom_diff_abun_plot, mas_diff_abun_plot, 
          ncol = 2,
          nrow = 1,
          labels = c("(A)", "(B)"))



################################# ALDEX2 #######################################

ps.genera_filtered_name # phyloseq object with genus+feature ID combined

library(ALDEx2)

aldex_genera <- aldex(
  data.frame(phyloseq::otu_table(ps.genera_filtered_name)),
             phyloseq::sample_data(ps.genera_filtered_name)$Status,
             test = "t",
             effect = T,
             denom = "iqlr")
)

fdr_aldex <- aldex_genera %>%
  dplyr::filter(we.eBH < 0.05)

fdr_aldex
dim(fdr_aldex)

# write.csv(fdr_aldex, "fdr_aldex.csv")

fdr_aldex_manual_edit <- read.csv("fdr_aldex.csv", header = T)


aldex_diff_abun_plot <- ggplot(fdr_aldex_manual_edit,
                             aes(x = reorder(Genus, -diff.btw),
                                 y = diff.btw, 
                                 color = ifelse(diff.btw < 0, "#E69F00", "#56B4E9"))) + # colour per group
  geom_point(size = 3) +
  coord_flip() +
  labs(x = "Genus-level Features", y = "Difference (AGE Cases vs. Healthy Controls)") +
  geom_hline(yintercept = 0, color = "black", size = 0.5, linetype = "dotted") +
  theme_bw() +
  theme(text = element_text(family = "Helvetica", size = 12)) +
  scale_color_identity()


aldex_diff_abun_plot
##################################################################################



################################# Selbal #######################################
ps.genera_filtered_name # physeq object with feature ID+genus names

library(selbal)

selbal_genera <- selbal.cv(x = data.frame(t(data.frame(phyloseq::otu_table(ps.genera_filtered_name)))),
                           y = phyloseq::sample_data(ps.genera_filtered_name)$Status,
                           n.fold = 5, 
                           n.iter = 20) # default is 10...default zero.rep used


selbal_variables <- selbal_genera$accuracy.nvar # 1. optimal number of variables = 4 (from using 20 iterations) 

selbal_barplot <- selbal_genera$var.barplot # 2. bar plot...Enterococcus is the most frequent taxa in CV-process

# making a grid for the variable number plot and bar plots (likely as supplementary data)
library(ggpubr)

ggarrange(selbal_variables, selbal_barplot, 
          ncol = 2,
          nrow = 1,
          labels = c("(A)", "(B)"))



plot.new()
selbal_global <- grid.draw(selbal_genera$global.plot) # 3. default global plot...boxplot, density plot, and AUC-ROC plots

plot.new()
grid.draw(selbal_genera$global.plot2) # a little less informative

plot.new()
plot.tab(selbal_genera$cv.tab) # 4. table summary for CV step

selbal_genera$cv.accuracy # 5. accuracy from CV (length = n.fold * n.iter; 100 here)
summary(selbal_genera$cv.accuracy) # mean cross-validation accuracy is 0.8479 (lower than that in the whole dataset)


selbal_genera$global.balance # 6. variables selected for balance

selbal_genera$glm # 7. regression model summary for balance selection
summary(selbal_genera$glm)

selbal_genera$opt.nvar # 8. number of variables selected for balance analysis (4)

ggarrange(selbal_genera$global.plot,
          ncol = 1,
          nrow = 1)
          #labels = c("(A)", "(B)"))





################################################################################
####################### exporting phyloseq data out of R #########################  
# https://forum.qiime2.org/t/importing-dada2-and-phyloseq-objects-to-qiime-2/4683

ps.genera # 216 taxa and 107 samples

# exporting taxonomy
tax <- as(tax_table(ps.genera), "matrix")

tax_cols <- colnames(tax)

tax <- as.data.frame(tax)

tax$taxonomy <- do.call(paste,
                        c(tax[tax_cols],
                          sep = ";"))
for(co in tax_cols) tax[co] <- NULL

write.table(tax,
            "tax-genera.txt",
            quote = F,
            col.names = F,
            sep = "\t")

# exporting genera feature table
library(biomformat)

otu <- as(otu_table(ps.genera), "matrix") # not transposed cos otus are rows

otu_biom <- make_biom(data = otu)

write_biom(otu_biom,
           "otu_biom.biom")

# exporting feature table as a text file (alternative)
write.table(
  otu_table(ps.genera),
  "otu-genera.txt",
  sep = "\t",
  row.names = T,
  col.names = NA,
  quote = F
)




################################################################################
##################### post network analysis in qiime2 ##########################

# importing files into R for ANCOM-BC analysis on network files

setwd("~/Desktop/Seqs/bacteria/R_analyses/R_to_QIIME2")

list.files("~/Desktop/Seqs/bacteria/R_analyses/R_to_QIIME2")

library(qiime2R) # already loaded for this session
ps.network <- qza_to_phyloseq(
  features = "feature-table-genera-collapsed.qza", 
  taxonomy = "taxonomy-genera.qza",
  metadata = "metadata.tsv"
)

ps.network # 44 taxa by 107 samples


##################### ANCOM-BC on network ps object ######################
library(ANCOMBC)


ancom_net <- ancombc(phyloseq = ps.network, 
                        formula = "Status",
                        p_adj_method = "BH",
                        zero_cut = 0.90, 
                        lib_cut = 1000, 
                        group = "Status", 
                        struc_zero = T, 
                        neg_lb = T, 
                        tol = 1e-5,
                        max_iter = 1000, # 1000 iterations
                        conserve = T,
                        alpha = 0.05, # alpha cut-off value of 5%
                        global = T)

# reformatting ancom results
ancom_net

ancom_net_df <- data.frame(
  Genus = row.names(ancom_net$res$beta),
  beta = unlist(ancom_net$res$beta),
  se = unlist(ancom_net$res$se),
  W = unlist(ancom_net$res$W),
  p_val = unlist(ancom_net$res$p_val),
  q_val = unlist(ancom_net$res$q_val),
  diff_abn = unlist(ancom_net$res$diff_abn)
)




fdr_ancom_net <- ancom_net_df %>%
  dplyr::filter(q_val < 0.05)

dim(fdr_ancom_net) # 15 taxa significant



View(fdr_ancom_net) # feature Id only

# exporting to manually add taxonomy details outside R

write.csv(fdr_ancom_net, "fdr_ancom_net.csv") # edited with the taxonomy-genera.qzv file

fdr_ancom_net_manual_edit <- read.csv("fdr_ancom_net.csv", header = T)

fdr_ancom_net_manual_edit



ancom_net_diff_abun_plot <- ggplot(fdr_ancom_net_manual_edit,
                               aes(x = reorder(Genus, -W),
                                   y = W, 
                                   color = ifelse(W < 0, "#E69F00", "#56B4E9"))) + # colour per group
  geom_point(size = 3) +
  coord_flip() +
 # ggtitle("ANCOM-BC Differential Abundance Test") +
  labs(x = "Genus-level Features", y = "W (AGE Cases vs. Healthy Controls)") +
  geom_hline(yintercept = 0, color = "black", size = 0.5, linetype = "dotted") +
  theme_bw() +
  theme(text = element_text(family = "Helvetica", size = 12)) +
  scale_color_identity()

ancom_net_diff_abun_plot



################################################################################
### splitting phyloseq object per group for network analysis in qiime 2

ps.genera

# healthy controls
ps.healthy <- phyloseq::subset_samples(
  ps.genera, Status == "Healthy")

ps.healthy # 216 taxa by 50 samples

# AGE cases
ps.AGE <- phyloseq::subset_samples(
  ps.genera, Status == "AGE")

ps.AGE # 216 taxa by 57 samples


# exporting taxonomy (heathy ctrls)
tax_healthy <- as(tax_table(ps.healthy), "matrix")

tax_healthy_cols <- colnames(tax_healthy)

tax_healthy <- as.data.frame(tax_healthy)

tax_healthy$taxonomy <- do.call(paste,
                        c(tax_healthy[tax_healthy_cols],
                          sep = ";"))
for(co in tax_healthy_cols) tax_healthy[co] <- NULL

write.table(tax_healthy,
            "tax-healthy-genera.txt",
            quote = F,
            col.names = F,
            sep = "\t")

# exporting genera feature table
library(biomformat)

otu_healthy <- as(otu_table(ps.healthy), "matrix") # not transposed cos otus are rows

otu_healthy_biom <- make_biom(data = otu_healthy)

write_biom(otu_healthy_biom,
           "otu_healthy_biom.biom")



# AGE cases
# exporting taxonomy (AGE)
tax_AGE <- as(tax_table(ps.AGE), "matrix")

tax_AGE_cols <- colnames(tax_AGE)

tax_AGE <- as.data.frame(tax_AGE)

tax_AGE$taxonomy <- do.call(paste,
                                c(tax_AGE[tax_AGE_cols],
                                  sep = ";"))
for(co in tax_AGE_cols) tax_AGE[co] <- NULL

write.table(tax_AGE,
            "tax-AGE-genera.txt",
            quote = F,
            col.names = F,
            sep = "\t")

# exporting genera feature table
library(biomformat)

otu_AGE <- as(otu_table(ps.AGE), "matrix") # not transposed cos otus are rows

otu_AGE_biom <- make_biom(data = otu_AGE)

write_biom(otu_AGE_biom,
           "otu_AGE_biom.biom")


########################################################################################
######## Network visualisation
library(igraph)

list.files(setwd("~/Desktop/Seqs/bacteria/R_analyses/R_to_QIIME2/module-exported-file"))

combined_net <- read.graph("network.gml", format = "gml")

summary(combined_net) # 65 nodes, 213 edges

vcount(combined_net) # 65 nodes
ecount(combined_net) # 213 edges

V(combined_net)$label

plot(combined_net,
     layout = layout.fruchterman.reingold(combined_net),
     vertex.size = 5,
     vertex.label = NA
     )

# NET for AGE cases
setwd("~/Desktop/Seqs/bacteria/R_analyses/R_to_QIIME2/Network_AGE_Healthy/modules-AGE-exported-file")

list.files(setwd("~/Desktop/Seqs/bacteria/R_analyses/R_to_QIIME2/Network_AGE_Healthy/modules-AGE-exported-file"))

AGE_net <- read.graph("network.gml", format = "gml")

summary(AGE_net) # 65 nodes, 213 edges

plot(AGE_net,
     layout = layout.fruchterman.reingold(AGE_net),
     vertex.size = 5,
     vertex.label = NA
)



# Healthy NET
setwd("~/Desktop/Seqs/bacteria/R_analyses/R_to_QIIME2/Network_AGE_Healthy/modules-healthy-exported-file")

list.files(setwd("~/Desktop/Seqs/bacteria/R_analyses/R_to_QIIME2/Network_AGE_Healthy/modules-healthy-exported-file"))

healthy_net <- read.graph("network.gml", format = "gml")

summary(healthy_net) # 65 nodes, 213 edges

plot(healthy_net,
     layout = layout.fruchterman.reingold(healthy_net),
     vertex.size = 5,
     vertex.label = NA
)
