## Decontam package to remove contaminants ##
# Parallel analysis to that on my macbook #


# Setting working directory
path <- setwd("C:/Users/Cobb2/Desktop/R_analyses") # on my windows pc

# Listing files in directory
list.files(path)

# Reading metadata and excluding all "NA"
metadata <- read.csv("metadata.csv", 
                     header = T, 
                     na.strings = c("N/A", "", " "))
# NB: second row of metadata with qiime2 format was deleted

View(metadata)

# Loading required packages
library(phyloseq) #1.26.1
library(qiime2R) #0.99.34
library(decontam) #1.2.1
library(ggplot2) #3.3.2
library(tidyverse) #1.3.0
library(microbiome) #1.4.2


# Creating phyloseq object (original code used)

ps <- qza_to_phyloseq(
  features = "filtered-table-bac.qza",
  taxonomy = "taxonomy-bac.qza",
  tree = "insertion-tree-bac.qza",
  metadata = "metadata.tsv"
)

ps


# Viewing the metadata details in the phyloseq object 
head(sample_data(ps))
tail(sample_data(ps))


library(showtext)

font_add(family = "Arial", regular = "Arial.ttf")

font_add(family = "Helvetica", regular = "Helvetica.ttc") # not found

showtext_auto()





## Contaminant removal with decontam ##

# Inspecting library sizes and plotting
df <- as.data.frame(sample_data(ps)) # data frame

df$LibrarySize <- sample_sums(ps)

df <- df[order(df$LibrarySize),]

df$Index <- seq(nrow(df))

ggplot(data = df,
       aes(x = Index, y = LibrarySize, color = Status)) +
  geom_point() +
  xlab("Index") + ylab("Library Size") + # can't verify if it's really Arial here...Modify in Inkscape
  theme_bw() +
  theme(text = element_text(size = 14, family = "Arial"),
        axis.text.x = element_text(size = 12, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        legend.title = element_text(size = 14, family = "Arial"),
        legend.text = element_text(size = 12, family = "Arial"))

# Identifying contaminants using Prevalence. Note the negative controls (NC)
sample_data(ps)$is.neg <- sample_data(ps)$Status == "NC"

contamdf.prev <- isContaminant(ps, 
                               method = "prevalence", 
                               neg = "is.neg")

table(contamdf.prev$contaminant) # number of contaminants was 17 

head(which(contamdf.prev$contaminant))

# Using a stringent threshold of 0.5
contamdf.prev05 <- isContaminant(ps,
                                 method = "prevalence",
                                 neg = "is.neg",
                                 threshold = 0.5)

table(contamdf.prev05$contaminant) # number of contaminants was 51

head(which(contamdf.prev05$contaminant))



## Removing contaminants from 0.5 prevalence object ##

# Have a look at phyloseq obj before removal
ps

ps.noncontam <- prune_taxa(!contamdf.prev05$contaminant, ps)

# noncontaminant phyloseq object
ps.noncontam


plot_bar(ps.noncontam, fill = "Phylum")


# Removing NCs (B represents 'Blank') from feature
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



# Removing PCs from features 
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


ps.noncontam_noctrls



# Calculating relative abundance data

ps.noncontam.noctrls.rel <- microbiome::transform(
  ps.noncontam_noctrls, "compositional")

phyloseq::otu_table(ps.noncontam.noctrls.rel)[1:10, 1:10]


phyloseq::plot_bar(ps.noncontam.noctrls.rel, # another code for plotting stacked bar plots
                   fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum),
           stat = "identity", position = "stack") +
  facet_wrap(~ Status, scales = "free")

# or
plot_bar(ps.noncontam.noctrls.phylum.rel, fill = "Phylum") +
  facet_wrap(~Status, scales = "free_x")


###### Boxplots of taxa ######

# Agglomerating to Phylum level
ps.noncontam.noctrls.phylum <- phyloseq::tax_glom(ps.noncontam_noctrls,
                                                  "Phylum")

phyloseq::taxa_names(ps.noncontam.noctrls.phylum) <- phyloseq::tax_table(ps.noncontam.noctrls.phylum)[, "Phylum"]

phyloseq::otu_table(ps.noncontam.noctrls.phylum)[1:5, 1:5] # viewing first five rows and columns




# Melt and plot
phyloseq::psmelt(ps.noncontam.noctrls.phylum) %>%
  ggplot(data = .,
         aes(x = Status, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  facet_wrap(~ Phylum, scales = "free", nrow = 3) +
  scale_color_viridis_d() +
  theme_bw() +
  theme(text = element_text(size = 14))



#### Trying a slightly different plot
phyloseq::psmelt(ps.noncontam.noctrls.phylum) %>%
  ggplot(data = .,
         aes(x = Status, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Status), height = 0, width = .2) + # this plot looks better
  facet_wrap(~ Phylum, scales = "free") +
  scale_color_viridis_d() +
  theme_bw()



#### Trying to test for phylum group differences

# subsetting groups
HCs <- phyloseq::subset_samples(ps.noncontam.noctrls.phylum,
                                Status == "Healthy")

AGE <- phyloseq::subset_samples(ps.noncontam.noctrls.phylum,
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




######## Testing for taxa differences ########

# Agglomerating to taxa level of choice #

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



library(magrittr)
library(broom)
library(Hmisc)


ps.noncontam.noctrls.phylum.rel <- microbiome::transform(
  ps.noncontam.noctrls.phylum, "compositional")


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


library(dplyr)

View(group_by(ps.phy.rel, Phylum) %>% # getting summary statistics for phyla abundances
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

ps.noncontam.noctrls.class.rel <- microbiome::transform(
  ps.noncontam.noctrls.class, "compositional")


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




###     Alpha diversity estimation     ###
library(vegan)
library(microbiome)


# finding the sequence counts
summary(sample_sums(ps.noncontam_noctrls)) # minimum is 4796 (Sample ID: AGE 33)

# plotting a rarefaction curve
asv_tab <- t(abundances(ps.noncontam_noctrls))


rarecurve(asv_tab,
          step = 1000,
          label = T,
          sample = min(rowSums(asv_tab)))

# rarefying to 10000 sequences
set.seed(9242) # for reproducibility

ps.noncontam_noctrls_rar <- rarefy_even_depth(ps.noncontam_noctrls, #AGE33 sample removed
                                              sample.size = 10000)

# checking phyloseq object to see how much data I have left
ps.noncontam_noctrls_rar #2338 taxa, 106 samples

barplot(sample_sums(ps.noncontam_noctrls_rar)) # all counts are 10,000

# Data frame of select alpha diversity measures

View(sample_data(ps.noncontam_noctrls_rar))

alpha.div <- data.frame(
  "Observed" = estimate_richness(ps.noncontam_noctrls_rar, measures = "Observed"),
  "Shannon" = estimate_richness(ps.noncontam_noctrls_rar, measures = "Shannon"),
  "PD" = picante::pd(samp = data.frame(t(data.frame(otu_table(ps.noncontam_noctrls_rar)))), #phylogenetic diversity
                     tree = phy_tree(ps.noncontam_noctrls_rar))[, 1],
  "Status" = sample_data(ps.noncontam_noctrls_rar)$Status,
  "Age" = sample_data(ps.noncontam_noctrls_rar)$Age_In_Months,
  "Gender" = sample_data(ps.noncontam_noctrls_rar)$Gender,
  "Breastfeeding" = sample_data(ps.noncontam_noctrls_rar)$Breastfeeding,
  "Family_Meal" = sample_data(ps.noncontam_noctrls_rar)$Family_Meal_Intake,
  "Artificial_Milk" = sample_data(ps.noncontam_noctrls_rar)$Artificial_Milk_Intake,
  "Barcode" = sample_data(ps.noncontam_noctrls_rar)$Barcode_Sequence_16S
)

head(alpha.div)
View(alpha.div)


library(RColorBrewer)
# plot of alpha.div measures

alpha.div %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "PD")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "PD"))) %>%
  ggplot(aes(x = Status, y = value, color = Status)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(aes(color = Status), alpha = 1.0, shape = 16, size = 3, height = 0.1, width = 0.3) +
  labs(x = "Status", y = "Alpha Diversity Measure") +
  facet_wrap(~ metric, scales = "free") +
  theme_bw() +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  theme(text = element_text(size = 14, family = "Arial"),
        axis.text.x = element_text(size = 12, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.title.x = element_blank(),
        legend.position = "none") # no legend needed
ggsave("Alpha-diversity-plot.pdf", width = 15, height = 10, dpi = 300)




# Summarise metrics

alpha.div %>% # getting a simple summary metric for alpha diversity measures
  group_by(Status) %>%
  summarise(
    median_observed = median(Observed),
    sd_observed = sd(Observed),
    median_shannon = median(Shannon),
    sd_shannon = sd(Shannon),
    median_pd = median(PD),
    sd_pd = sd(PD)
  )

# two-group test for alpha diversity measures


wilcox.test(data = alpha.div, Observed ~ Status, conf.int = T)

wilcox.test(data = alpha.div, Shannon ~ Status, conf.int = T)

wilcox.test(data = alpha.div, PD ~ Status, conf.int = T)


# Linear model to incorporate likely confounders/covariates???

summary(glm(data = alpha.div, 
            Observed ~ Status + Age + Artificial_Milk + Family_Meal, 
            family = "gaussian")) # only Age was significant with Status , so exclude the others


summary(glm(data = alpha.div, 
            Observed ~ Status + Age, 
            family = "gaussian")) #confounders


summary(lm(data = alpha.div,
           Observed ~ Status + Age + Breastfeeding))


summary(lm(data = alpha.div,
           Observed ~ Status + Age))


summary(lm(data = alpha.div,
           Observed ~ Status * Age))


plot(lm(data = alpha.div,
        Observed ~ Status + Age))



####### Beta diversity ##########

ord.unwt.uni <- ordinate(ps.noncontam_noctrls_rar,
                         "PCoA", 
                         "unifrac", 
                         weighted = FALSE) # unweighted unifrac


# checking eigen values
barplot(ord.unwt.uni$values$Eigenvalues[1:10])

plot_ordination(ps.noncontam_noctrls_rar, 
                ord.unwt.uni, 
                color = "Status") +
  ggtitle("Unweighted UniFrac") + 
  geom_point(size = 3, alpha = 1.0) + 
  stat_ellipse(aes(group = Status),
               linetype = 5) +
  theme_bw() + 
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +  # horse shoe effect observed here
  theme(text = element_text(size = 14, family = "Arial"),
        axis.text.x = element_text(size = 12, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        legend.title = element_text(size = 14, family = "Arial"),
        legend.text = element_text(size = 12, family = "Arial")) # Helvetica font not on this Windows pc. Download neccessary package


# PERMANOVA on Unweighted unifrac distance
unwt.unifrac <- distance(ps.noncontam_noctrls_rar,
                        method = "unifrac",
                        weighted = F)

adonis(unwt.unifrac ~ sample_data(ps.noncontam_noctrls_rar)$Status, 
       permutations = 999)

# adonis with covariates
adonis(unwt.unifrac ~ sample_data(ps.noncontam_noctrls_rar)$Status + sample_data(ps.noncontam_noctrls_rar)$Age_In_Months, 
       permutations = 999)


adonis(unwt.unifrac ~ sample_data(ps.noncontam_noctrls_rar)$Status + 
         sample_data(ps.noncontam_noctrls_rar)$Age_In_Months + 
         sample_data(ps.noncontam_noctrls_rar)$Gender +
         sample_data(ps.noncontam_noctrls_rar)$Breastfeeding +
         sample_data(ps.noncontam_noctrls_rar)$Family_Meal_Intake +
         sample_data(ps.noncontam_noctrls_rar)$Artificial_Milk_Intake, 
       permutations = 999)

  


# Weighted UniFrac with rarified data
ps.noncontam_noctrls_rar


ps.rar.rel <- microbiome::transform(ps.noncontam_noctrls_rar, "compositional") # rarified

ord.wt.uni.rar <- ordinate(ps.rar.rel, 
                       "PCoA", 
                       "unifrac", 
                       weighted = TRUE)

barplot(ord.wt.uni.rar$values$Eigenvalues[1:10])


plot_ordination(ps.rar.rel, 
                ord.wt.uni.rar, 
                color = "Status") + 
  ggtitle("Weighted UniFrac") + 
  geom_point(size = 3, alpha = 1.0) + 
  stat_ellipse(aes(group = Status),
               linetype = 5) +
  theme_bw() + 
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +  
  theme(text = element_text(size = 14, family = "Arial"),
        axis.text.x = element_text(size = 12, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        legend.title = element_text(size = 14, family = "Arial"),
        legend.text = element_text(size = 12, family = "Arial"))


# PERMANOVA on Weighted Unifrac

wt.unifrac <- distance(ps.noncontam_noctrls_rar,
                         method = "unifrac",
                         weighted = T)


adonis(wt.unifrac ~ sample_data(ps.noncontam_noctrls_rar)$Status, 
       permutations = 999)


adonis(wt.unifrac ~ sample_data(ps.noncontam_noctrls_rar)$Status + # with covariates
         sample_data(ps.noncontam_noctrls_rar)$Age_In_Months + 
         sample_data(ps.noncontam_noctrls_rar)$Gender +
         sample_data(ps.noncontam_noctrls_rar)$Breastfeeding +
         sample_data(ps.noncontam_noctrls_rar)$Family_Meal_Intake +
         sample_data(ps.noncontam_noctrls_rar)$Artificial_Milk_Intake, 
       permutations = 999)



# Bray-Curtis 
ord.bray <- ordinate(ps.noncontam_noctrls_rar, # use relative abundance data
                         "PCoA",  
                         "bray") # bray curtis 


# checking eigen values
barplot(ord.bray$values$Eigenvalues[1:10])

plot_ordination(ps.noncontam_noctrls_rar, 
                ord.bray, 
                color = "Status") +
  ggtitle("Bray-Curtis") + 
  geom_point(size = 3, alpha = 1.0) + 
  stat_ellipse(aes(group = Status),
               linetype = 5) +
  theme_bw() + 
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +  
  theme(text = element_text(size = 14, family = "Arial"),
        axis.text.x = element_text(size = 12, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        legend.title = element_text(size = 14, family = "Arial"),
        legend.text = element_text(size = 12, family = "Arial"))


# bray-curtis distance
ps.rar.rel

bray_dist <- distance(ps.rar.rel, method = "bray")

bray_dist <- as.matrix(bray_dist)

head(bray_dist)[, 1:5]

bray_dist %>%  # modified from Pat Schloss video but couldn't use further
  as_tibble(rownames = "Samples") %>%
  pivot_longer(- Samples) %>%
  filter(Samples < name)


sub_dist <- list()

groups_all <- sample_data(ps.rar.rel)$Status 


for(group in levels(groups_all)){
  row_group <- which(groups_all == group)
  sample_group <- sample_names(ps.rar.rel)[row_group]
  sub_dist[[group]] <- bray_dist[sample_group, sample_group]
  sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
}

library(reshape)
library(Matrix)

bray_groups <- melt(sub_dist)

df_bray <- bray_groups[complete.cases(bray_groups), ]

df_bray$L1 <- factor(df_bray$L1, levels = names(sub_dist))

head(df_bray)
tail(df_bray)


ggplot(df_bray,
       aes(x = L1, y = value, color = L1)) +
  geom_boxplot() +
scale_color_manual(values = c("#E69F00", "#56B4E9"))

wilcox.test(df_bray$value ~ df_bray$L1)



# PERMANOVA on Bray-Curtis
bray.dis <- distance(ps.noncontam_noctrls_rar,
                       method = "bray")


adonis(bray.dis ~ sample_data(ps.noncontam_noctrls_rar)$Status, 
       permutations = 999)

adonis(bray.dis ~ sample_data(ps.noncontam_noctrls_rar)$Status + # with covariates
         sample_data(ps.noncontam_noctrls_rar)$Age_In_Months + 
         sample_data(ps.noncontam_noctrls_rar)$Gender +
         sample_data(ps.noncontam_noctrls_rar)$Breastfeeding +
         sample_data(ps.noncontam_noctrls_rar)$Family_Meal_Intake +
         sample_data(ps.noncontam_noctrls_rar)$Artificial_Milk_Intake, 
       permutations = 999)




######################### taxa-specific analysis (request from JGC) ####################################

library(microbiome)
library(magrittr)
library(broom)
library(Hmisc)



ps.noncontam_noctrls

ps.noncontam.noctrls.rel <- microbiome::transform(
  ps.noncontam_noctrls, "compositional")

# creating df so I can manually search for the presence of specific taxa before subsetting
ps.noncontam.noctrls.rel.df <- psmelt(ps.noncontam.noctrls.rel)

View(ps.noncontam.noctrls.rel.df)


# Enterococcus
ps.Enterococcus <- subset_taxa(ps.noncontam.noctrls.rel,
                               Genus == "Enterococcus")

ps.Enterococcus.df <- psmelt(ps.Enterococcus) ## Converting to data frame

wilcox.test(Abundance ~ Status, data = ps.Enterococcus.df, confint = T) # Wilcox test

View(group_by(ps.Enterococcus.df, Genus, Status) %>% # getting summary statistics for status
       summarise(
         count= n(), 
         mean = mean(Abundance), 
         sd = sd(Abundance), 
         median = median(Abundance), 
         IQR = IQR(Abundance),
         min = min(Abundance),
         max = max(Abundance)
       ))


# Streptococcus
ps.Streptococcus <- subset_taxa(ps.noncontam.noctrls.rel,
                               Genus == "Streptococcus")

ps.Streptococcus.df <- psmelt(ps.Streptococcus)

View(ps.Streptococcus.df)

wilcox.test(Abundance ~ Status, data = ps.Streptococcus.df, confint = T) # Wilcox test

View(group_by(ps.Streptococcus.df, Genus, Status) %>% # getting summary statistics for status
       summarise(
         count= n(), 
         mean = mean(Abundance), 
         sd = sd(Abundance), 
         median = median(Abundance), 
         IQR = IQR(Abundance),
         min = min(Abundance),
         max = max(Abundance)
       ))


# Enterobacteriaceae
ps.Enterobacteriaceae <- subset_taxa(ps.noncontam.noctrls.rel,
                                Family == "Enterobacteriaceae")

ps.Enterobacteriaceae.df <- psmelt(ps.Enterobacteriaceae)

View(ps.Enterobacteriaceae.df)

wilcox.test(Abundance ~ Status, data = ps.Enterobacteriaceae.df, confint = T) # Wilcox test

View(group_by(ps.Enterobacteriaceae.df, Family, Status) %>% # getting summary statistics for status
       summarise(
         count= n(), 
         mean = mean(Abundance), 
         sd = sd(Abundance), 
         median = median(Abundance), 
         IQR = IQR(Abundance),
         min = min(Abundance),
         max = max(Abundance)
       ))


# Bacteroides
ps.Bacteroides <- subset_taxa(ps.noncontam.noctrls.rel,
                                Genus == "Bacteroides")

ps.Bacteroides.df <- psmelt(ps.Bacteroides)

View(ps.Bacteroides.df)

wilcox.test(Abundance ~ Status, data = ps.Bacteroides.df, confint = T) # Wilcox test

View(group_by(ps.Bacteroides.df, Genus, Status) %>% # getting summary statistics for status
       summarise(
         count= n(), 
         mean = mean(Abundance), 
         sd = sd(Abundance), 
         median = median(Abundance), 
         IQR = IQR(Abundance),
         min = min(Abundance),
         max = max(Abundance)
       ))


# Granulicatella
ps.Granulicatella <- subset_taxa(ps.noncontam.noctrls.rel,
                              Genus == "Granulicatella")

ps.Granulicatella.df <- psmelt(ps.Granulicatella)

View(ps.Granulicatella.df)

wilcox.test(Abundance ~ Status, data = ps.Granulicatella.df, confint = T) # Wilcox test

View(group_by(ps.Granulicatella.df, Genus, Status) %>% # getting summary statistics for status
       summarise(
         count= n(), 
         mean = mean(Abundance), 
         sd = sd(Abundance), 
         median = median(Abundance), 
         IQR = IQR(Abundance),
         min = min(Abundance),
         max = max(Abundance)
       ))











