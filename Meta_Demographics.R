# changing directory
path <- setwd("~/Desktop/Seqs/bacteria/R_analyses")

# listing files in the directory
list.files(path)

# reading in metadata
metadata <- read.csv("metadata_dem.csv", 
                     header = T, 
                     na.strings = c("N/A", "", " ")) # excluding cells with no data

View(metadata)

# NB: second row of metadata with qiime2 format was deleted


# checking the Gender and finding proportions
table(metadata$Status, metadata$Gender)

library(tidyverse)

metadata %>%
  count(Status, Gender) %>%
  group_by(Status) %>%
  mutate(percent = n/sum(n) * 100)


chisq.test(metadata$Status, 
           metadata$Gender, 
           correct = F) # excluding the Yates' correction
?chisq.test

fisher.test(metadata$Status, metadata$Gender)
 

# checking age of participants
shapiro.test(metadata$Age_In_Months)

# data is non-parameteric (not normally distributed) [i.e. I reject the null hypothesis (data is normally distributed) of the Shapiro-Wilk normality test and accept the alternative (data is not normally distributed), based on p-val < 0.05]

# summary stats for age
library(dplyr) # tidyverse

age <- group_by(metadata, Status) %>% 
  summarise(
    count= n(), 
    mean = mean(Age_In_Months), 
    sd = sd(Age_In_Months), 
    median = median(Age_In_Months), 
    IQR = IQR(Age_In_Months),
    min = min(Age_In_Months),
    max = max(Age_In_Months))

as.data.frame(age)

# analysing age with Mann-Whitney test
wilcox.test(Age_In_Months ~ Status, 
            data = metadata)

# finding levels under vomiting
levels(as.factor(metadata$Vomitting)) 

table(metadata$Status, metadata$Vomitting)

metadata %>%
  #filter(!is.na(Vomitting)) %>% filter excludes missing responses...not needed for summary, as we need to know the # missing
  count(Status, Vomitting) %>%
  group_by(Status) %>%
  mutate(percent = n/sum(n) * 100)

chisq.test(metadata$Status, 
           metadata$Vomitting,
           correct = F)

fisher.test(metadata$Status, metadata$Vomitting)


# finding levels under fever variable
table(metadata$Status, metadata$Fever)

metadata %>%
  #filter(!is.na(Fever)) %>% filter excludes missing responses...not needed for summary, as we need to know the # missing
  count(Status, Fever) %>%
  group_by(Status) %>%
  mutate(percent = n/sum(n) * 100)

chisq.test(metadata$Status, 
           metadata$Fever,
           correct = F)

fisher.test(metadata$Status, metadata$Fever)

# finding levels under mode of feeding
# breastfeeding
levels(as.factor(metadata$Breastfeeding))

table(metadata$Status, metadata$Breastfeeding)

metadata %>%
  #filter(!is.na(Breastfeeding)) %>% filter excludes missing responses...not needed for summary, as we need to know the # missing
  count(Status, Breastfeeding) %>%
  group_by(Status) %>%
  mutate(percent = n/sum(n) * 100)

chisq.test(metadata$Status, 
           metadata$Breastfeeding,
           correct = F)

fisher.test(metadata$Status, metadata$Breastfeeding)


# artificial milk
table(metadata$Status, metadata$Artificial_Milk_Intake)

metadata %>%
  #filter(!is.na(Artificial_Milk_Intake)) %>% filter excludes missing responses...not needed for summary, as we need to know the # missing
  count(Status, Artificial_Milk_Intake) %>%
  group_by(Status) %>%
  mutate(percent = n/sum(n) * 100)

chisq.test(metadata$Status, 
           metadata$Artificial_Milk_Intake,
           correct = F)

fisher.test(metadata$Status, metadata$Artificial_Milk_Intake)

# formula
table(metadata$Status, metadata$Formula_Intake)

metadata %>%
  #filter(!is.na(Formula_Intake)) %>% filter excludes missing responses...not needed for summary, as we need to know the # missing
  count(Status, Formula_Intake) %>%
  group_by(Status) %>%
  mutate(percent = n/sum(n) * 100)

chisq.test(metadata$Status, 
           metadata$Formula_Intake,
           correct = F)

# family meal
table(metadata$Status, metadata$Family_Meal_Intake)

metadata %>%
  #filter(!is.na(Family_Meal_Intake)) %>% filter excludes missing responses...not needed for summary, as we need to know the # missing
  count(Status, Family_Meal_Intake) %>%
  group_by(Status) %>%
  mutate(percent = n/sum(n) * 100)

chisq.test(metadata$Status, 
           metadata$Family_Meal_Intake,
           correct = F)


# Rotavirus vaccination
levels(as.factor(metadata$Rotavirus_Vaccination))

table(metadata$Status, metadata$Rotavirus_Vaccination)

metadata %>%
  #filter(!is.na(Rotavirus_Vaccination)) %>% filter excludes missing responses...not needed for summary, as we need to know the # missing
  count(Status, Rotavirus_Vaccination) %>%
  group_by(Status) %>%
  mutate(percent = n/sum(n) * 100)

chisq.test(metadata$Status, 
           metadata$Rotavirus_Vaccination,
           correct = F)

######## end for summarising demographic data of participants ##########