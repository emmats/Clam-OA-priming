#packages
  library(dplyr)
  library(tidyverse)
  library(data.table)
  library(ggplot2)
  library(RColorBrewer)
  library(cowplot)
  library(ggpubr)
  library(rstatix)
  library(paletteer)


#read in the count table from STAR (or DESEq2)<-ask Giles
#which are the mitochondrial genes? - it's notated by the GeneID (MT start with BJM09_g.. p = protein coding, r = ribo, t = tRNA)

counts<-fread("Transcriptomics/data/STAR_count_data.csv")
head(counts)

#correct the treatment for sample 193 (states T (low pH treatment), but should be C (ambient pH))
names(counts)[names(counts) == "M-T-193"] <- "M-C-193"
head(counts)

#get the total counts for each sample
counts.total <- colSums(counts[,2:31])
head(counts.total)

#sum the counts for each sample by "Type"
counts.agg <- aggregate(counts[,2:31], by=list(counts$Type), FUN=sum)
head(counts.agg)

counts.agg_per <- t(t(counts.agg[,2:31])/counts.total)
row.names(counts.agg_per) <- counts.agg[,1]

#You can check this by doing the following

x <- counts.agg[,2]
y <- counts.total[1]
x/y

head(counts.agg_per)

#transpose the table

t.counts.agg_per<-t(counts.agg_per)
t.counts.agg_per<-as.data.frame(t.counts.agg_per)
head(t.counts.agg_per)
nrow(t.counts.agg_per)

#slice the ID column to get treatment in a separate column
#1. make sample ID a new column
t.counts.agg_per <- tibble::rownames_to_column(t.counts.agg_per, "sampleID")
#2. split into multiple columns using hyphen as delimiter
t.counts.agg_per<-t.counts.agg_per %>% separate_wider_delim(sampleID, "-", names = c("species", "Tx","ID"))
head(t.counts.agg_per)
#fwrite(t.counts.agg_per,"Transcriptomics/output/RNAclassification_compare_distributions_updated050925/Type_breakdown.csv")
#get a mean or median per treatment for each Type
t.counts.agg_per_mean<-t.counts.agg_per %>% 
  group_by(Tx) %>%
  summarise(lncRNA=mean(lncRNA),ncRNA=mean(ncRNA),
            protein_coding=mean(protein_coding),rRNA=mean(rRNA),
            snoRNA=mean(snoRNA),snRNA=mean(snRNA),
            transcribed_pseudogene=mean(transcribed_pseudogene), tRNA=mean(tRNA))
head(t.counts.agg_per_mean)
#visualize
head(t.counts.agg_per_mean)
flip <- setNames(data.frame(t(t.counts.agg_per_mean[,-1])), t.counts.agg_per_mean[,1])
head(flip)
colnames(flip)<-c("ambient","low pH")
head(flip)
flip <- tibble::rownames_to_column(flip, "Category")
head(flip)
tibble::rownames_to_column(flip, "Category")

flip_long <- flip %>%
  pivot_longer(cols = c("ambient","low pH"), 
               names_to = "Treatment", 
               values_to = "Proportion")
head(flip_long)
ggplot(flip_long, aes(x = Treatment, y = Proportion, fill = Category)) +
  geom_bar(stat = "identity") +
  theme(text = element_text(size = 20)) +
  labs(x = "Category", y = "Proportion", title = "Stacked Bar Plot of Treatments") +
  scale_fill_brewer(palette = "Dark2")

###############################################################let's look at this grouped by MT or nuclear
head(counts)

#need to add column that is "mt_orNot" (mitochondrial)
#which are the mitochondrial genes? - it's notated by the GeneID (MT start with BJM09_g.. p = protein coding, r = ribo, t = tRNA)#
#make a new column representing if the gene is on the  MT (mitochondrial) genome or the nuclear genome based on the GeneID
counts <-counts %>% mutate(mt_orNot = ifelse(grepl("^BJM09",GeneID),"MT","nuclear"))

#sum the counts for each sample by mitochondrial or not
counts.agg <- aggregate(counts[,2:31], by=list(counts$mt_orNot), FUN=sum)
head(counts.agg)

counts.agg_per <- t(t(counts.agg[,2:31])/counts.total)
row.names(counts.agg_per) <- counts.agg[,1]

#You can check this by doing the following

x <- counts.agg[,2]
y <- counts.total[1]
x/y

head(counts.agg_per)


#transpose the table

t.counts.agg_per<-t(counts.agg_per)
t.counts.agg_per<-as.data.frame(t.counts.agg_per)
head(t.counts.agg_per)
nrow(t.counts.agg_per)

#slice the ID column to get treatment in a separate column
#1. make sample ID a new column
t.counts.agg_per <- tibble::rownames_to_column(t.counts.agg_per, "sampleID")
#2. split into multiple columns using hyphen as delimiter
t.counts.agg_per<-t.counts.agg_per %>% separate_wider_delim(sampleID, "-", names = c("species", "Tx","ID"))
head(t.counts.agg_per)
#fwrite(t.counts.agg_per,"Transcriptomics/output/RNAclassification_compare_distributions_updated050925/MTorNuclear.csv")

#get a mean or median per treatment for each category
t.counts.agg_per %>% 
  group_by(Tx) %>%
  summarise(MT=mean(MT),nuclear=mean(nuclear))

#proportion MT
head(t.counts.agg_per)

ggplot(t.counts.agg_per, aes(fill=Tx, y=MT, x=Tx)) + 
  geom_boxplot() +
  theme(text = element_text(size = 20)) +
  labs(x = "Category", y = "Proportion", title = "Stacked Bar Plot of Treatments") +
  scale_fill_manual(values=c("darkblue", "goldenrod1"))
 
#outliers?
t.counts.agg_per %>%
  group_by(Tx) %>%
  identify_outliers(MT)
#normality?
t.counts.agg_per %>%
  group_by(Tx) %>%
  shapiro_test(MT)
#qqplot
ggqqplot(t.counts.agg_per, x = "MT", facet.by = "Tx")
#equal variances
t.counts.agg_per %>% levene_test(MT ~ Tx)
#Welch
stat.test <- t.counts.agg_per %>% 
  t_test(MT ~ Tx) %>%
  add_significance()
stat.test
#Students
stat.test2 <- t.counts.agg_per %>%
  t_test(MT ~ Tx, var.equal = TRUE) %>%
  add_significance()
stat.test2

####################################what does the distribution looks like inside mito?

#let's cut the table to only include MT
head(counts)

counts.MT<-counts %>% filter(mt_orNot == "MT")
head(counts.MT)
#get the total counts for each sample
counts.total <- colSums(counts.MT[,2:31])
head(counts.total)

#update the 2 rRNA genes to specify either 12s or 16s( BJM09_gr01 is 12S, BJM09_gr01 is 16S)
counts.MT$Type[14]="rRNA_12S"
counts.MT$Type[15]="rRNA_16S"


#sum the counts for each sample by "Type"
counts.agg <- aggregate(counts.MT[,2:31], by=list(counts.MT$Type), FUN=sum)
head(counts.agg)

counts.agg_per <- t(t(counts.agg[,2:31])/counts.total)
row.names(counts.agg_per) <- counts.agg[,1]


#transpose the table

t.counts.agg_per<-t(counts.agg_per)
t.counts.agg_per<-as.data.frame(t.counts.agg_per)
head(t.counts.agg_per)
nrow(t.counts.agg_per)

#slice the ID column to get treatment in a separate column
#1. make sample ID a new column
t.counts.agg_per <- tibble::rownames_to_column(t.counts.agg_per, "sampleID")
#2. split into multiple columns using hyphen as delimiter
t.counts.agg_per<-t.counts.agg_per %>% separate_wider_delim(sampleID, "-", names = c("species", "Tx","ID"))
head(t.counts.agg_per)


#get a mean or median per treatment for each category
t.counts.agg_per_mean<-t.counts.agg_per %>% 
  group_by(Tx) %>%
  summarise(protein_coding=mean(protein_coding),rRNA_12S=mean(rRNA_12S), rRNA_16S=mean(rRNA_16S), tRNA=mean(tRNA))
head(t.counts.agg_per_mean)
View(t.counts.agg_per_mean)
