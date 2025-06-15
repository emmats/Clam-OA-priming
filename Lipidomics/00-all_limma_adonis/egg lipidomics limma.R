#load libraries
library(limma)
library(tidyr)
library(ggplot2)
library(reshape)
source('biostats.R')

############lipid class##########
class.quant<-read.csv('class quant.csv', row.names=1)
class.t<-data.frame(t(class.quant))
class.t$prepared_blank[is.na(class.t$prepared_blank)] <- 0

#cut-off 80% missingness for samples
class.rmv.msg <- class.t[rowSums(is.na(class.t[,2:31])) < ncol(class.t)*.8,]

#subtract the values in first column (blanks) from values in all other columns
class.rmv.msg[-1] <- lapply(class.rmv.msg[-1], function(x) x - class.rmv.msg$prepared_blank)

#set negative quant values to 0
class.rmv.msg[class.rmv.msg < 0] <- 0

#Normalize?
class.melt<-melt(class.rmv.msg[-1])
ggplot(data=class.melt, aes(x=variable, y=value))+
  geom_boxplot()

# Compute sample medians (row-wise) 
sample_medians <- apply(class.rmv.msg[-1], 2, median, na.rm = TRUE)

# Normalize by sample median
class_normalized <- sweep(class.rmv.msg[-1], 2, sample_medians, '/')

boxplot(class_normalized)

hist(t(class_normalized), breaks=100)

#not normally distributed. Log transform.
class.tra<-data.trans((class_normalized), method='log', plot=F)

boxplot(class.tra)
hist(t(class.tra))

#design matrix
treat.dat<-read.csv('lipidomics treatment.csv')
treat.dat$XName<-colnames(class.rmv.msg[-1])
treat.dat<-treat.dat[-1]
colnames(treat.dat)<-c('Condition', 'Replicate')

design.treat <- model.matrix(~Condition, data = treat.dat)
colnames(design.treat) = gsub("Condition", "", colnames(design.treat))

fit.class <- limma::lmFit(class.tra, design = design.treat)

fit.class1 <- limma::eBayes(fit.class, robust = TRUE, trend = TRUE)

diffTab.class = limma::topTable(fit.class1,
                                  n=Inf, adjust = "fdr", p.value = 0.05, sort.by = "B")
diffTab.class.p1 = limma::topTable(fit.class1,
                                n=Inf, adjust = "fdr", p.value = 1, sort.by = "B")
write.csv(diffTab.class.p1, "lipidomics/Limma output lipid classes p1.csv", quote=F)
#0

volcanoplot(fit.class1, coef=2, main='Lipid Classes')

###########Lipid Species##########
species.quant<-read.csv('species quant.csv', row.names=1)
species.t<-data.frame(t(species.quant))
species.t$prepared_blank[is.na(species.t$prepared_blank)] <- 0

#cut-off 80% missingness for samples
species.rmv.msg <- species.t[rowSums(is.na(species.t[,2:31])) < ncol(species.t)*.8,]

#subtract the values in first column (blanks) from values in all other columns
species.rmv.msg[-1] <- lapply(species.rmv.msg[-1], function(x) x - species.rmv.msg$prepared_blank)

#set negative quant values to 0
species.rmv.msg[species.rmv.msg < 0] <- 0

#Normalize?
boxplot(species.rmv.msg[-1])

# Compute sample medians (row-wise) 
sample_medians2 <- apply(species.rmv.msg[-1], 2, median, na.rm = TRUE)

# Normalize by sample median 
species_normalized <- sweep(species.rmv.msg[-1], 2, sample_medians2, '/')

boxplot(species_normalized)

hist(t(species_normalized), breaks=100)

#not normally distributed. Log transform.
sp.tra<-data.trans((species_normalized), method='log', plot=F)

boxplot(sp.tra)
hist(t(sp.tra))

fit.sp <- limma::lmFit(sp.tra, design = design.treat)

fit.sp1 <- limma::eBayes(fit.sp, robust = TRUE, trend = TRUE)

diffTab.sp = limma::topTable(fit.sp1,
                                n=Inf, adjust = "fdr", p.value = 0.05, sort.by = "B")

diffTab.sp.p1 = limma::topTable(fit.sp1,
                                n=Inf, adjust = "fdr", p.value = 1, sort.by = "B")
write.csv(diffTab.sp.p1, "lipidomics/Limma output lipid species p1.csv", quote=F)
#0
volcanoplot(fit.sp1, coef=2, main='Lipid Species')

#Lipid species volcano plot with colors
#read in data
lipsp.annot<-read.csv('lipid species categories.csv')
diff.sp.annot<-merge(x=diffTab.sp.p1, y=lipsp.annot, by.x='row.names', by.y='Species.1', all.x=T)

#assigned broad category colors to lipid class
lipid.colors<-c("Cholesteryl Ester"= 'black', "Ceramide" = 'purple3', "Diglyceride"='tan1', "Fatty acid" = 'violetred1', "Hexosylceramide" = 'plum', "Lysophosphatidylcholine"='dodgerblue3', "Lysophosphatidylethanolamine" = 'springgreen4', "Lysohposphatidylglycerol" = 'gold3', "Lysophosphatidylinositol" = 'darkorange',"Lipopolysaccharide" = 'red', "Palmitic acid"= 'violetred3', "Phosphatidylcholine" = 'dodgerblue', "Phosphatidylethaolamine" = 'mediumspringgreen', "Ether-linked phosphatidyl-ethanolamine" = 'springgreen3',"Phosphatidylethaolamine plasmalogen" = 'springgreen2', "Phosphatidylglycerol" = 'gold',"Phosphatidylinositol" = 'orange', "Phosphatidylserine"= 'pink', "Sphingomyelin" = 'turquoise1', "Triglyceride"='sienna3' )

ggplot(dat=diff.sp.annot) +
  geom_point(aes(x=logFC, y=-log(P.Value), color=Category), alpha=0.5) +
  theme_bw() +
  scale_color_manual(values=lipid.colors) +
  labs(x='Log Fold Change', y='-Log(p-value)')

############Fatty Acids############
FA.quant<-read.csv('FA quant.csv', row.names=1)
FA.t<-data.frame(t(FA.quant))
FA.t$prepared_blank[is.na(FA.t$prepared_blank)] <- 0

#cut-off 80% missingness for samples
FA.rmv.msg <- FA.t[rowSums(is.na(FA.t[,2:31])) < ncol(FA.t)*.8,]

#subtract the values in first column (blanks) from values in all other columns
FA.rmv.msg[-1] <- lapply(FA.rmv.msg[-1], function(x) x - FA.rmv.msg$prepared_blank)

#set negative quant values to 0
FA.rmv.msg[FA.rmv.msg < 0] <- 0

#Normalize?
boxplot(FA.rmv.msg[-1])

# Compute sample medians (row-wise) 
sample_medians3 <- apply(FA.rmv.msg[-1], 2, median, na.rm = TRUE)

# Normalize by sample median 
FA_normalized <- sweep(FA.rmv.msg[-1], 2, sample_medians3, '/')

boxplot(FA_normalized)

hist(t(FA_normalized), breaks=100)

#not normally distributed. Log transform.
FA.tra<-data.trans((FA_normalized), method='log', plot=F)

boxplot(FA.tra)
hist(t(FA.tra))

fit.fa <- limma::lmFit(FA.tra, design = design.treat)

fit.fa1 <- limma::eBayes(fit.fa, robust = TRUE, trend = TRUE)

diffTab.fa = limma::topTable(fit.fa1,
                             n=Inf, adjust = "fdr", p.value = 0.05, sort.by = "B")

diffTab.fa.p1 = limma::topTable(fit.fa1,
                             n=Inf, adjust = "fdr", p.value = 1, sort.by = "B")
write.csv(diffTab.fa.p1, "lipidomics/Limma output fatty acids p1.csv", quote=F)

#0
volcanoplot(fit.fa1, coef=2, main='Fatty Acids')

#Fatty acids volcano plot with colors
#read in data
fa.annot<-read.csv('FA classification and dbl bonds.csv')
diff.fa.annot<-merge(x=diffTab.fa.p1, y=fa.annot, by.x='row.names', by.y='MatchFA', all.x=T)

#assigned broad category colors to lipid class
lipid.colors<-c("Cholesteryl Ester"= 'black', "Ceramide" = 'purple3', "Diglyceride"='tan1', "Fatty acid" = 'violetred1', "Hexosylceramide" = 'plum', "Lysophosphatidylcholine"='dodgerblue3', "Lysophosphatidylethanolamine" = 'springgreen4', "Lysohposphatidylglycerol" = 'gold3', "Lysophosphatidylinositol" = 'darkorange', "Lipopolysaccharide" = 'red', "Palmitic acid"= 'violetred3', "Phosphatidylcholine" = 'dodgerblue', "Phosphatidylethaolamine" = 'mediumspringgreen', "Phosphatidylglycerol" = 'gold',"Phosphatidylinositol" = 'orange', "Phosphatidylserine"= 'pink', "Triglyceride"='sienna3' )

ggplot(dat=diff.fa.annot) +
  geom_point(aes(x=logFC, y=-log(P.Value), color=Category), alpha=0.5) +
  theme_bw() +
  scale_color_manual(values=lipid.colors) +
  labs(x='Log Fold Change', y='-Log(p-value)')

  
