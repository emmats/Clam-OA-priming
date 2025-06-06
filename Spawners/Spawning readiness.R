#load libraries
library(tidyverse)
library(summarytools)
library(ggplot2)
library(GGally)

## input data  ##
data <-read.csv("spawned_2023_mockIDs_chisq.csv")


## reformat ##

colnames(data)[1] <- gsub('^...','',colnames(data)[1])


## separate out males and females  ##

manila_f <- subset(data, sex == "f")
manila_m <- subset(data, sex == "m")


## run chi squared test on each ##

chisq.test(manila_f$treatment, manila_f$spawned)
#Pearson's Chi-squared test with Yates' continuity correction

chisq.test(manila_m$treatment, manila_m$spawned)
#Pearson's Chi-squared test with Yates' continuity correction

#Stacked bar chart showing proportion spawners by sex and treatment
spawn.dat2<-read.csv('MAC SPAWN DAT.csv')
level_order<-c('Females.C', 'Females.OA', 'Males.C', 'Males.OA')
  
  ggplot(spawn.dat2, aes(x=interaction(Sex,Treatment), y=spawn, fill=Treatment, alpha=Spawn))+
    geom_bar(stat='identity',position='fill') +
    theme_bw() +
    scale_fill_manual(values=c('darkblue', 'goldenrod1')) +
    facet_wrap(~Sex, scales='free_x') +
    scale_alpha_manual(values=c(0.7, 1))+
    ylab('Proportion Spawners') +
    xlab('Sex and Treatment')
