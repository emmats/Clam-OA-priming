#load libraries
library(ggplot2)
library(nlme)
library(stats)
library(dplyr)

#read in data
larv.dat<-read.csv('shell_data.csv')

#linear mixed effects model to determine effect of parental primin, larval OA treatment, and time on larval size
larv.size<-lme(length~B_treat*L_treat+day, random=~1|ID, data=larv.dat)

summary(larv.size)

#Density plot of larval shell size by treatment; mean shell length vertical lines
ggplot(data=larv.dat, aes(x=length, fill=B_treat)) +
  geom_density(alpha=0.7) +
  theme_bw() +
  stat_summary(aes(xintercept = ..x.., y = 0, color=B_treat), fun = mean, geom = "vline", orientation = "y") +
  facet_wrap(~L_treat+day, scales='free_y') +
  scale_fill_manual(values=c('darkblue', 'goldenrod1')) +
  scale_color_manual(values=c('darkblue', 'goldenrod1'))
