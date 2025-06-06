#load libraries
library(ggplot2)
library(stats)
library(dplyr)

#read in data
larv.dat<-read.csv('shell_data.csv')

#replace a (abnormal) with 0 and n (normal) with 1
larv.dat$normality<-gsub('a', 0, larv.dat$normality)
larv.dat$normality<-gsub('n', 1, larv.dat$normality)
larv.dat$normality<-as.factor(larv.dat$normality)

#generalized linear model of shell normality with possible effects of parent priming, larval treatment, and time
glm.larv<-glm(normality ~ as.factor(B_treat) * as.factor(L_treat) + day + (1|as.numeric(ID)), family=binomial, data=larv.dat)
summary(glm.larv)

#reformat data so normality is a percentage
res_perc<- larv.dat %>%
  group_by(L_treat, B_treat, day) %>%
  summarise(
    total = n(),
    normal_count=sum(normality == 1),
    percent_normal=100*normal_count/total
  )

level_order<-c('C.2', 'C.7', 'C.14', 'T.2', 'T.7', 'T.14')

#bar plot of percent normal shells by parenta and larval treatment
ggplot(data=res_perc, aes(x=factor(interaction(B_treat,day), level=level_order), y=percent_normal, fill=B_treat)) +
  geom_bar(stat='identity', position=position_dodge()) +
  theme_bw() +
  facet_wrap(~L_treat) +
  ylab('Percent Normal Shells') +
  xlab('Broodstock Treatment') +
  scale_fill_manual(values=c('darkblue', 'goldenrod1')) 

#consider size distribution and shell normality
ggplot(data=larv.dat, aes(x=length, fill=interaction(day,B_treat))) +
  geom_density(alpha=0.7) +
  #facet_wrap(~normality+L_treat, scales='free_y', strip.position='right') +
  facet_grid(normality~L_treat)+
  theme_bw() +
  scale_fill_manual(values=c(rep('darkblue',3), rep('goldenrod1',3))) +
  theme(legend.position="none")
