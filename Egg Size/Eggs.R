#Load libraries
library(ggplot2)
library(nlme)

#read in data
egg.dat<-read.csv('clam egg size.csv')

#linear mixed effects model of the effect of OA treatment on egg surface area (individual eggs)
egg.size<-lme(Area~as.factor(Treatment), random=~1|Clam_Num, data=egg.dat, na.action=na.omit)

summary(egg.size)

#consider egg size as average per spawner
egg.sub<-subset(egg.dat, select=c(Clam_Num, Area, Treatment))
egg.avg<-aggregate(Area~Clam_Num, data=egg.sub, FUN=mean)
egg.treatment<-unique(subset(egg.sub, select=c(Clam_Num, Treatment)))
egg.avg2<-merge(x=egg.avg, y=egg.treatment, by='Clam_Num', all.x=T)

#lme (average egg size per female)
egg.size2<-lme(Area~as.factor(Treatment), random=~1|Clam_Num, data=egg.avg2, na.action=na.omit)

summary(egg.size2)

#histogram of egg surface area by treatment
ggplot(data=egg.dat, aes(x=Area, fill=Treatment)) +
  geom_histogram(alpha=0.7, binwidth=20) +
  geom_vline(xintercept=3755.66, color='turquoise')+
  geom_vline(xintercept=3686.09, color='coral')+
  theme_bw() +
  scale_fill_manual(values=c('turquoise', 'coral')) +
  xlab('Oocyte Surface Area') +
  ylab('Frequency of Oocyte Size')

