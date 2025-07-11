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

#plot of egg size by individual clam and treatment
ggplot(data=egg.dat, aes(x=interaction(as.factor(Clam_Num), Treatment), y=Area, fill=Treatment)) +
  geom_boxplot() +
  theme_bw() +
  geom_hline(yintercept=3755.66, color='darkblue', linetype='dashed')+
  geom_hline(yintercept=3686.09, color='goldenrod1', linetype='dashed')+
  scale_fill_manual(values=c('darkblue', 'goldenrod1')) +
  scale_x_discrete(labels = c('199', '211', '218', '226', '310', '319', '334', '341', '363', '376', '460', '482', '488', '7', '24', '30', '33', '43', '44', '83', '88', '383')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab('Clam ID') +
  ylab('Oocyte Surface Area')
