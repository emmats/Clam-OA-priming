#load libraries
library(ggplot2)
library(reshape)
library(nlme)

#read in data
surv.dat<-read.csv('survival_data.csv')

#@ 0 hpf put approximately same number of unfertilized eggs in each container, but eggs were no accurately counted
#@ 48 hpf, "out" is a rough approximate of fertilization success

#create a metric of parent x larval treatment
surv.dat$combined.treatment<-paste(surv.dat$brood_treatment, surv.dat$larvae_treatment)

#Calculate percent survival and 7 and 14 dpf
surv.dat$perc.surv7dpf<-surv.dat$X7d_out/surv.dat$X7d_in
surv.dat$perc.surv14dpf<-surv.dat$X14_out/surv.dat$X14d_in

surv.dat[is.na(surv.dat)] <- 0

#LME to determine effect of treatment at days 7 and 14
#largest correlative effects were from silo, not brood group

#lme for 7 dpf
larv.surv7<-lme(perc.surv7dpf~brood_treatment*larvae_treatment, random = ~1|silo, data=surv.dat, na.action=na.omit)
summary(larv.surv7)

surv.dat$perc.surv14dpf[is.na(surv.dat$perc.surv14dpf)] <- 0

#lme for 14dpf
larv.surv14<-lme(perc.surv14dpf~brood_treatment*larvae_treatment, random = ~1|silo, data=surv.dat, na.action=na.omit)
summary(larv.surv14)

#Reaction norm plots of survival by parent and larval treatments
#day 7 post-fertilization

surv.7dpf<-subset(surv.dat, select=c(brood_treatment, larvae_treatment, perc.surv7dpf))

melt.7dpf<-melt(surv.7dpf, id.vars=c('brood_treatment', 'larvae_treatment'))

ggplot(data=melt.7dpf, aes(x=larvae_treatment, y=value)) +
  stat_summary(aes(group=brood_treatment), fun=mean,geom='path')+
  stat_summary(aes(color=brood_treatment), fun.data=mean_cl_boot, geom='errorbar', width=0.1)+
  stat_summary(aes(color=brood_treatment), fun=mean, geom='point', size=4)+
  theme_bw() +
  scale_color_manual(values=c('darkblue', 'goldenrod1')) +
  xlab('Larval Treatment') +
  ylab('Percent Survival of Larvae at 7dpf')

#day 14 post-fertilization
surv.14dpf<-subset(surv.dat, select=c(brood_treatment, larvae_treatment, perc.surv14dpf))

melt.14dpf<-melt(surv.14dpf, id.vars=c('brood_treatment', 'larvae_treatment'))

ggplot(data=melt.14dpf, aes(x=larvae_treatment, y=value)) +
  stat_summary(aes(group=brood_treatment), fun=mean,geom='path')+
  stat_summary(aes(color=brood_treatment), fun.data=mean_cl_boot, geom='errorbar', width=0.1)+
  stat_summary(aes(color=brood_treatment), fun=mean, geom='point', size=4)+
  theme_bw() +
  scale_color_manual(values=c('darkblue', 'goldenrod1')) +
  xlab('Larval Treatment') +
  ylab('Percent Survival of Larvae at 14dpf')
