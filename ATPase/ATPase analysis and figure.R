#load package libraries
library(dplyr)
library(nlme)
library(ggplot2)
library(bestNormalize)
library(car)

atp.dat<-read.csv('ATPase Data (T0 and Tfinal).csv')

#Assess Data normality and heteroscedasticity
shapiro.test(atp.dat$ATPase)
#p<0.05 so cannot assume normality

#transform data to achieve normality
hist(atp.dat$ATPase)
BN.atp<-bestNormalize(atp.dat$ATPase)
#box cox transformation chosen
x<-atp.dat$ATPase
MASS::boxcox(lm(x ~ 1))
#optimal lamda is near 0.5 so do squareroot transformation

atp.dat$atp.trans<-sqrt(atp.dat$ATPase)

shapiro.test(atp.dat$atp.trans)
Shapiro-Wilk normality test
#p-value = 0.7472; can assume normality

#Levene's test
#check by time point
lev.test1<-leveneTest(atp.trans~timepoint, atp.dat)
print(lev.test1)
#p=0.2165 (>0.05) so variance across samples is equal at 0.05 significance

#check by sex
lev.test2<-leveneTest(atp.trans~Sex, atp.dat)
print(lev.test2)
#p=0.3755

#check by treatment
lev.test3<-leveneTest(atp.trans~Treatment, atp.dat)
print(lev.test3)
#p=0.4055

#Linear Mixed Effects Model
#time only, without sex
atp.dat$NewTank<-as.factor(atp.dat$Tank..new.)
ATPase.time<-lme(atp.trans~timepoint, random=~1|NewTank, data=atp.dat)
summary(ATPase.time)

#treatment and sex
atp.dat.Tf<-subset(atp.dat, timepoint=='Tf')
ATPase.timesex<-lme(atp.trans~Treatment*Sex, random=~1|NewTank, data=atp.dat.Tf)
summary(ATPase.timesex)

#plot ATPase activity by OA treatment and sex
ggplot(atp.dat, aes(x=timepoint, y=ATPase, fill=Treatment, color=Sex)) +
  geom_boxplot() +
  theme_bw() +
  ylab(expression(paste('µmol ADP/mg protein ', h^-1))) +
  xlab('Sampling Date') +
  scale_fill_manual(values=c('darkblue', 'gray', 'goldenrod1')) +
  scale_color_manual(values=c(rep('black', 3)))
