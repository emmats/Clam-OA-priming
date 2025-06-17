# load libraries

library(seacarb)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
library(data.table)
library(lubridate)

# get data

setwd("~/Documents/OAP_shellfish/2023_Manila_chemistry_data/")
data <-fread("updated environmentaldata.csv")

# organize and clean up the data a bit
# 1)insitu temp as numeric, 2)make a new variable to group type and treatment, 
# 3)get date format, 4)remove NA rows (since we can calc with those anyway)

data$insitu_temperature<-as.numeric(data$insitu_temperature)
data$group <- paste0(data$treatment,"_",data$rep_type)
data$date <- mdy(data$date)
data_clean <- na.omit(data)
str(data_clean)

### use seacarb for alkalinity calculation

# calculate dissolved carbon and create new row
data_clean$dic <- carb(flag = 8, data_clean$ph_with_correction, data_clean$alkalinity/1000000, T = data_clean$spectrometery_temperature, S = data_clean$salinity)$DIC * 1000000

#calculate in situ pH
data_clean$insitu_pH <- carb(flag = 15, data_clean$alkalinity/1000000, data_clean$dic/1000000, T = data_clean$insitu_temperature, S = data_clean$salinity)$pH

#calculate dissolved CO2 partial pressure
data_clean$pCO2 <- carb(flag = 15, data_clean$alkalinity/1000000, data_clean$dic/1000000, T = data_clean$insitu_temperature, S = data_clean$salinity)$pCO2

#calculate aragonite saturation
data_clean$aragonite <- carb(flag = 15, data_clean$alkalinity/1000000, data_clean$dic/1000000, T = data_clean$insitu_temperature, S = data_clean$salinity)$OmegaAragonite

#calculate calcite saturation
data_clean$calcite <- carb(flag = 15, data_clean$alkalinity/1000000, data_clean$dic/1000000, T = data_clean$insitu_temperature, S = data_clean$salinity)$OmegaCalcite

#######split data into broodstock and larval treatments
##? should I remove acclimation? 

ADLT<- data_clean %>%
  filter(date > "2023-01-24" & date <"2023-04-18")

LRV<- data_clean %>%
  filter(date > "2023-04-18")

#######summary data for adult table (mean and SD)
#need to subset date to start when TARGET pH was reached (date >)
#pull: pH, salinity and all the calculated stuff
ADLT_post_targetPH<- ADLT %>%
  filter(date >= "2023-02-08")
head(ADLT_post_targetPH)

ADLTsummary<-ADLT %>% 
              group_by(treatment,tank) %>%
              summarise(across(c(insitu_pH, pCO2,salinity,aragonite,dic), 
                  list(mean=mean, sd=sd), na.rm=TRUE))
ADLTsummary

#save the summary
#write.table(ADLTsummary,file = "Manila_2023_adult_chem_summary.txt", sep = '\t')
#######summary data for larval table (mean and SD)
#pull: temp pH, salinity and all the calculated stuff
head(LRV)

LRVsummary<-LRV %>% 
  group_by(treatment,tank) %>%
  summarise(across(c(insitu_temperature,insitu_pH, pCO2,salinity,aragonite,dic), 
                   list(mean=mean, sd=sd), na.rm=TRUE))
LRVsummary
#save the summary
#write.table(LRVsummary,file = "Manila_2023_larval_chem_summary.txt", sep = '\t')

#figures for adult chemistry
a_pH <- ADLT %>% filter(rep_type == "B") %>%
  ggplot(., aes(x=date,insitu_pH,color = group, fill = group)) +
  geom_smooth() + 
  theme(legend.position="none") +
  xlab("date") +
  ylab ("pH")
a_pH

a_pCO2 <- ADLT %>% filter(rep_type == "B") %>%
  ggplot(., aes(x=date,pCO2,color = group, fill = group)) +
  geom_smooth() + 
  scale_y_reverse() +
  theme(legend.position="none") +
  xlab("") +
  ylab ("pCO2")
a_pCO2

a_arag <- ADLT %>% filter(rep_type == "B") %>%
  ggplot(., aes(x=date,aragonite,color = group, fill = group)) +
  geom_smooth() + 
  theme(legend.position="none") +
  xlab("") +
  ylab ("aragonite")
a_arag

a_temp <- ADLT %>% filter(rep_type == "B") %>%
  ggplot(., aes(x=date,insitu_temperature,color = group, fill = group)) +
  geom_smooth() + 
  theme(legend.position="none") +
  xlab("") +
  ylab ("temperature")
a_temp

a_sal <- ADLT %>% filter(rep_type == "B") %>%
  ggplot(., aes(x=date,salinity,color = group, fill = group)) +
  geom_smooth() + 
  ylim(29, 30) +
  theme(legend.position="none") +
  xlab("") +
  ylab ("salinity")
a_sal

#####facet the adult plots - may only need temp and pH for the final TBH
#the other data can be in a summary table#

#grid the plots
a_faceted<-plot_grid(a_temp, a_sal, a_pCO2, a_arag,a_pH, ncol=1, label_size = 12, rel_heights = c(2,2,1,1,1))
a_faceted
# extract the legend from well-formatted mock plot
#mock
newlabels = c(C_B = "ambient", T_B = "primed (lowpH)")
a_mock <- ADLT %>% filter(rep_type == "B") %>%
  ggplot(., aes(x=date,salinity,color = group, fill = group)) +
  geom_smooth() +
  labs(color = "Treatment",fill = "Treatment") +
  scale_color_discrete(labels = newlabels) +
  scale_fill_discrete(labels = newlabels) 
a_mock
# extract a legend that is laid out horizontally
legend <- get_legend(
  a_mock + theme(legend.box.margin = margin(0, 0, 0, 12)))

#plot faceted plot and legend together
plot_grid(a_faceted, legend, ncol = 2, rel_widths = c(3, .4))

#figures for larval chemistry 
#--larval data is pretty straightforward (no temp changes) maybe no figure needed?####
#leaving in script for now, but prob not needed

l_pCO2 <- LRV %>% filter(rep_type == "S") %>%
  ggplot(., aes(x=date,pCO2,color = group, fill = group)) +
  geom_smooth() + 
  scale_y_reverse() +
  theme(legend.position="none") +
  xlab("") 
l_pCO2

l_arag <- LRV %>% filter(rep_type == "S") %>%
  ggplot(., aes(x=date,aragonite,color = group, fill = group)) +
  geom_smooth() +
  theme(legend.position="none") +
  xlab("") 
l_arag

l_pH <- LRV %>% filter(rep_type == "S") %>%
  ggplot(., aes(x=date,insitu_pH,color = group, fill = group)) +
  geom_smooth() +
  xlab("date") +
  theme(legend.position="none") 
l_pH

l_temp <- LRV %>% filter(rep_type == "S") %>%
  ggplot(., aes(x=date,insitu_temperature,color = group, fill = group)) +
  geom_smooth() +
  theme(legend.position="none") +
  xlab("") 
l_temp

l_sal <- LRV %>% filter(rep_type == "S") %>%
  ggplot(., aes(x=date,salinity,color = group, fill = group)) +
  geom_smooth() + 
  ylim(28, 30) +
  theme(legend.position="none") + 
  xlab("") 
l_sal

#facet - update to no title and one legend when you get a minute..
#try this: https://wilkelab.org/cowplot/articles/shared_legends.html
l_faceted<-plot_grid(l_temp, l_sal, l_pCO2, l_arag,l_pH , ncol=1, label_size = 12)
l_faceted

#plot faceted plot and legend together
plot_grid(l_faceted, legend, ncol = 2, rel_widths = c(3, .4))

