# load libraries
library(car)
library(ggplot2)


# get data
data <-read.csv("condition_index.csv")

# calculate dry weight
#subtract mass of foil from dry mass
data$dry_tissue <- (data$dry - data$foil)

# calculate body condition
data$condition <- (data$dry_tissue/data$shell)*100

# set treatment as factor
data$treatment <- as.factor(data$treatment)


### Testing for differences in means with ANOVA ####
#test for data normality
test_me <- data$condition

#produce q-q plot and est line fit
qqnorm(test_me, main = "Q-Q Plot: untransformed") # check linearity
qqline(test_me)
norm_test <- shapiro.test(test_me) # p-value > 0.05 = good, don't need transformation
print(paste("shapiro test p-value, untransformed:", norm_test$p.value))


### Testing for differences in means with ANOVA ####
#Look at just T1 and T2 to determine effect of OA treatment
data_t1_2 <- subset(data, date != "20623")

# test for homoscedasticity
leveneTest(data_t1_2$condition, data_t1_2$treatment)

#anova test
my_test <- aov(condition ~ treatment * date, data = data_t1_2)
summary(my_test)

#consider sex (T2 only)
data_t2 <- subset(data, date == "40323")

#anova test
my_test2 <- aov(condition ~ treatment * Sex, data = data_t2)
summary(my_test2)

#Boxplot of Condition index by sex and treatment
ggplot(data=data, aes(x=as.factor(date), y=condition, fill=interaction(treatment,Corr.Sex))) +
  geom_boxplot(show.legend=F) +
  theme_bw() +
  ylab("Condition Index") +
  xlab('Sampling Date') +
  scale_x_discrete(labels=c('T0', '1 month', '2 months')) +
  scale_fill_manual(values=c('darkblue', 'goldenrod1','darkblue','goldenrod1','darkblue', 'gray', 'goldenrod1'))
