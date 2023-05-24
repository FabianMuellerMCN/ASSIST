## ADDITIONAL ANALYSIS FOR REVISION ###
### Fabian Mueller, 2023, 05, 22. ####

library(dplyr)
library(effectsize)
library(ggpubr)
library(nlme)
library(gtools)
library(lme4)
library(scatr)
library(jmv)
library(ggplot2)
library(reshape)
library(gridExtra) 
library(tidyr)
library(stats)
library(emmeans) 
library(Rmisc) 
library(r2glmm)


#---- set system language to "en", clean workspace, scientific notation ----------------------------------------------------
rm(list = ls())
options("scipen"=100, "digits"=5)
options(stringsAsFactors = F)
#options(contrasts=c(unordered="contr.sum", ordered="contr.poly"))


load("/Users/fmueller/Documents/GitHub/ASSIST/DATA/dat_Analyse_PROTOCOLL.RData")  # MAKE SURE TO ADJUST PATH! 
#load("/Users/interface/Documents/GitHub/ASSIST/DATA/dat_Analyse_PROTOCOLL.RData")  # MAKE SURE TO ADJUST PATH! 


# CHECKS ###------------------------------------------------------------------------
noncompleters <-  subset(datAnalyse, Dropout_visit_2 == 1) # dropouts in phase 2!
length(unique(noncompleters$VP_Nr))


TG_noncompleters <- subset(noncompleters, Condition == 1)
#CG_noncompleters <- subset(noncompleters, Condition == 0)
TG_completers <- subset(noncompleters, Condition == 1)

#CHECKS: Good 13 dropouts in study phase 2 in TG and 1 dropout in CG
length(unique(TG_noncompleters$VP_Nr))
length(unique(CG_noncompleters$VP_Nr))



# Compare within treatment group completers vs noncompleters: 

# To compare whether IPQ assessment at timepoint post intervention 1 may account for differences, in who dropped out or not.

# exclude all assessments of IPQ and Stuff at timepoint post 2: 

filt_data <- datAnalyse[datAnalyse$Intervention_Phase != "post2", ]

# filter for Treatmentgroup only
filt_Treatmentgroup <- subset(filt_data, Condition == 1)


# Assuming that the two groups completers and noncompleters are independent, => Welch Two Sample t-test (two-sided (default))
completers <- filter(filt_Treatmentgroup, Dropout_visit_2 == 0)
noncompleters <- filter(filt_Treatmentgroup, Dropout_visit_2 == 1)

# IPQ Scores 

summarySE(filt_Treatmentgroup, measurevar="IPQ_Score", groupvars=c("Condition"), na.rm = TRUE) ### Here reported IPQ mean 
summarySE(filt_Treatmentgroup, measurevar="IPQ_Score", groupvars=c("Dropout_visit_2"), na.rm = TRUE)

shapiro.test(completers$IPQ_Score) 
shapiro.test(noncompleters$IPQ_Score)

ttest_result <- t.test(completers$IPQ_Score, noncompleters$IPQ_Score, alternative = "two.sided")


# USX Scores 
summarySE(filt_Treatmentgroup, measurevar="USX_Score", groupvars=c("Condition"), na.rm = TRUE) ### Here reported mean across both interventionphases 
summarySE(filt_Treatmentgroup, measurevar="USX_Score", groupvars=c("Dropout_visit_2"), na.rm = TRUE)

shapiro.test(completers$USX_Score) 
shapiro.test(noncompleters$USX_Score)

ttest_result <- t.test(completers$USX_Score, noncompleters$USX_Score, alternative = "two.sided")

# Subjective Improvement Angst 
summarySE(filt_Treatmentgroup, measurevar="Angst_Verbesserung", groupvars=c("Condition"), na.rm = TRUE)
summarySE(filt_Treatmentgroup, measurevar="Angst_Verbesserung", groupvars=c("Dropout_visit_2"), na.rm = TRUE)

shapiro.test(completers$Angst_Verbesserung) 
shapiro.test(noncompleters$Angst_Verbesserung)

ttest_result <- t.test(completers$Angst_Verbesserung, noncompleters$Angst_Verbesserung, alternative = "two.sided")


# Subjective Improvement Blickverhalten 
summarySE(filt_Treatmentgroup, measurevar="Blickv_Verbesserung", groupvars=c("Condition"), na.rm = TRUE) ### Here reported IPQ mean 
summarySE(filt_Treatmentgroup, measurevar="Blickv_Verbesserung", groupvars=c("Dropout_visit_2"), na.rm = TRUE)

shapiro.test(completers$Blickv_Verbesserung) 
shapiro.test(noncompleters$Blickv_Verbesserung)

ttest_result <- t.test(completers$Blickv_Verbesserung, noncompleters$Blickv_Verbesserung, alternative = "two.sided")


# Subjective Improvement Performance 
summarySE(filt_Treatmentgroup, measurevar="Improvisation_Verbesserung", groupvars=c("Condition"), na.rm = TRUE)
summarySE(filt_Treatmentgroup, measurevar="Improvisation_Verbesserung", groupvars=c("Dropout_visit_2"), na.rm = TRUE)

shapiro.test(completers$Improvisation_Verbesserung) 
shapiro.test(noncompleters$Improvisation_Verbesserung)

ttest_result <- t.test(completers$Improvisation_Verbesserung, noncompleters$Improvisation_Verbesserung, alternative = "two.sided")
mann_whitney_test <- wilcox.test(completers$Improvisation_Verbesserung, noncompleters$Improvisation_Verbesserung)


######### CHapter 2: Effects of N Training on Delta #######################################################################

### 2.1 Data wrangling ####
# Select Data for Treatment condition only
Treatmentgroup <- subset(datAnalyse, Condition == 1)

# Reduce Dataframe to relevant variables in a way that am able to the convert to wide format (hereby I dont yet take the n trainings into account, cuase it gets ungäbig)
red_dat <- 
  Treatmentgroup %>% 
  select(c("VP_Nr","Intervention_Phase","Angst_im_moment", "Age", "Sex")) #"Angst_Augenkontakt", "relativeDwelltimeAOI"#,"App_ASSIST-nTrainings"))
 # "App_ASSIST_fearMaximumLevel-mean", "App_ASSIST_training_min-mean", "App_ASSIST_level-mean", "App_ASSIST_fear-mean", "Age", "Sex"))
  
# Convert the dataset from long format to wide format
red_dat_wide <- red_dat %>%
  pivot_wider(names_from = Intervention_Phase, values_from = Angst_im_moment)

# Convert the dataset from long format to wide format
red_dat_wide <- red_dat %>%
  spread(key = Intervention_Phase, value = Angst_im_moment)

# Remove rows with NAs using complete.cases()
CC <- red_dat_wide[complete.cases(red_dat_wide), ]
# CHECK its now 28 participants, thats correct!

# Calculate the delta between timepoint 2 and timepoint 3
CC$delta_fear <- CC$post2 - CC$post1


# Now creat a reduced dataframe that contains the n trainings from timepoint 2 to timepoint 3

dat_n <- 
  Treatmentgroup %>% 
  select(c("VP_Nr","Intervention_Phase","App_ASSIST_training_min-mean")) #"Angst_Augenkontakt", "relativeDwelltimeAOI"#,"App_ASSIST-nTrainings"))
# "App_ASSIST_fearMaximumLevel-mean", "App_ASSIST_training_min-mean", "App_ASSIST_level-mean", "App_ASSIST_fear-mean", "Age", "Sex"))


# Filter dataset to include only rows with Intervention_Phase == "post2"
dat_n_post2 <- filter(dat_n, Intervention_Phase == "post2")


# merge the dataframes CC with the dat_n_post2
df <- merge(CC, dat_n_post2, by = "VP_Nr")

dat_fin <- df[, c("VP_Nr", "delta_fear", "App_ASSIST_training_min-mean", "Age", "Sex")]

# rename App_ASSIST-nTrainings to n_train 
colnames(dat_fin)[colnames(dat_fin) == "App_ASSIST_training_min-mean"] <- "n_train"

## Linear Model on the Number of Trainings on Delta Angst

D_train <- lm(delta_fear ~ n_train + Age + Sex, data = dat_fin)
summary(D_train) # n_train: - 
# Pro 1 Training mehr, eine Zunahme von Mood um 0.006


### Ploting a Linear Regression of Number of Trainings on Mood Delta

D_train %>%
  ggplot(aes(x=n_train,y=delta_fear)) +
  geom_point(alpha=0.5) +
  labs(x= "App_ASSIST_training_min-mean", y="Fear Delta")+
  geom_smooth(method=lm)







