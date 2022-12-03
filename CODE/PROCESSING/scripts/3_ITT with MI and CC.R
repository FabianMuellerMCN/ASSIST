#### Intention to Treat Analyses ####################################################################################################

# Fabian - 2022.12.03

#---- set system language to "en", clean workspace, scientific notation -----------------------------------------------------------
rm(list = ls())
#options("scipen"=100, "digits"=5)
options(stringsAsFactors = F)
#options(contrasts=c(unordered="contr.sum", ordered="contr.poly"))

#---- load all the required packages ----------------------------------------------------------------------------------------------
library(reshape)
library(tidyr) # short-to-long format manipulations
library(shape)
library(lsr) # to calculate effect sizes out of t-tests
library(r2glmm) #  to obtain estimates for specific fixed effects
library(ggpubr) # to arrange multiple plots with one common legend
library(Rmisc)
library(nlme)
library(car)
library(plyr)
library(emmeans) #for adjusted means
library(magrittr) #for confidence intervals of adjusted means group comparision
library(misty)
library(dplyr)
library(mice)
library(broom.mixed)

#---- Import Data relevant for the ITT (dat.Analyse) --------------------------------------------------------------------------------------
load("/Users/fmueller/Documents/GitHub/ASSIST/DATA/dat_Analyse_PROTOCOLL.RData")  # MAKE SURE TO ADJUST PATH! 



#---- 1.0 add Variable: "baseline measurement of fear" to dataframe # to control for baseline at PostHoc Tests ----------------------------
datPre <- datAnalyse[datAnalyse$Intervention_Phase == "pre",]
v1 = "Angst_im_moment"
c1 = "VP_Nr"
mean_baseline_fear  <- aggregate(datPre[v1],by=datPre[c1],FUN=mean, na.rm=TRUE, na.action=na.omit)
names(mean_baseline_fear) [which(names(mean_baseline_fear) %in% "Angst_im_moment")] <- "Angst_im_moment_baseline" #replace name
datAnalyse <-merge(datAnalyse, mean_baseline_fear, by = "VP_Nr")

#Reduce the dataset to the required data for the main analysis for ITT
redDat <- datAnalyse  %>% select(Angst_im_moment, Sex, Age, Condition, Intervention_Phase, VP_Nr, Angst_im_moment_baseline)

#---- 2.0 MCAR-Test: test if NA's are completely at random: if sign => suggests they are not completely at random  -----------------------------------------------------------------

na.test(redDat, digits = 2, p.digits = 3, as.na = NULL, check = TRUE, output = TRUE)


#---- 3.0 ITT Analysis based on Multiple Imputations ----------------------------------------------------------------------------------------------------------------------

## 3.1 Determine m für mice function: How many NA in % of my data? Based on White's suggestion > % of NA's
totalobs = nrow(redDat)
print("Total number of cells ")
print(totalobs)

# calculating the number of cells with na
missingobs = sum(is.na(redDat$Angst_im_moment))
print("Missing value observations")
print(missingobs)

# calculating percentage of missing values
percentage = (missingobs * 100 )/(totalobs)
print("Percentage of missing values' cells")
print (percentage)

# CHECK: 
#sum(is.na(redDat))

# 3.2 Impute 10 datasets using the Mice function ------------------------------------------------------------------------------------------------------------------
MI_Data <- mice(redDat, m = 10, maxit = 20, seed = 420) # 

MI_Data_com <- complete(MI_Data,"long", inc = F)

# To be able to average, all strings need to be numeric: 
MI_Data_com$Intervention_Phase <- revalue(MI_Data_com$Intervention_Phase, c("pre"="1", "post1"="2", "post2"="3"))
MI_Data_com$Intervention_Phase <- as.numeric(MI_Data_com$Intervention_Phase)
MI_Data_com$VP_Nr <- as.numeric(MI_Data_com$VP_Nr)

# Average the models using the mean function across the 10 imputations as grouped by the number of the .id (which is the id of the row (same across each dataset))
model_MI_average <- MI_Data_com  %>% group_by(-.id) %>%    
  summarise_all(.funs = mean) %>%
  select(-.id, -.imp)                                     # .id, .imp # remove meta-info columns with - select

# refactor what needs to be a factor for the model 
model_MI_average$Intervention_Phase <- as.factor(model_MI_average$Intervention_Phase)
model_MI_average$VP_Nr <- as.factor(model_MI_average$VP_Nr)
model_MI_average$Condition <- as.factor(model_MI_average$Condition)
model_MI_average$Sex <- as.factor(model_MI_average$Sex)

# 3.3 #calculate the main model based on the averaged dataset ------------------------------------------------------------------------------------------------------------------

lme_ave<- lme(Angst_im_moment ~  Age + factor(Sex) + factor(Condition) * factor(Intervention_Phase), random=~1|VP_Nr, data=model_MI_average,control=lmeControl(opt="optim"), na.action="na.omit",method = "REML" )
MI_ave <- anova(lme_ave, type = "marginal")
Res_ave <- as.data.frame(Anova(lme_ave, type = "II"))                   
sum_MI <- summary(lme_ave)

#---- 4.0 POST HOC TESTS  ---------------------------------------------------------------------------------------------------------------------------------------------------------------

# 4.1 MI PostHoc Test at T2; Intervention_Phase == post2 (HomeTreatment) --------------------------------------------------------------------------------------
MIpost2 <- model_MI_average[model_MI_average$Intervention_Phase == "3",]
MIpost2.model<- lme(Angst_im_moment ~ factor(Condition) + Age + factor(Sex) + Angst_im_moment_baseline, random=~1|VP_Nr, data=MIpost2, na.action= "na.omit", control=lmeControl(opt="optim"))

MI_res_post2.model_Anova <- Anova(MIpost2.model, type = "II")                   
MI_res_post2.model_Anova <- as.data.frame(MI_res_post2.model_Anova) 
MI_res_post2.model_summary <- summary(MIpost2.model)
MI_res_post2.model_summary <- MI_res_post2.model_summary[[20]]
MI_res_tValuePost2 <- MI_res_post2.model_summary[2,4]

# #---- get adjusted means for group difference Home ----
emm <- emmeans(MIpost2.model, specs = revpairwise ~ as.factor(Condition),adjust = "none",type = "response") #adjust for changing type of multiple comparison adjustment
emmgroup <- as.data.frame(emm)
emm2 <- emm$contrasts %>% summary(infer = TRUE) 

# Cohen's D T2
used_N <- length(MIpost2$VP_Nr)
used_N_c <- length(which(MIpost2$Condition == 0))
used_N_t <- length(which(MIpost2$Condition == 1))
tval <- MI_res_post2.model_summary["factor(Condition)1", "t-value"]

cohensHome <- tval *sqrt(((used_N_t + used_N_c)/(used_N_t*used_N_c)) *((used_N_t + used_N_c)/(used_N_t+used_N_c-2)))

# summary stats T2 for SI Table S3
MIpost2 %>%
  group_by(Condition, Intervention_Phase) %>%
  summarise(
    count = n(),
    mean = mean(Angst_im_moment),
    sd = sd(Angst_im_moment)
  )

# 4.2. MI PostHoc Test at T1; Intervention_Phase == post1 (Acute Effects)------------------------------------------------------------------------------------
MIpost1 <- model_MI_average[model_MI_average$Intervention_Phase == "2",]
MIpost1.model<- lme(Angst_im_moment ~ factor(Condition) + Age + factor(Sex) + Angst_im_moment_baseline, random=~1|VP_Nr, data=MIpost1, na.action= "na.omit", control=lmeControl(opt="optim"))

# get the p value from posthoc at T1:
res_MI_post1_summary <- summary(MIpost1.model)    # p = 0.2430

# summary stats for SI Table S3
MIpost1 %>%
  group_by(Condition, Intervention_Phase) %>%
  summarise(
    count = n(),
    mean = mean(Angst_im_moment),
    sd = sd(Angst_im_moment)
  )

# 4.3 MI PostHoc Test at Baseline; Intervention_Phase == pre (1) ---------------------------------------------------------------------------------------------------------------------------------------------------------
MIdatpre <- model_MI_average[model_MI_average$Intervention_Phase == "1",]
MIpre.model<- lme(Angst_im_moment ~ factor(Condition) + Age + factor(Sex), random=~1|VP_Nr, data=MIdatpre, na.action= "na.omit", control=lmeControl(opt="optim"))

# get the p value from posthoc at T0:
res_MI_pre_summary <- summary(MIpre.model)    # p = 0.80 => must discuss, confused for Beni's LOCF which should actually be similar?

res_MI_pre_Anova <- Anova(MIpre.model, type = "II")                   
MI_res_post2.model_Anova <- as.data.frame(MI_res_post2.model_Anova) 

# summary stats Baseline for SI Table S3 
MIdatpre %>%
  group_by(Condition, Intervention_Phase) %>%
  summarise(
    count = n(),
    mean = mean(Angst_im_moment),
    sd = sd(Angst_im_moment)
  )


#---- 5.0 Missing complete at random #############################################################################################################################################################################################

# 5.1 transform dataset to wide
library(tidyr)
data_wide <- spread(redDat, Intervention_Phase, Angst_im_moment)
data_wide <- data_wide %>% drop_na(post1)
data_wide <- data_wide %>% drop_na(post2)

data_long <- gather(data_wide, Intervention_Phase, Angst_im_moment, pre:post2, factor_key=TRUE)

# 5.2 calculate model with all missings excluded # complete cases analysis, under the missing completely at random assumption
lme_CC<- lme(Angst_im_moment ~  Age + factor(Sex) + factor(Condition) * factor(Intervention_Phase), random=~1|VP_Nr, data=data_long,control=lmeControl(opt="optim"))

resCC <- anova(lme_CC, type = "marginal")
resCC <- as.data.frame(Anova(lme_CC, type = "II"))                   
sum_CC <- summary(lme_CC)
sum_CC <- summary(lme_CC)[[20]]
res_sum_CC <- sum_CC[2,4]

###-- PostHoc Test at T2; Intervention_Phase == post2 (HomeTreatment)---------------------------------------------------------------------------------------------------------------------------------------------------------
datpost2 <- data_long[data_long$Intervention_Phase == "post2",]
post2.model<- lme(Angst_im_moment ~ factor(Condition) + Age + factor(Sex) + Angst_im_moment_baseline, random=~1|VP_Nr, data=datpost2, na.action= "na.omit", control=lmeControl(opt="optim"))

res_post2.model_Anova <- Anova(post2.model, type = "II")                   
res_post2.model_Anova <- as.data.frame(res_post2.model_Anova) 
res_post2.model_summary <- summary(post2.model)
res_post2.model_summary <- res_post2.model_summary[[20]]
res_tValuePost2 <- res_post2.model_summary[2,4]
res_fValuePost2 <- anova(post2.model, type = "marginal")

# #---- get adjusted means for group difference Home 
emm <- emmeans(post2.model, specs = revpairwise ~ as.factor(Condition),adjust = "none",type = "response") #adjust for changing type of multiple comparison adjustment
emmgroup <- as.data.frame(emm)
emm2 <- emm$contrasts %>% summary(infer = TRUE) 

# Cohen's D HT
used_N <- length(datpost2$VP_Nr)
used_N_c <- length(which(datpost2$Condition == 0))
used_N_t <- length(which(datpost2$Condition == 1))
tval <- res_post2.model_summary["factor(Condition)1", "t-value"]

cohensHome <- res_tValuePost2 *sqrt(((used_N_t + used_N_c)/(used_N_t*used_N_c)) *((used_N_t + used_N_c)/(used_N_t+used_N_c-2)))

# 2.0 summary stats Post T2 for SI table S3
datpost2 %>%
  group_by(Condition, Intervention_Phase) %>%
  summarise(
    count = n(),
    mean = mean(Angst_im_moment),
    sd = sd(Angst_im_moment)
  )

# PostHoc Test at T1; Intervention_Phase == post1 (HomeTreatment)---------------------------------------------------------------------------------------------------------------------------------------------------------

datpost1 <- data_long[data_long$Intervention_Phase == "post1",]
post1.model<- lme(Angst_im_moment ~ factor(Condition) + Age + factor(Sex) + Angst_im_moment_baseline, random=~1|VP_Nr, data=datpost1, na.action= "na.omit", control=lmeControl(opt="optim"))

res_post1.model_Anova <- Anova(post1.model, type = "II")                   
res_post1.model_Anova <- as.data.frame(res_post1.model_Anova) 
res_post1.model_summary <- summary(post1.model)
res_post1.model_summary <- res_post1.model_summary[[20]]  # 4.514596e-01
res_tValuePost1 <- res_post1.model_summary[2,4]
res_fValuePost1 <- anova(post1.model, type = "marginal")

# N of participants included in the ITT: hereby checking if its the same as should be due to CC
used_N <- length(datpost1$VP_Nr)
used_N_c <- length(which(datpost1$Condition == 0))
used_N_t <- length(which(datpost1$Condition == 1))

# 2.0 summary stats Post Hoc Post1 

datpost1 %>%
  group_by(Condition, Intervention_Phase) %>%
  summarise(
    count = n(),
    mean = mean(Angst_im_moment),
    sd = sd(Angst_im_moment)
  )

# PostHoc Test at Baseline; Intervention_Phase == pre ------------------------------------------------------------------------------------------------------

datpre <- data_long[data_long$Intervention_Phase == "pre",]
pre.model<- lme(Angst_im_moment ~ factor(Condition) + Age + factor(Sex), random=~1|VP_Nr, data=datpre, na.action= "na.omit", control=lmeControl(opt="optim"))

res_pre.model_Anova <- Anova(pre.model, type = "II")                   
res_pre.model_Anova <- as.data.frame(res_pre.model_Anova) 
res_pre.model_summary <- summary(pre.model)
res_pre.model_summary <- res_pre.model_summary[[20]]  # 0.985
res_tValuePre <- res_pre.model_summary[2,4]
res_fValuePre <- anova(pre.model, type = "marginal")

# N of participants included in the ITT: hereby checking if its the same as should be due to CC
used_N <- length(datpre$VP_Nr)
used_N_c <- length(which(datpre$Condition == 0))
used_N_t <- length(which(datpre$Condition == 1))

# 2.0 summary stats T1 for SI Table S3 

datpre %>%
  group_by(Condition, Intervention_Phase) %>%
  summarise(
    count = n(),
    mean = mean(Angst_im_moment),
    sd = sd(Angst_im_moment)
  )

