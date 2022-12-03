#### Main script to infer treatment effects of ASSIST study ################################################################
# chapters: 
# 1.0 get the data
# 2.0 summary stats
# 3.0 calculate the models
# 4.0 plot the data
# 5.0 Posthoc tests per condition
# 6.0 Influence of how the training was done on outomce (treatment group only)
# Beni - 2019.08.15
#### 0.0 start header ######################################################################################################

#---- set system language to "en", clean workspace, scientific notation ----------------------------------------------------
rm(list = ls())
options("scipen"=100, "digits"=5)
options(stringsAsFactors = F)
#options(contrasts=c(unordered="contr.sum", ordered="contr.poly"))

#---- load all the required packages ---------------------------------------------------------------------------------------
library(gtools)
library(ggplot2)
library(reshape)
library(gridExtra) # for saving two ggplots together
library(tidyr) # short-to-long format manipulations
library(shape)
library(grid)
library(gridGraphics)
library(d3heatmap)
library(lattice)
library(openxlsx) # to write excel-files
library(lsr) # to calculate effect sizes out of t-tests
library(r2glmm) #  to obtain estimates for specific fixed effects
library(ggpubr) # to arrange multiple plots with one common legend
library(Rmisc)

# ---- Nathalie's packages -------------------------------------------------------------------------------------------------
library(nlme)
library(car)
library(plyr)
library(emmeans) #for adjusted means
library(magrittr) #for confidence intervals of adjusted means group comparision


#---- set directories ------------------------------------------------------------------------------------------------------
# Choose where to write results to + give the folder a meaningful name
PROCESSING_PATH = '/Volumes/dcn_assist-8/PROCESSING/'
odir =  paste0(PROCESSING_PATH, 'output/', gsub("-", "_", Sys.Date()), "_", "influence_of_S8logfileParameters_on_outcome/")# output-directory
dir.create(odir)
sdir = paste0(PROCESSING_PATH, 'scripts/')

# ATTENTION WHERE TO GET DATA FROM!!!
idir_LABKEY = "/Volumes/dcn_assist-8/Daten/Labkey/"
idir_MAIN = "/Volumes/dcn_assist-8/PREPROCESSING/output/main/2020_03_03_primary_variables_only/"  # data from preprocessing main-pipeline
idir_RATINGS = "/Volumes/dcn_assist-8/PREPROCESSING/output/ratings/2020_03_05_ratings/"  # data from preprocessing ratings-pipeline
idir_QUESTIONNAIRES_CORT = "/Volumes/dcn_assist-8/PREPROCESSING/output/questionnaires/2020_03_06_includeIPQ_USX_socImm_selfExp/" # data from preprocessing ratings-pipeline
idir_LOGFILES = "/Volumes/dcn_assist-8/PREPROCESSING/output/logfiles/2020_03_05_logfiles/" # data from s8-logfiles (treatment group only)
idir_EYETRACKING = "/Volumes/dcn_assist-8/PREPROCESSING/output/eyetracking/2020_03_12_all_VP_with_AOI/" 

alldirs = c(odir, sdir, idir_LABKEY, idir_MAIN, idir_RATINGS, idir_QUESTIONNAIRES_CORT, idir_LOGFILES, idir_EYETRACKING)
setwd(odir)

# ---- source the scripts --------------------------------------------------------------------------------------------------
source(paste0(sdir, "Rallfun-v31.txt"))# robust-statistic-functions of Rand Wilcox
source (paste0(sdir, "rcontrast.R"))
source(paste0(sdir, "Plot_BF.R"))

#---- define the variables -------------------------------------------------------------------------------------------------
# 1.0
names_MAIN <- c("Angst_im_moment", "Angst_Augenkontakt") 
names_RATINGS_committee <- c("Time_s.committee", "Redefluss.committee", "Ausdruck.committee", "Sprechweise.committee","Sprechtempo.committee",      
                             "Gestik_Haltung.committee", "Mimik.committee","Blickkontakt.committee","Nervositaet.committee","Gesamturteil.committee")
names_RATINGS_VP <- c("Redefluss.VP", "Ausdruck.VP", "Sprechweise.VP","Sprechtempo.VP",      
                      "Gestik_Haltung.VP", "Mimik.VP","Blickkontakt.VP","Nervositaet.VP","Gesamturteil.VP")
names_QUESTIONNAIRES <- c("FNEK_Score","SPIN_Score","IPQ_Score","USX_Score","selfExp_totalDuration","socImm_Score", "Angst_Verbesserung","Blickv_Verbesserung","Improvisation_Verbesserung")
names_CORT <- c("deltaCort")
names_EYETRACKING <- c("pupil_mean","pupil_median", "pupil_sd", 
                       "dwelltime","dwelltimeAOI","dwelltimeNonAOI","relativeDwelltimeAOI",
                       "nfix","nfixAOI","nfixNonAOI","relativeNfixAOI",
                       "AOI_visits", "AOIdurationInequality", "speechtime",
                       "blinktime","relativeBlinktime","nBlinks","relativeNBlinks")

names_DV <- c(names_MAIN,names_RATINGS_committee,names_RATINGS_VP, names_QUESTIONNAIRES, names_CORT,names_EYETRACKING) 

completelyExcludedSubs <- "102"   # subjects that are completely excluded, meaning their data should also not be filled in by BOCF or LOCF (exclusion due to protocoll violation)
# choose type of analysis
# "PROTOCOLL": NAs left as they are
# "BOCF": baseline observation carried forward
# "LOCF": last observation carried forward
# for all analysis options: If a participant has values that were generated after dropout (f.e. if training was later identified as not done correctly),
# these values are carried forward instead of baseline/last observation
# IMPORTANT: in order to decide if BOCF or LOCF should be preferred, Nathalies idea was to look at the protocoll analysis. 
# This analysis reveals that the control group seems to improve from "pre" to "post1" and at post2 fall back to "pre" again. 
# This makes the LOCF more plausible, because by carrying forward the baseline, one would discard the actual improvement 
# of the control group, which might artificially smaller the effects --> LOCF is preferred, but BOFC would be more conservative.
# LOCF resolves an additional issue: Values only measured at "post1" and "post2" ("Verbesserung..."-measures) can't be filled in by a baseline, but by LOCF
# In these cases, BOCF is similar to PROTOCOLL, but not identical since f.e. treatment dropouts of visit 2 keep their "post2"-values
# In these cases, BOCF is similar to LOCF, but not identical since f.e. NA values at "post2" are filled in in the latter case
analysis <- "PROTOCOLL"       # alternatives: c("PROTOCOLL", "BOCF", "LOCF")

# 2.0
plotStyle <- "PointLinePlot"   # alternatives: c("PointLinePlot", "Boxplot")

#---- create log file ------------------------------------------------------------------------------------------------------
sink('PROCESSING_MAIN_ASSIST_logfile.txt', append = F)
print(paste0("Date: ", Sys.Date()))
print(paste0("input_directory LABKEY: ", idir_LABKEY))
print(paste0("input_directory MAIN: ", idir_MAIN))
print(paste0("input_directory RATINGS: ", idir_RATINGS))
print(paste0("input_directory QUETIONNAIRES and CORTISOL: ", idir_QUESTIONNAIRES_CORT))
print(paste0("output_directory: ", odir))
print(paste0("scripts_directory: ", sdir))
print(paste0("names DV: ", paste(names_DV, collapse= ", ")))
print(paste0("completelyExcludedSubs: ", completelyExcludedSubs))
print(paste0("analysis: ", analysis))
print(paste0("plotStyle: ", plotStyle))
sink()
#### end header ############################################################################################################


#### 1.0 get the data ######################################################################################################
sml = read.csv(file  = paste0(idir_LABKEY, "99_SML_v1.csv"),  sep = ',', header = T) 
load(paste0(idir_MAIN, "dat_aggreg.RData"))
dat_MAIN <- dat_aggreg
dat_EMPTY <- data.frame("VP_Nr" =  unique(dat_MAIN$VP_Nr), "pre" = "pre", "post1" = "post1","post2" = "post2")
dat_EMPTY <- merge(dat_EMPTY, sml, by = "VP_Nr")
dat_EMPTY <- reshape(dat_EMPTY,  idvar = c("VP_Nr"),varying = list(2:4),v.names = "Intervention_Phase", direction = "long")
dat_EMPTY <- merge(dat_EMPTY, unique(dat_MAIN[,c("VP_Nr", "Age", "Sex")]), by = c("VP_Nr"))  
load(paste0(idir_RATINGS, "dat_aggreg_wide.RData"))
dat_RATINGS <- dat_aggreg_wide
load(paste0(idir_QUESTIONNAIRES_CORT, "dat_aggreg.RData"))
dat_QUESTIONNAIRES_CORT <- dat_aggreg[,-c(which(names(dat_aggreg) %in% c("Angst_im_moment","Angst_Augenkontakt")))]
load(paste0(idir_LOGFILES, "dat_aggreg.RData"))
dat_LOGFILES <- dat_aggreg
load(paste0(idir_EYETRACKING, "dat_aggreg.RData"))
dat_EYETRACKING <- dat_aggreg
dat_Analyse2 <- merge(dat_EMPTY, dat_MAIN[,-c(which(names(dat_MAIN) %in% c("Age", "Sex", "Condition")))], by = c("VP_Nr", "Intervention_Phase"), all = T)
dat_Analyse1.5 <- merge(dat_Analyse2, dat_LOGFILES,  by = c("VP_Nr", "Intervention_Phase"),all = T)
dat_Analyse1.4 <- merge(dat_Analyse1.5, dat_EYETRACKING,  by = c("VP_Nr", "Intervention_Phase"),all = T)
datAnalyse1 <- merge(dat_Analyse1.4, dat_RATINGS, by = c("VP_Nr", "Intervention_Phase", "Age", "Sex", "Condition"), all = T)
datAnalyse <- merge(datAnalyse1, dat_QUESTIONNAIRES_CORT, by = c("VP_Nr", "Intervention_Phase", "Age", "Sex", "Condition"), all = T)
datAnalyse <- datAnalyse[!datAnalyse$VP_Nr %in% completelyExcludedSubs,]
datAnalyse <- datAnalyse[with(datAnalyse, order(VP_Nr, time)), ]
datAnalyse$Intervention_Phase <- factor(datAnalyse$Intervention_Phase, levels = c("pre", "post1", "post2"))
datAnalyse[datAnalyse$Condition == 2, "Condition"] <- 0 # change Controlgroup to 0 in order to make t-values more easy to interpret
datAnalyse[datAnalyse$Dropout_visit_1 == 1 & datAnalyse$Intervention_Phase %in% c("post1","post2"), names_DV] <- NA  # remove values of "post1" and "post2" for dropouts of visit 1
datAnalyse[datAnalyse$Dropout_visit_2 == 1 & datAnalyse$Intervention_Phase == "post2", names_DV] <- NA  # remove values of "post2" for dropouts of visit 2

# Intention to treat options
if (analysis == "PROTOCOLL"){
  #datAnalyse <- datAnalyse[datAnalyse$Dropout_visit_1 != 1,]  # remove dropouts of visit 1 completely (????)
} else if (analysis == "BOCF"){
  for (n in 1:nrow(datAnalyse)){
    #n = 9
    toFill <- which(is.na(datAnalyse[n,]))
    if (datAnalyse[n,"Intervention_Phase" ] %in% c("post1", "post2")){
      datAnalyse[n,toFill] <- datAnalyse[which(datAnalyse$VP_Nr == datAnalyse[n, "VP_Nr"] & datAnalyse$Intervention_Phase == "pre"), toFill]
    }# end of "post"-checks
  } # end of n-loop
} else if (analysis == "LOCF"){
  for (n in 1:nrow(datAnalyse)){
    #n = 5
    toFill <- which(is.na(datAnalyse[n,]))
    if (datAnalyse[n,"Intervention_Phase" ]=="post1"){
      datAnalyse[n,toFill] <- datAnalyse[which(datAnalyse$VP_Nr == datAnalyse[n, "VP_Nr"] & datAnalyse$Intervention_Phase == "pre"), toFill]
    } else if  (datAnalyse[n,"Intervention_Phase" ]=="post2"){
      datAnalyse[n,toFill] <- datAnalyse[which(datAnalyse$VP_Nr == datAnalyse[n, "VP_Nr"] & datAnalyse$Intervention_Phase == "post1"), toFill]
    }# end of post2-check
  } # end of n-loop
} # end of analysis-check
# ATTENTION: after carrying forward previously observed values, the ones specified below should be put back to NA
# these represent cases where the values aren't missing, but have not been measured at all (f.e SPIN scores after first PST)
datAnalyse[which(datAnalyse$Intervention_Phase == "post1"), c("SPIN_Score","FNEK_Score","BDI_Score","deltaCort")] <- NA
datAnalyse[which(datAnalyse$Intervention_Phase == "post2" & datAnalyse$Condition==0), c("USX_Score")] <- NA # no post2 assessment of usability for control group
datAnalyse[which(datAnalyse$Intervention_Phase == "post2"), c("IPQ_Score")]  <- NA # no post2 assessment of IPQ
save (datAnalyse, file = paste0(odir,"dat_Analyse_",analysis, ".RData"))

#### 2.0 summary stats ######################################################################################################
c1=c("Condition", "Intervention_Phase")
v1=c(names_DV)
datAnalyse_mean <- data.frame(t(aggregate(datAnalyse[v1],by=datAnalyse[c1],FUN=mean, na.rm=TRUE, na.action=na.omit)))
datAnalyse_sd <- as.data.frame(t(aggregate(datAnalyse[v1],by=datAnalyse[c1],FUN=sd, na.rm=TRUE)))
descriptives_DF  <- cbind(datAnalyse_mean, datAnalyse_sd)
names(descriptives_DF) <- c(rep("mean",ncol(descriptives_DF)/2),rep("sd",ncol(descriptives_DF)/2))
save(descriptives_DF, file = paste0(odir, "descriptives_DF.RData"))
write.xlsx(descriptives_DF, file = paste0(odir, "descriptives_DF.xlsx"), row.names = T)


#### 3.0 Influence of how the training was done on outomce (treatment group only) ###########################################
names(datAnalyse)
names_IV = c("App_ASSIST_level-mean", "App_ASSIST_training_min-mean","App_ASSIST_fear-mean","App_ASSIST_timeOffTarget-mean","App_ASSIST-nTrainings") # ATTENTION: change this to f.e. time in App or amout of Trainings as soon as available!!!!
names_DV = names_DV[!names_DV %in% names_IV]
res_all_DV <- c()
res2_all_DV <- c()
res_sub_all_DV <- c()
main_effects_all <- c()
sub_effects_all <- c()

for (j in names_DV){
  #j =  "FNEK_Score" 
  name_DV = j
  print(paste0("DV = ", name_DV))
  
  for(i in names_IV){
    #i = "App_ASSIST-nTrainings"
    name_IV = i
    if (name_DV != "IPQ_Score" | name_IV != "App_ASSIST-nTrainings"){  #for this combination, analyses does not make sense because variance for in ntrainings at post 1 is 0
      
      print(paste0("IV = ", name_IV))
      datAnalyse_j <- datAnalyse[datAnalyse$Condition==1, c("VP_Nr", "Age", "Sex", "Condition", "Intervention_Phase", name_DV, name_IV)]
      
      if (length(na.omit(datAnalyse_j[datAnalyse_j$Intervention_Phase=="pre", name_DV]))>0){
        datAnalyse_baseline <- datAnalyse_j[which(!is.na(datAnalyse_j[name_DV]) & datAnalyse_j$Intervention_Phase == "pre"),]
        datAnalyse_j <- datAnalyse_j[which(!is.na(datAnalyse_j[name_DV]) & datAnalyse_j$Intervention_Phase != "pre"),]
        datAnalyse_j$DV_baseline <- apply(datAnalyse_j [,],1, function(x) datAnalyse_baseline[datAnalyse_baseline$VP_Nr== x["VP_Nr"], name_DV] )
        DV = datAnalyse_j[,name_DV]
        IV = datAnalyse_j[,name_IV]
        
        # ---- calculate ----
        # is there still data for more than one intervention phase? If not, build model without factor(Intervention_Phase)
        if (length(unique(datAnalyse_j$Intervention_Phase))!=1){
          if (i == "App_ASSIST_fear-mean" & j == "Nervositaet.VP"| i == "App_ASSIST-nTrainings"){ # special case with following error: "Singularity in backsolve at level 0, block 1" --> calculate model without interaction term (???)
            lme_mod <- lme(DV ~  Age + factor(Sex) + DV_baseline + IV + factor(Intervention_Phase) , random=~1|VP_Nr, data=datAnalyse_j,control=lmeControl(opt="optim"), na.action="na.omit")
          } else {
            lme_mod <- lme(DV ~  Age + factor(Sex) + DV_baseline + IV * factor(Intervention_Phase) , random=~1|VP_Nr, data=datAnalyse_j,control=lmeControl(opt="optim"), na.action="na.omit")
          }
        } else if (length(unique(datAnalyse_j$Intervention_Phase))==1){
          lme_mod <- lme(DV ~  Age + factor(Sex) + DV_baseline + IV, random=~1|VP_Nr, data=datAnalyse_j,control=lmeControl(opt="optim"), na.action="na.omit")
        }
      } else if (length(na.omit(datAnalyse_j[datAnalyse_j$Intervention_Phase=="pre", name_DV])) == 0){
        datAnalyse_j <- datAnalyse_j[which(!is.na(datAnalyse_j[name_DV]) & datAnalyse_j$Intervention_Phase != "pre"),]
        DV = datAnalyse_j[,name_DV]
        IV = datAnalyse_j[,name_IV]
        # ---- calculate ----
        # is there still data for more than one intervention phase? If not, build model without factor(Intervention_Phase)
        if (length(unique(datAnalyse_j$Intervention_Phase))!=1){
          if (i == "App_ASSIST_fear-mean" & j == "Nervositaet.VP" | i == "App_ASSIST-nTrainings"){ # special case with following error: "Singularity in backsolve at level 0, block 1" --> calculate model without interaction term (???)
            lme_mod <- lme(DV ~  Age + factor(Sex) + IV + factor(Intervention_Phase) , random=~1|VP_Nr, data=datAnalyse_j,control=lmeControl(opt="optim"), na.action="na.omit")
          } else{
            lme_mod <- lme(DV ~  Age + factor(Sex) + IV * factor(Intervention_Phase), random=~1|VP_Nr, data=datAnalyse_j,control=lmeControl(opt="optim"), na.action="na.omit")
          }
        }else if (length(unique(datAnalyse_j$Intervention_Phase))==1){
          lme_mod <- lme(DV ~  Age + factor(Sex) +IV, random=~1|VP_Nr, data=datAnalyse_j,control=lmeControl(opt="optim"), na.action="na.omit")
        }
      } # end of check if there is a baseline to include into the model
      res2 <- summary(lme_mod)
      res2_anova <- anova(lme_mod, type = "marginal")
      res2 <- as.data.frame(res2[[20]])
      res2$DV = name_DV
      res2$IVs = row.names(res2)
      res2$r <- rcontrast(res2[,"t-value"], res2[,"DF"])
      row.names(res2) <- NULL
      res2 <- res2[,c("DV","IVs", "Value","Std.Error", "DF" ,"t-value","p-value", "r")]
      res2_dfl <- res2
      
      # main analyses
      res <- Anova(lme_mod, type = "II") # specified model
      res$DV = name_DV
      res$IVs = row.names(res)
      row.names(res) <- NULL
      
      # write columns from res2 to res 
      for (i in c("Value", "Std.Error", "t-value", "r")){ 
        #print(i)
        res[,i] <- NA
        res[res$IVs == "Age",i] <- res2[res2$IVs == "Age", i]
        res[res$IVs == "factor(Sex)",i] <- res2[res2$IVs == "factor(Sex)2", i]
        res[res$IVs == "IV",i] <- res2[res2$IVs == "IV", i]
        res[res$IVs == "factor(Intervention_Phase)",i] <- res2[res2$IVs == "factor(Intervention_Phase)post2", i]
      }
      
      res <- res[,c("DV","IVs", "Chisq", "Df", "Pr(>Chisq)", "Value", "t-value", "r")]
      res_dfl <- cbind(res, res2_anova[-1,])
      res <- res_dfl
      res[res$IVs == "IV", "IVs"] <- paste0("IV = ", name_IV)
      interaction_p = res[res$IVs == "IV:factor(Intervention_Phase)", "Pr(>Chisq)"]
      posthoc <-  length(interaction_p) != 0 
      F_value = NA
      if (posthoc == T){
        posthoc <- interaction_p < 0.05 
        F_value = res[res$IVs == "IV:factor(Intervention_Phase)", "F-value"]
      }
      
      # side analyses
      if (posthoc == F) {
        print(paste0("no interactions (p = ", interaction_p, "; F = ", F_value," : calculate main effects" ))
        if ("DV_baseline" %in% names(datAnalyse_j)){
          if (length(unique(datAnalyse_j$Intervention_Phase))!=1){
            lme_mod2 <- lme(DV ~  Age + factor(Sex) + DV_baseline + IV + factor(Intervention_Phase) , random=~1|VP_Nr, data=datAnalyse_j,control=lmeControl(opt="optim"), na.action="na.omit")
          } else if (length(unique(datAnalyse_j$Intervention_Phase))==1) {
            lme_mod2 <- lme(DV ~  Age + factor(Sex) + DV_baseline + IV, random=~1|VP_Nr, data=datAnalyse_j,control=lmeControl(opt="optim"), na.action="na.omit")
          }
        } else {
          if (length(unique(datAnalyse_j$Intervention_Phase))!=1){
            lme_mod2 <- lme(DV ~  Age + factor(Sex) + IV + factor(Intervention_Phase) , random=~1|VP_Nr, data=datAnalyse_j,control=lmeControl(opt="optim"), na.action="na.omit")
          } else if (length(unique(datAnalyse_j$Intervention_Phase))==1) {
            lme_mod2 <- lme(DV ~  Age + factor(Sex) + IV, random=~1|VP_Nr, data=datAnalyse_j,control=lmeControl(opt="optim"), na.action="na.omit")
          }
        } # end of check if there is a baseline to include into the model
        
        res2 <- summary(lme_mod2)
        res2_anova <- anova(lme_mod2, type = "marginal")
        res2 <- as.data.frame(res2[[20]])
        res2$DV = name_DV
        res2$IVs = row.names(res2)
        res2$r <- rcontrast(res2[,"t-value"], res2[,"DF"])
        row.names(res2) <- NULL
        res2 <- res2[,c("DV","IVs", "Value","Std.Error", "DF" ,"t-value","p-value", "r")]
        
        # main analyses
        res <- Anova(lme_mod2, type = "II") # specified model
        res$DV = name_DV
        res$IVs = row.names(res)
        row.names(res) <- NULL
        
        #----get effect sizes ----
        res3 <- lme_mod2
        res3_sgv <- r2beta(res3, method = 'sgv')
        effect_IV <- res3_sgv[which(res3_sgv$Effect == "IV"),]
        effect_DV <- cbind(name_DV,effect_IV)
        main_effects_all <- rbind(main_effects_all, effect_DV)
        
        # write columns from res2 to res (for ET-variable) 
        for (i in c("Value", "Std.Error", "t-value", "r")){ 
          #print(i)
          res[,i] <- NA
          res[res$IVs == "Age",i] <- res2[res2$IVs == "Age", i]
          res[res$IVs == "factor(Sex)",i] <- res2[res2$IVs == "factor(Sex)2", i]
          res[res$IVs == "IV",i] <- res2[res2$IVs == "IV", i]
          res[res$IVs == "factor(Intervention_Phase)",i] <- res2[res2$IVs == "factor(Intervention_Phase)post2", i]
        }
        
        res <- res[,c("DV","IVs", "Chisq", "Df", "Pr(>Chisq)", "Value", "t-value", "r")]
        res <- cbind(res, res2_anova[-1,])
        which_r_to_add <- c(which(!res_dfl$IVs %in% res$IVs))
        res<- rbind(res,res_dfl[which_r_to_add,])
        which_r2_to_add <- c(which(!res2_dfl$IVs %in% res2$IVs))
        res2 <- rbind(res2,res2_dfl[which_r2_to_add,])
        res[res$IVs == "IV", "IVs"] <- paste0("IV = ", name_IV)
        
        
      } else if (posthoc == T) {# end no intercations-loop 
        print(paste0("post-hoc tests for Condition-Invervention_Phase-interaction"))
        res_sub_all_phases <- c()
        for (y in c(unique(as.character(datAnalyse_j$Intervention_Phase)))){
          #print(y)
          #y = "post2"
          res_sub <- c()
          dat_sub = datAnalyse_j[datAnalyse_j$Intervention_Phase %in% y,]
          DV_sub = dat_sub [,name_DV]
          IV_sub = dat_sub [,name_IV]
          
          # in case of existant pre-values: include pre-values as a covariate to correct for baseline differences
          if ("DV_baseline" %in% names(datAnalyse_j)){
            lme_mod_sub <- lme(DV_sub ~  Age + factor(Sex) + DV_baseline + IV_sub , random=~1|VP_Nr, data=dat_sub,control=lmeControl(opt="optim"), na.action="na.omit")
          } else {
            lme_mod_sub <- lme(DV_sub ~  Age + factor(Sex) + IV_sub, random=~1|VP_Nr, data=dat_sub,control=lmeControl(opt="optim"), na.action="na.omit")
          } # end of check if there is a baseline to include into the model
          
          res2_sub <- summary(lme_mod_sub)
          res2_sub <- as.data.frame(res2_sub[[20]])
          res2_sub$DV = name_DV
          res2_sub$IVs = row.names(res2_sub)
          res2_sub$r <- rcontrast(res2_sub[,"t-value"], res2_sub[,"DF"])
          row.names(res2_sub) <- NULL
          res2_sub <- res2_sub[,c("DV","IVs", "Value","Std.Error", "DF" ,"t-value","p-value", "r")]
          
          # main analyses
          res_sub <- Anova(lme_mod_sub, type = "II") # specified model
          get_df <- anova(lme_mod_sub, type = "marginal") # specified model
          res_sub$DV = name_DV
          res_sub$IVs = row.names(res_sub)
          row.names(res_sub) <- NULL
          
          #---- get effect sizes ----
          res4_sgv <- r2beta(lme_mod_sub, method = 'sgv')
          effect_IV_sub <- res4_sgv[which(res4_sgv$Effect == "IV_sub"),]
          effect_DV_sub <- cbind(name_DV,y, effect_IV_sub) 
          sub_effects_all <- rbind(sub_effects_all, effect_DV_sub)
          
          # write columns from res2 to res (for ET-variable) 
          for (i in c("Value", "Std.Error", "t-value", "r")){ 
            #print(i)
            res_sub[,i] <- NA
            res_sub[res_sub$IVs == "Age",i] <- res2_sub[res2_sub$IVs == "Age", i]
            res_sub[res_sub$IVs == "factor(Sex)",i] <- res2_sub[res2_sub$IVs == "factor(Sex)2", i]
            if (length(grep("DV_baseline",res2_sub$IVs))>0){
              res_sub[res_sub$IVs == "DV_baseline",i] <- res2_sub[res2_sub$IVs == "DV_baseline", i]
            }
            res_sub[res_sub$IVs == "IV_sub",i] <- res2_sub[res2_sub$IVs == "IV_sub", i]
          }
          res_sub <- res_sub[,c("DV","IVs", "Chisq", "Df", "Pr(>Chisq)", "Value", "t-value", "r")]
          res_sub[res_sub$IVs == "IV_sub","IVs"] <- paste0("IV_sub = ", name_IV)
          res_sub$Intervention_Phase = y
          res_sub$Df <- get_df$denDF[1]
          res_sub_all_phases <- rbind (res_sub_all_phases,res_sub)
        }# end Intervention_Phase-loop
        res_sub_all_DV <- rbind(res_sub_all_DV, res_sub_all_phases)
      } # end significant interactions-loop
      
      res_all_DV <- rbind(res_all_DV, res)
      res2_all_DV <- rbind(res2_all_DV, res2)
    }   # end check if (name_DV != "IPQ_Score" | name_IV != "App_ASSIST-nTrainings")
  } # end of IV-loop
} # end DV-loop
write.xlsx(x = res_all_DV, file = paste0(odir,"MAIN_all_DV_trainingInfluence_",analysis, ".xlsx"), row.names = FALSE)
save(res_all_DV, file = paste0(odir,"MAIN_all_DV_trainingInfluence_",analysis, ".RData"))
write.xlsx(x = res2_all_DV, file = paste0(odir,"TVAL_all_DV_trainingInfluence_",analysis, ".xlsx"), row.names = FALSE)
save(res2_all_DV, file = paste0(odir,"TVAL_all_DV_trainingInfluence_",analysis, ".RData"))
write.xlsx(x = res_sub_all_DV, file = paste0(odir,"POSTHOC_all_DV_trainingInfluence_",analysis, ".xlsx"), row.names = FALSE)
save(res_sub_all_DV, file = paste0(odir,"POSTHOC_all_DV_trainingInfluence_",analysis, ".RData"))
write.xlsx(main_effects_all, file = paste0(odir,"MAIN_effects_all_trainingInfluence_",analysis, ".xlsx"), row.names = FALSE)
save (main_effects_all, file = paste0(odir,"MAIN_effects_all_trainingInfluence_",analysis, ".RData"))
write.xlsx(sub_effects_all, file = paste0(odir,"sub_effects_all_trainingInfluence_",analysis, ".xlsx"),row.names = FALSE)
save (sub_effects_all, file = paste0(odir,"sub_effects_all_trainingInfluence_",analysis, ".RData"))




