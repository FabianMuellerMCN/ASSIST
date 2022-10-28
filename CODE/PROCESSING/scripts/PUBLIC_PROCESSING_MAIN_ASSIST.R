#### Main public script to infer treatment effects of ASSIST study ################################################################
# chapters: 
# 1.0 get the data
# 2.0 summary stats
# 3.0 calculate the models
# 4.0 plot the data
# 5.0 Posthoc tests per condition
# 6.0 Influence of how the training was done on outomce (treatment group only)
# Fabian - 2022.10.19
#### 0.0 start header #############################################################################################################

#---- set system language to "en", clean workspace, scientific notation -----------------------------------------------------------
rm(list = ls())
options("scipen"=100, "digits"=5)
options(stringsAsFactors = F)
#options(contrasts=c(unordered="contr.sum", ordered="contr.poly"))

#---- load all the required packages ----------------------------------------------------------------------------------------------
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
library(nlme)
library(car)
library(plyr)
library(emmeans) #for adjusted means
library(magrittr) #for confidence intervals of adjusted means group comparision


#---- set directories ------------------------------------------------------------------------------------------------------------
# Choose where to write results to + give the folder a meaningful name
PROCESSING_PATH = "/Users/fmueller/Documents/GitHub/ASSIST/CODE/PROCESSING/" # MAKE SURE TO ADJUST PATH!
odir =  paste0(PROCESSING_PATH, "output/", gsub("-", "_", Sys.Date()))# output-directory
dir.create(odir)
sdir = paste0(PROCESSING_PATH, "scripts/")
setwd(odir)

# ---- source the scripts --------------------------------------------------------------------------------------------------
source(paste0(sdir, "Rallfun-v31.txt"))# robust-statistic-functions of Rand Wilcox
source (paste0(sdir, "rcontrast.R"))
source(paste0(sdir, "Plot_BF.R"))
analysis <- "PROTOCOLL"       # alternatives: c("PROTOCOLL", "BOCF", "LOCF"), no included in this script


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

# 2.0
plotStyle <- "PointLinePlot"   # alternatives: c("PointLinePlot", "Boxplot")

#### end header ############################################################################################################


#### 1.0 get the data ######################################################################################################
load("/Users/fmueller/Documents/GitHub/ASSIST/DATA/dat_Analyse_PROTOCOLL.RData")  # MAKE SURE TO ADJUST PATH!

#### 2.0 summary stats ######################################################################################################
c1=c("Condition", "Intervention_Phase")
v1=c(names_DV)
datAnalyse_mean <- data.frame(t(aggregate(datAnalyse[v1],by=datAnalyse[c1],FUN=mean, na.rm=TRUE, na.action=na.omit)))
datAnalyse_sd <- as.data.frame(t(aggregate(datAnalyse[v1],by=datAnalyse[c1],FUN=sd, na.rm=TRUE)))
descriptives_DF  <- cbind(datAnalyse_mean, datAnalyse_sd)
names(descriptives_DF) <- c(rep("mean",ncol(descriptives_DF)/2),rep("sd",ncol(descriptives_DF)/2))
save(descriptives_DF, file = paste0(odir, "descriptives_DF.RData"))
write.xlsx(descriptives_DF, file = paste0(odir, "descriptives_DF.xlsx"), row.names = T)

#### 3.0 calculate models ###################################################################################################
res_all_DV <- c()
res2_all_DV <- c()
res_sub_all_DV <- c()
main_effects_all <- c()
sub_effects_all <- c()
summary_N_all <- c()

for (j in names_DV[!names_DV %in% "socImm_Score"]){   # do not run the models for social immersion Score because this is only available at one time point and only for treatment group
  #j = "Angst_im_moment"   
  name_DV = j
  print(paste0("DV = ", name_DV))
  datAnalyse_j <- datAnalyse[, c("VP_Nr", "Age", "Sex", "Condition", "Intervention_Phase", name_DV)]
  datAnalyse_j <- datAnalyse_j[which(!is.na(datAnalyse_j[name_DV])),]
  DV = datAnalyse_j[,name_DV]
  
  # ---- N subjects available for primary outcome measure ----
  n_pre1_cond1 <- length(unique(datAnalyse_j[datAnalyse_j$Condition == 1 & datAnalyse_j$Intervention_Phase == "pre" & !is.na(datAnalyse_j[, j]),"VP_Nr"]))
  n_pre1_cond0 <- length(unique(datAnalyse_j[datAnalyse_j$Condition == 0 & datAnalyse_j$Intervention_Phase == "pre" & !is.na(datAnalyse_j[, j]),"VP_Nr"]))
  n_post1_cond1 <- length(unique(datAnalyse_j[datAnalyse_j$Condition == 1 & datAnalyse_j$Intervention_Phase == "post1" & !is.na(datAnalyse_j[, j]),"VP_Nr"]))
  n_post1_cond0 <- length(unique(datAnalyse_j[datAnalyse_j$Condition == 0 & datAnalyse_j$Intervention_Phase == "post1" & !is.na(datAnalyse_j[, j]),"VP_Nr"]))
  n_post2_cond1 <- length(unique(datAnalyse_j[datAnalyse_j$Condition == 1 & datAnalyse_j$Intervention_Phase == "post2" & !is.na(datAnalyse_j[, j]),"VP_Nr"]))
  n_post2_cond0 <- length(unique(datAnalyse_j[datAnalyse_j$Condition == 0 & datAnalyse_j$Intervention_Phase == "post2" & !is.na(datAnalyse_j[, j]),"VP_Nr"]))
  summary_N <- c(j,n_pre1_cond1,n_post1_cond1,n_post2_cond1, n_pre1_cond0,n_post1_cond0,n_post2_cond0)
  summary_N_all <- rbind(summary_N_all, summary_N)
  
  if (!name_DV %in% c("IPQ_Score", "USX_Score", "selfExp_totalDuration")){    # measures only assessed once --> directly proceed to post-hoc-testing
    # ---- calculate ----
    lme_mod<- lme(DV ~  Age + factor(Sex) + factor(Condition) * factor(Intervention_Phase), random=~1|VP_Nr, data=datAnalyse_j,control=lmeControl(opt="optim"), na.action="na.omit")
    res2 <- summary(lme_mod)
    res2_anova <- anova(lme_mod, type = "marginal")
    
    
    res2 <- as.data.frame(res2[[20]])
    res2$DV = name_DV
    res2$IVs = row.names(res2)
    #res2$`p-value`<-round(res2$`p-value`, digits = 7)
    res2$r <- rcontrast(res2[,"t-value"], res2[,"DF"])
    row.names(res2) <- NULL
    res2 <- res2[,c("DV","IVs", "Value","Std.Error", "DF" ,"t-value","p-value", "r")]
    res2_dfl <- res2
    
    # main analyses
    res <- Anova(lme_mod, type = "II") # specified model
    res$DV = name_DV
    res$IVs = row.names(res)
    #res$`Pr(>Chisq)`<- round(res$`Pr(>Chisq)`, digits = 7)
    row.names(res) <- NULL
    
    # write columns from res2 to res 
    for (i in c("Value", "Std.Error", "t-value", "r")){ 
      #print(i)
      res[,i] <- NA
      res[res$IVs == "Age",i] <- res2[res2$IVs == "Age", i]
      res[res$IVs == "factor(Sex)",i] <- res2[res2$IVs == "factor(Sex)2", i]
      res[res$IVs == "factor(Condition)",i] <- res2[res2$IVs == "factor(Condition)1", i]
    }
    
    res <- res[,c("DV","IVs", "Chisq", "Df", "Pr(>Chisq)", "Value", "t-value", "r")]
    res_dfl <- cbind(res, res2_anova[-1,])
    res <- res_dfl
    # get effect sizes (vgl. Dominiques recommendation to indicate main effects, although interactions might be significant...)
    #res1000_sgv <- r2beta(lme_mod, method = 'sgv')
    
    # ---- 3.1 if no interactions: calculate model without interaction term ----
    interaction_p = res[res$IVs == "factor(Condition):factor(Intervention_Phase)", "Pr(>Chisq)"]
    posthoc <- interaction_p < 0.05
    F_value = res[res$IVs == "factor(Condition):factor(Intervention_Phase)", "F-value"]
  } else {
    posthoc = T
    res <- c()
    res2 <- c()
  } # end exclusion of variables only measured once
  
  if (posthoc == F) {
    print(paste0("no interactions (p = ", interaction_p, "; F = ", F_value," : calculate main effects" ))
    
    # side analyses
    lme_mod2 <- lme(DV ~  Age + factor(Sex) + factor(Condition) + factor(Intervention_Phase), random=~1|VP_Nr, data=datAnalyse_j ,control=lmeControl(opt="optim"), na.action="na.omit")
    res2 <- summary(lme_mod2)
    res2_anova <- anova(lme_mod2, type = "marginal")
    #res2_anova <-rbind(res2_anova, NA)
    res2 <- as.data.frame(res2[[20]])
    res2$DV = name_DV
    res2$IVs = row.names(res2)
    #res2$`p-value`<-round(res2$`p-value`, digits = 7)
    res2$r <- rcontrast(res2[,"t-value"], res2[,"DF"])
    row.names(res2) <- NULL
    res2 <- res2[,c("DV","IVs", "Value","Std.Error", "DF" ,"t-value","p-value", "r")]
    #res2  <- rbind(res2, res2_dfl [!res2_dfl$IVs %in% res2$IVs,])
    
    # main analyses
    res <- Anova(lme_mod2, type = "II") # specified model
    res$DV = name_DV
    res$IVs = row.names(res)
    #res$`Pr(>Chisq)`<-round(res$`Pr(>Chisq)`, digits = 7)
    row.names(res) <- NULL
    #res <- rbind(res, res_dfl[!res_dfl$IVs %in% res$IVs,which(names(res_dfl) %in% names(res))])
    
    #----get effect sizes ----
    res3 <- lme_mod2
    res3_sgv <- r2beta(res3, method = 'sgv')
    effect_Condition <- res3_sgv[which(res3_sgv$Effect == "factor(Condition)1"),]
    effect_DV <- cbind(name_DV,effect_Condition)
    main_effects_all <- rbind(main_effects_all, effect_DV)
    
    # write columns from res2 to res (for ET-variable) 
    for (i in c("Value", "Std.Error", "t-value", "r")){ 
      #print(i)
      res[,i] <- NA
      res[res$IVs == "Age",i] <- res2[res2$IVs == "Age", i]
      res[res$IVs == "factor(Sex)",i] <- res2[res2$IVs == "factor(Sex)2", i]
      res[res$IVs == "factor(Condition)",i] <- res2[res2$IVs == "factor(Condition)1", i]
      #res[res$IVs == "factor(Intervention_Phase)",i] <- res2[res2$IVs == "factor(Intervention_Phase)pre", i]
      res[res$IVs == "factor(Condition):factor(Intervention_Phase)",i] <- res2_dfl[res2_dfl$IVs == "factor(Condition)1:factor(Intervention_Phase)post2", i]
    }
    
    res <- res[,c("DV","IVs", "Chisq", "Df", "Pr(>Chisq)", "Value", "t-value", "r")]
    res <- cbind(res, res2_anova[-1,])
    which_r_to_add <- c(which(!res_dfl$IVs %in% res$IVs))
    res <- rbind(res,res_dfl[which_r_to_add,])
    #res[which_r_to_add,]<-NA
    which_r2_to_add <- c(which(!res2_dfl$IVs %in% res2$IVs))
    res2 <- rbind(res2,res2_dfl[which_r2_to_add,])
    #res2[which_r2_to_add,]<-NA
    
  } else if (posthoc == T) {# end no intercations-loop 
    # ---- 3.2: if significant interactions: calculate post-hoc-tests ---- 
    print(paste0("post-hoc tests for Condition-Invervention_Phase-interaction"))
    res_sub_all_phases <- c()
    
    for (y in c(unique(as.character(datAnalyse_j$Intervention_Phase)))){
      res_sub <- c()
      #y = "post2"
      # in case of post1 or post2-analysis and existant pre-values: include pre-values as a covariate to correct for baseline differences
      if (y %in% c("post1", "post2") & "pre" %in% unique(datAnalyse_j$Intervention_Phase)){
        dat_sub = datAnalyse_j[datAnalyse_j$Intervention_Phase %in% c(y, "pre"),]
        dat_sub <- reshape(dat_sub, 
                           v.names = c(name_DV),
                           idvar = c("VP_Nr"),
                           timevar = "Intervention_Phase",
                           direction = "wide")
        names(dat_sub)[names(dat_sub) %in% paste0(name_DV, ".", y)] <- name_DV
        DV_sub = dat_sub [,name_DV]
        IV_baseline = dat_sub[, paste0(name_DV, ".pre")]
        # side analyses
        lme_mod_sub <- lme(DV_sub ~  Age + factor(Sex) + factor(Condition) + IV_baseline, random=~1|VP_Nr, data=dat_sub ,control=lmeControl(opt="optim"), na.action="na.omit")
      } else if (y == "post2" & name_DV == "USX_Score") {  # in case of usability score: post2 is only assessed for treatment group and therefore not comparable between groups. In that case, do nothing
      } else { 
        dat_sub = datAnalyse_j[datAnalyse_j$Intervention_Phase %in% y,]
        DV_sub = dat_sub [,name_DV]
        # side analyses
        lme_mod_sub <- lme(DV_sub ~  Age + factor(Sex) + factor(Condition), random=~1|VP_Nr, data=dat_sub ,control=lmeControl(opt="optim"), na.action="na.omit")
      } # end of post1 or post2 or pre-check
      
      res2_sub <- summary(lme_mod_sub)
      res2_sub <- as.data.frame(res2_sub[[20]])
      res2_sub$DV = name_DV
      res2_sub$IVs = row.names(res2_sub)
      #res2_sub$`p-value`<-round(res2_sub$`p-value`, digits = 7)
      res2_sub$r <- rcontrast(res2_sub[,"t-value"], res2_sub[,"DF"])
      row.names(res2_sub) <- NULL
      res2_sub <- res2_sub[,c("DV","IVs", "Value","Std.Error", "DF" ,"t-value","p-value", "r")]
      
      
      # Cohen's D
      used_N <- length(dat_sub$VP_Nr)
      used_N_c <- length(which(dat_sub$Condition == 0))
      used_N_t <- length(which(dat_sub$Condition == 1))
      tval <- res2_sub[res2_sub$IVs == "factor(Condition)1", "t-value"]
      d <- tval*sqrt(((used_N_t + used_N_c)/(used_N_t*used_N_c)) *((used_N_t + used_N_c)/(used_N_t+used_N_c-2)))
      
      # main analyses
      res_sub <- Anova(lme_mod_sub, type = "II") # specified model
      get_df <- anova(lme_mod_sub, type = "marginal") # specified model
      res_sub$DV = name_DV
      res_sub$IVs = row.names(res_sub)
      #res_sub$`Pr(>Chisq)`<-round(res_sub$`Pr(>Chisq)`, digits = 7)
      row.names(res_sub) <- NULL
      
      #---- get adjusted means for group difference ----
      emm <- emmeans(lme_mod_sub, specs = revpairwise ~ as.factor(Condition),adjust = "none",type = "response") #adjust for changing type of multiple comparison adjustment
      emmgroup <- as.data.frame(emm)
      emm1 <- emm$contrasts %>%summary(infer = TRUE) 
      
      #---- get effect sizes ----
      res4_sgv <- r2beta(lme_mod_sub, method = 'sgv')
      effect_Condition_sub <- res4_sgv[which(res4_sgv$Effect == "factor(Condition)1"),]
      effect_DV_sub <- cbind(name_DV,y, effect_Condition_sub, emm1,d) #  emmgroup[1,], emmgroup[2,]
      sub_effects_all <- rbind(sub_effects_all, effect_DV_sub)
      
      # write columns from res2 to res (for ET-variable) 
      for (i in c("Value", "Std.Error", "t-value", "r")){ 
        #print(i)
        res_sub[,i] <- NA
        res_sub[res_sub$IVs == "Age",i] <- res2_sub[res2_sub$IVs == "Age", i]
        res_sub[res_sub$IVs == "factor(Sex)",i] <- res2_sub[res2_sub$IVs == "factor(Sex)2", i]
        if (length(grep("IV_baseline",res2_sub$IVs))>0){
          res_sub[res_sub$IVs == "IV_baseline",i] <- res2_sub[res2_sub$IVs == "IV_baseline", i]
        }
        res_sub[res_sub$IVs == "factor(Condition)",i] <- res2_sub[res2_sub$IVs == "factor(Condition)1", i]
      }
      
      res_sub <- res_sub[,c("DV","IVs", "Chisq", "Df", "Pr(>Chisq)", "Value", "t-value", "r")]
      res_sub$Intervention_Phase = y
      res_sub$Df <- get_df$denDF[1]
      res_sub[res_sub$IVs == "factor(Condition)","CohensD"] <- d 
      
      res_sub_all_phases <- rbind (res_sub_all_phases,res_sub)
    }# end Intervention_Phase-loop
    res_sub_all_DV <- rbind(res_sub_all_DV, res_sub_all_phases)
  } # end significant interactions-loop
  
  res_all_DV <- rbind(res_all_DV, res)
  res2_all_DV <- rbind(res2_all_DV, res2)
  
} # end DV-loop
summary_N_all <- as.data.frame(summary_N_all)
names(summary_N_all) <- c("DV","N_treatment_pre","N_treatment_post1","N_treatment_post2", "N_control_pre","N_control_post1", "N_control_post2")

write.xlsx(x = res_all_DV, file = paste0(odir,"MAIN_all_DV_",analysis, ".xlsx"), row.names = FALSE)
save(res_all_DV, file = paste0(odir,"MAIN_all_DV_",analysis, ".RData"))
write.xlsx(x = res2_all_DV, file = paste0(odir,"TVAL_all_DV_",analysis, ".xlsx"), row.names = FALSE)
save(res2_all_DV, file = paste0(odir,"TVAL_all_DV_",analysis, ".RData"))
write.xlsx(x = res_sub_all_DV, file = paste0(odir,"POSTHOC_all_DV_",analysis, ".xlsx"), row.names = FALSE)
save(res_sub_all_DV, file = paste0(odir,"POSTHOC_all_DV_",analysis, ".RData"))
write.xlsx(main_effects_all, file = paste0(odir,"MAIN_effects_all_",analysis, ".xlsx"), row.names = FALSE)
save (main_effects_all, file = paste0(odir,"MAIN_effects_all_",analysis, ".RData"))
write.xlsx(sub_effects_all, file = paste0(odir,"sub_effects_all_",analysis, ".xlsx"),row.names = FALSE)
save (sub_effects_all, file = paste0(odir,"sub_effects_all_",analysis, ".RData"))
write.xlsx(summary_N_all, file = paste0(odir,"summary_N_all_",analysis, ".xlsx"), row.names = FALSE)
save (summary_N_all, file = paste0(odir,"summary_N_all_",analysis, ".RData"))


#### 4.0 plot the data ######################################################################################################
# plotting subjective fear data 
if (plotStyle == "Boxplot"){
  plot_BF <- boxplot_BF
  getLegend_BF <- legendBoxplot_BF
} else if (plotStyle == "PointLinePlot"){
  plot_BF <- pointLinePlot_BF
  getLegend_BF <- legendPointLinePlot_BF
}

# plotting subjective fear (of eyecontact) data 
datPlot <- datAnalyse
whichVars_plot <- c("Angst_im_moment","Angst_Augenkontakt")  
#whichVars_plot <- c("Angst_im_moment")  # manual settings for paper
p <- lapply(c(whichVars_plot), plot_BF)
legend <- getLegend_BF (whichVars_plot[1]) 
p <- append (p,list(legend))
#p <- append(list(bf1),list(legend))  # manual settings for paper
bf2 <- grid.arrange(grobs=p, 
                    ncol=2, 
                    top=textGrob("Subj. empfundene Angst (vor Augenkontakt) in Abhängigkeit von Zeitpunkt und Gruppenzuteilung", gp=gpar(fontsize=12)))

# plotting eyetracking measures
datPlot <- datAnalyse
whichVars_plot <- names_EYETRACKING    #"Angst_im_moment"   ##"relativeDwelltimeAOI"
p <- lapply(c(whichVars_plot),plot_BF)
p <- append(p,list(legend))
bf3 <- grid.arrange(grobs=p, 
                    ncol = 2, 
                    top=textGrob("Eyetracking measures", gp=gpar(fontsize=10)))

# plotting global external and self-assessment of performance
load("/Users/fmueller/Documents/GitHub/ASSIST/DATA/dat_ratings_aggreg.RData") # ADAPT PATH


datPlot <- dat_aggreg
whichVars_plot <- gsub(".committee","", names_RATINGS_committee)
p <- lapply(c(whichVars_plot),plot_BF)
p <- append(p,list(legend))
bf4 <- grid.arrange(grobs=p, 
                    ncol = 4, 
                    top=textGrob("Bewertung der Ansprache durch das Gremium und die Probanden", gp=gpar(fontsize=10)))

# plotting FNEK and SPIN, CORT, and GlobalImprovements
datPlot <- datAnalyse
whichVars_plot <- c(names_QUESTIONNAIRES, names_CORT)
p <- lapply(c(whichVars_plot), plot_BF)
p <- append (p,list(legend))
bf5 <- grid.arrange (grobs=p, 
                     ncol = 3, 
                     top = textGrob("Fragebogen und Kortisol",gp=gpar(fontsize=10)))

pdf(file = paste0(odir, "plots_", analysis, ".pdf"))
plot(bf2)
plot(bf3)
plot(bf4)
plot(bf5)
dev.off()

#### 5.0 POST-HOC Per Condition ######################################################################################################################################
# sub_effects_all_Intervention_Phase <- c()
# res_sub_all_DV_Intervention_Phase <- c()
# res_sub_all_Conditions <- c()
# res_sub_all_phases <- c()
# 
# for (j in names_DV){
#   #j = "Angst_Augenkontakt"
#   #j = "Angst_im_moment"
#   name_DV = j
#   datAnalyse_j <- datAnalyse[, c("VP_Nr", "Age", "Sex", "Condition", "Intervention_Phase", name_DV)]
#   print(paste0("DV = ", name_DV))
#   DV = datAnalyse_j[,name_DV]
#   
#   for (y in c(0, 1)){
#     res_sub <- c()
#     #y = 1
#     print(paste0("Condition = ", y))
#     dat_sub = datAnalyse_j[datAnalyse_j$Condition == y & datAnalyse_j$Intervention_Phase %in% c("pre", "post2"),]
#     DV_sub = dat_sub [,name_DV]
#     
#     # side analyses
#     lme_mod_sub <- lme(DV_sub ~  Age + factor(Sex) + factor(Intervention_Phase), random=~1|VP_Nr, data=dat_sub ,control=lmeControl(opt="optim"), na.action="na.omit")
#     res2_sub <- summary(lme_mod_sub)
#     res2_sub <- as.data.frame(res2_sub[[20]])
#     res2_sub$DV = name_DV
#     res2_sub$IVs = row.names(res2_sub)
#     res2_sub$`p-value`<-round(res2_sub$`p-value`, digits = 7)
#     res2_sub$r <- rcontrast(res2_sub[,"t-value"], res2_sub[,"DF"])
#     row.names(res2_sub) <- NULL
#     res2_sub <- res2_sub[,c("DV","IVs", "Value","Std.Error", "DF" ,"t-value","p-value", "r")]
#     
#     # main analyses
#     res_sub <- Anova(lme_mod_sub, type = "II") # specified model
#     get_df <- anova(lme_mod_sub, type = "marginal") # specified model
#     
#     res_sub$DV = name_DV
#     res_sub$IVs = row.names(res_sub)
#     res_sub$`Pr(>Chisq)`<-round(res_sub$`Pr(>Chisq)`, digits = 7)
#     row.names(res_sub) <- NULL
#     
#     
#     # ---- get effect sizes ----
#     
#     res4_sgv <- r2beta(lme_mod_sub, method = 'sgv')
#     #effect_InterventionPhase_post1_sub <- res4_sgv[which(res4_sgv$Effect == "factor(Intervention_Phase)post1"),]
#     effect_InterventionPhase_post2_sub <- res4_sgv[which(res4_sgv$Effect == "factor(Intervention_Phase)post2"),]
#     effect_DV_sub <- cbind(name_DV,y,  effect_InterventionPhase_post2_sub) # effect_InterventionPhase_post1_sub
#     sub_effects_all_Intervention_Phase <- rbind(sub_effects_all_Intervention_Phase , effect_DV_sub)
#     
#     # write columns from res2 to res (for ET-variable)
#     for (i in c("Value", "Std.Error", "t-value", "r")){
#       #print(i)
#       res_sub[,i] <- NA
#       res_sub[res_sub$IVs == "Age",i] <- res2_sub[res2_sub$IVs == "Age", i]
#       res_sub[res_sub$IVs == "factor(Sex)",i] <- res2_sub[res2_sub$IVs == "factor(Sex)2", i]
#       res_sub[res_sub$IVs == "factor(Condition)",i] <- res2_sub[res2_sub$IVs == "factor(Condition)1", i]
#     }
#     
#     res_sub <- res_sub[,c("DV","IVs", "Chisq", "Df", "Pr(>Chisq)", "Value", "t-value", "r")]
#     res_sub$Condition = y
#     res_sub$Df <- get_df$denDF[1]
#     res_sub_all_Conditions <- rbind (res_sub_all_Conditions,res_sub)
#   }# end Condition-loop
#   
#   res_sub_all_DV_Intervention_Phase <- rbind(res_sub_all_DV_Intervention_Phase, res_sub_all_Conditions)
# }# end of j-loop


#### 6.0 Report habituation effects from Baseline to time point 1 analyses (added, 28.10.2022) ####################################################################################################
# Significant decrease in fear in both the VR gaze and the control group from baseline to post intervention phase 1 (F [1,84] = 74.53, p < 0.0001), 
# independent of group allocation with no interaction between group and time point (F [1,84] = 1.17, p = 0.28).   

# dataPrePost1 <- datAnalyse[datAnalyse$Intervention_Phase != "post2",]
# lme_modt1_2<- lme(Angst_im_moment ~  Age + factor(Sex) + factor(Condition) * factor(Intervention_Phase), random=~1|VP_Nr, data=dataPrePost1,control=lmeControl(opt="optim"), na.action="na.omit")
# res2 <- summary(lme_modt1_2)
# res2_anova <- anova(lme_modt1_2, type = "marginal")
# res_lme_modt1_2 <- Anova(lme_modt1_2, type = "II")                   
# res_lme_modt1_2 <- as.data.frame(res_lme_modt1_2) 
# res_lme_modt1_2_summary <- summary(lme_modt1_2)
# res_lme_modt1_2_summary <- res_lme_modt1_2_summary[[20]]
# res_lme_modt1_2 <- res_lme_modt1_2_summary[2,4]
# res_lme_modt1_2 <- anova(lme_modt1_2, type = "marginal")


#### 7.0 side analyses ####################################################################################################
# #is dwell time negatively related to dwell time at baseline (irrespecitve of group)?
#cor.test( ~ relativeDwelltimeAOI + Angst_im_moment, data = datAnalyse[datAnalyse$Intervention_Phase == "pre",])
# #p = 0.44 --> no association; also take into account covariates:
# lme_mod <- lme(Angst_im_moment ~  Age + factor(Sex) + factor(Condition) + relativeDwelltimeAOI, random=~1|VP_Nr, data=datAnalyse[datAnalyse$Intervention_Phase == "pre",],control=lmeControl(opt="optim"), na.action="na.omit")
# res2 <- summary(lme_mod)
# res2_anova <- Anova(lme_mod, type = "II")
# #p = 0.48 (0.49 without Condition) --> no association; formally test Intervention_Phase-Interaction as well:
# lme_mod <- lme(Angst_im_moment ~  Age + factor(Sex) + factor(Condition) * factor(Intervention_Phase) + relativeDwelltimeAOI, random=~1|VP_Nr, data=datAnalyse,control=lmeControl(opt="optim"), na.action="na.omit")
# res2 <- summary(lme_mod)
# res2_anova <- Anova(lme_mod, type = "3")
# #all p > 0.16 is is dwell time negatively related to dwell time, ONLY AT TOPIC 1 OF BASELINE (irrespecitve of group --> vgl. Mail Dominique, 22.05.2020)?
# ET_Fear_comparison <- "beforeSpeech"  # options: "afterSpeech", "beforeSpeech" should ET in speech 1 be compared to fear rating AFTER OR BEFORE speech 1
# load (paste0(idir_EYETRACKING, "dat_long.RData"))
# dat_EYETRACKING_long <- dat_long
# load (paste0(idir_MAIN, "dat_long.RData"))
# dat_MAIN_long <- dat_long
# 
# #assign fear-rating to a speech topic
# if (ET_Fear_comparison == "afterSpeech"){
#   dat_MAIN_long[dat_MAIN_long$timepoint == 2,"topic_Nr"] <- "T1"  
#   dat_MAIN_long[dat_MAIN_long$timepoint == 3,"topic_Nr"] <- "T2"   
#   dat_MAIN_long[dat_MAIN_long$timepoint == 4,"topic_Nr"] <- "T3"  
# } else if (ET_Fear_comparison == "beforeSpeech") { 
#   dat_MAIN_long[dat_MAIN_long$timepoint == 1,"topic_Nr"] <- "T1"  
#   dat_MAIN_long[dat_MAIN_long$timepoint == 2,"topic_Nr"] <- "T2"  
#   dat_MAIN_long[dat_MAIN_long$timepoint == 3,"topic_Nr"] <- "T3"   
# }# end of ET_Fear_comparison-check
# 
# datAnalyse_long <- merge(dat_MAIN_long, dat_EYETRACKING_long, by= c("VP_Nr", "Intervention_Phase", "topic_Nr"))
# cor.test( ~ relativeDwelltimeAOI + Angst_im_moment, data = datAnalyse_long[which(datAnalyse_long$Intervention_Phase == "pre" & datAnalyse_long$topic_Nr == "T1"),])
# #p = 0.51/0.68 (afterSpeech-analysis) --> no association; also take into account covariates:
# lme_mod <- lme(Angst_im_moment ~  Age + factor(Sex) +  factor(Condition) + relativeDwelltimeAOI, random=~1|VP_Nr, data=datAnalyse_long[which(datAnalyse_long$Intervention_Phase == "pre" & datAnalyse_long$topic_Nr == "T1"),],control=lmeControl(opt="optim"), na.action="na.omit")
# res2 <- summary(lme_mod)
# res2_anova <- Anova(lme_mod, type = "2")
# #p = 0.51/0.68 (afterSpeech-analysis)  (the same without Condition) --> no association; formally test Intervention_Phase-Interaction as well:
# lme_mod <- lme(Angst_im_moment ~  Age + factor(Sex) + factor(Condition) * factor(Intervention_Phase) * relativeDwelltimeAOI, random=~1|VP_Nr, data=datAnalyse_long[which(datAnalyse_long$topic_Nr == "T1"),],control=lmeControl(opt="optim"), na.action="na.omit")
# res2 <- summary(lme_mod)
# res2_anova <- Anova(lme_mod, type = "2")
# #all p > 0.21/0.34(afterSpeech-analysis)