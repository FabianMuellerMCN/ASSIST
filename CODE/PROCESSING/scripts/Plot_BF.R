get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legendBoxplot_BF <- function (i){
  #i = "Angst_im_moment"
  bf1 <- ggplot(datAnalyse, aes(col = as.factor(Intervention_Phase),x=factor(Intervention_Phase), y = datAnalyse[,i], fill = as.factor(Condition)))+
    geom_boxplot()+
    labs (title = "Subjektiv empfundene Angst in Abhängigkeit von Zeitpunkt und Gruppenzuteilung", x = "Zeitpunkt", y  = "Angstrating")+
    scale_color_manual(values=c("#C59900", "#39B600", "#00ABFD"), name = "Zeitpunkt", labels = c("vor Training 1", "nach Training 1 (3x 20 min)", "nach Training 2 (9x 20 min)"), drop=F)+
    scale_fill_manual(values=c("black","white"), name = "Gruppenzuteilung", labels = c("Treatmentgruppe", "Kontrollgruppe"), drop=F)+
    # scale_x_continuous(breaks = c(-range, -range+ticdist, -range+2*ticdist,-range+3*ticdist, -range+4*ticdist, -range+5*ticdist),
    #                    label = c("1", "2", "3", "1", "2", "3"))+
    theme_classic()
  legend<- get_legend(bf1)
  return(legend)
}# end of function
  

legendPointLinePlot_BF <- function (i){
  if ("whoRated" %in% names(datAnalyse)){
    summary_i <- summarySE(datAnalyse, measurevar= i, groupvars=c("Condition","Intervention_Phase", "whoRated"), na.rm =T)  # datAnalyse[!is.na(datAnalyse[,i]),]
  } else{
    summary_i <- summarySE(datAnalyse, measurevar= i, groupvars=c("Condition","Intervention_Phase"), na.rm =T) #  datAnalyse[!is.na(datAnalyse[,i]),]
  }
  pd <- position_dodge(0.1)
  bf1 <- ggplot(summary_i, aes(x=factor(Intervention_Phase), y= summary_i[,i], colour = as.factor(Condition), group= as.factor(Condition))) + 
    geom_point(position = pd) +
    geom_line(aes(linetype=as.factor(Condition))) +
    geom_errorbar(aes(ymin=summary_i[,i]-summary_i[,"se"], ymax=summary_i[,i]+summary_i[,"se"]), width=.1, position = pd, show.legend = F) +
    labs (y=i, x = "")+
    scale_linetype_manual(values = c(1,5), name = "Gruppenzuteilung", labels = c( "Kontrollgruppe", "Treatmentgruppe"), drop=F)+
    scale_color_manual(values=c("darkblue", "darkorange"), name = "Gruppenzuteilung", labels = c("Kontrollgruppe","Treatmentgruppe"), drop=F)+ 
    theme_classic()
  legend<- get_legend(bf1)
  return(legend)
}# end of function


boxplot_BF <- function(i){
  #range = 0.312
  #ticdist<- (range*2/5)
  #i = whichVars_plot[1]
  if ("whoRated" %in% names(datPlot) ){ # & i !="Time_s"
    summary_i <- summarySE(datPlot, measurevar= i, groupvars=c("Condition","Intervention_Phase", "whoRated"), na.rm =T)
    c1=c("Intervention_Phase", "whoRated")
    createMergeName <- T
  } else{
    summary_i <- summarySE(datPlot, measurevar= i, groupvars=c("Condition","Intervention_Phase"), na.rm =T)
    c1=c("Intervention_Phase")
    createMergeName <- F
  }
  v1=i
  datPlot_max <- aggregate(datPlot[v1],by=datPlot[c1],FUN=max, na.rm=TRUE)
  datPlot_max$Intervention_Phase_numeric <- as.numeric(datPlot_max$Intervention_Phase)
  toMerge <- res_sub_all_DV[res_sub_all_DV$IVs == "factor(Condition)" &   grepl(i, res_sub_all_DV$DV), ]
  if (createMergeName == T) {
    datPlot_max$DV = paste0(i,".",datPlot_max$whoRated)
    datPlot_max <- merge(datPlot_max, toMerge[,c("Pr(>Chisq)","Intervention_Phase", "DV")], by = c("Intervention_Phase", "DV"), all = T)
  } else {
    datPlot_max <- merge(datPlot_max, toMerge[, c("Pr(>Chisq)","Intervention_Phase")], by = "Intervention_Phase", all = T)
  }# end of createMergeName-check
  datPlot_max$label = NA
  datPlot_max[which(datPlot_max$`Pr(>Chisq)`> 0.05),"label"] <- ""
  datPlot_max[which(datPlot_max$`Pr(>Chisq)`< 0.05),"label"] <- "*"
  datPlot_max[which(datPlot_max$`Pr(>Chisq)`< 0.01),"label"] <- "**"
  datPlot_max[which(datPlot_max$`Pr(>Chisq)`< 0.001),"label"] <- "***"
  bf1 <- ggplot(datPlot, aes(col = as.factor(Intervention_Phase), x=factor(Intervention_Phase), y = datPlot[,i], fill = as.factor(Condition)))+
    geom_boxplot()+
    labs (y=i, x = "")+   # title = "Bewertung durch Gremium und Probanden", x = "Zeitpunkt",
    #ylim(min(na.omit(summary_i[,i])),max(datPlot_max[,i])+4)+
    scale_color_manual(values=c("#C59900", "#39B600", "#00ABFD"), drop=F)+
    scale_fill_manual(values=c("black","white"),drop=F)+   # name = "Gruppenzuteilung", labels = c("Treatmentgruppe", "Kontrollgruppe")
    theme_classic()+
    theme(legend.position="none")+
    geom_text(datPlot_max, mapping = aes( x=Intervention_Phase_numeric, y=datPlot_max[,i]+4, label=label, fill = NULL), col = "black")
  if ("whoRated" %in% names(datPlot)){ # & i !="Time_s"
    bf1 <- bf1 + facet_wrap(~whoRated)
  } # end of whoRated-check
  #bf1
  return(bf1)
}  # end of function


pointLinePlot_BF <- function(i){
  
  #i = "Angst_im_moment"  # primary outcome
  #i = "Time_s"
  #i= "FNEK_Score"  
  #i = "Blickkontakt"
  #i = "Angst_Verbesserung"
  #i = "Gesamturteil"
  #i = "relativeDwelltimeAOI"
  #print(i)
  if ("whoRated" %in% names(datPlot) ){  # & i !="Time_s"
    summary_i <- summarySE(datPlot, measurevar= i, groupvars=c("Condition","Intervention_Phase", "whoRated"), na.rm =T)
    c1=c("Intervention_Phase", "whoRated", "Condition")
    c2=c("Intervention_Phase", "whoRated")
    createMergeName <- T
  } else{
    summary_i <- summarySE(datPlot, measurevar= i, groupvars=c("Condition","Intervention_Phase"), na.rm =T)
    c1=c("Intervention_Phase", "Condition")
    c2=c("Intervention_Phase")
    createMergeName <- F
  }
  v1=i
  datPlot_mean <- aggregate(datPlot[v1],by=datPlot[c1],FUN=mean, na.rm=TRUE, na.action = na.omit)
  datPlot_max_mean <- aggregate(datPlot_mean[v1],by=datPlot_mean[c2],FUN=max, na.rm=TRUE)
  
  datPlot_max_mean$Intervention_Phase_numeric <- as.numeric(datPlot_max_mean$Intervention_Phase)
  toMerge <- res_sub_all_DV[res_sub_all_DV$IVs == "factor(Condition)" &   grepl(i, res_sub_all_DV$DV), ]
  if (createMergeName == T) {
    datPlot_max_mean$DV = paste0(i,".",datPlot_max_mean$whoRated)
    datPlot_max_mean <- merge(datPlot_max_mean, toMerge[,c("Pr(>Chisq)","Intervention_Phase", "DV")], by = c("Intervention_Phase", "DV"), all = T)
    } else {
      datPlot_max_mean <- merge(datPlot_max_mean, toMerge[, c("Pr(>Chisq)","Intervention_Phase")], by = "Intervention_Phase", all = T)
    }# end of createMergeName-check
  datPlot_max_mean <- merge(datPlot_max_mean, summary_i[, c(i, "se")], by = i)
  datPlot_max_mean$label = NA
  datPlot_max_mean[which(datPlot_max_mean$`Pr(>Chisq)`> 0.05),"label"] <- ""
  datPlot_max_mean[which(datPlot_max_mean$`Pr(>Chisq)`< 0.05),"label"] <- "*"
  datPlot_max_mean[which(datPlot_max_mean$`Pr(>Chisq)`< 0.01),"label"] <- "**"
  datPlot_max_mean[which(datPlot_max_mean$`Pr(>Chisq)`< 0.001),"label"] <- "***"
  
  pd <- position_dodge(0.1)
  bf1 <- ggplot(summary_i, aes(x=factor(Intervention_Phase), y= summary_i[,i], colour=as.factor(Condition), group= as.factor(Condition))) + 
    geom_errorbar(aes(ymin=summary_i[,i]-summary_i[,"se"], ymax=summary_i[,i]+summary_i[,"se"]), width=.1, position = pd) +
    geom_line(data=summary_i[!is.na(summary_i[, i]),], aes(y = summary_i[!is.na(summary_i[,i]),i],linetype=as.factor(Condition)), position = pd) +
    geom_point(position = pd) +
    labs (y=i, x = "")+
    ylim(min(na.omit(summary_i[,i]-summary_i[,"se"])), max(na.omit(summary_i[,i]+summary_i[,"se"]))*1.2)+
    #labs (y="Relative Dwell Time on Faces", x = "")+   # manual settings for paper
    scale_color_manual(values=c("darkblue", "darkorange"), drop=F)+ 
    theme_classic()+
    theme(legend.position="none")+
    #scale_y_continuous(expand = c(0, 0), limits = c(0, 0.35), breaks = c(0, 0.07, 0.14,0.21,0.28,0.35)) + # manual settings for paper
    geom_text(datPlot_max_mean, mapping = aes(x=Intervention_Phase_numeric, y=datPlot_max_mean[,i]+datPlot_max_mean[,"se"]+mean(datPlot_max_mean[,"se"], na.rm = T), label=label, fill = NULL,linetype = NULL, group= NULL), col = "black")
  if ("whoRated" %in% names(datPlot)){  # & i !="Time_s"
    bf1 <- bf1 + facet_wrap(~whoRated)
  } # end of whoRated-check
  bf1
  return(bf1)
} # end of function



