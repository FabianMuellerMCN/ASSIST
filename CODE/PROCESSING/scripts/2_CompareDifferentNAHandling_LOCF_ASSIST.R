#script to compare output of different intention to treat analysis-options and see if there are changes in results (significance, effect-size etc.)

rm(list=ls())
library(openxlsx) # to write excel-files
odir = ("/Volumes/dcn_assist-8/PROCESSING/output/")

# choose folders to compare
list.files(odir)
analysis_folders <- c("2020_02_24_ITT_PROTOCOLL_IPQ_etc_added", "2020_02_24_ITT_BOCF_IPQ_etc_added", "2020_02_24_ITT_LOCF_IPQ_etc_added")


rm(posthoc_res_merged)
rm(main_res_merged)
for (f in analysis_folders){
  #f = "2020_02_24_ITT_PROTOCOLL_IPQ_etc_added"
  logfile_f <- read.csv(paste0(odir, f, "/PROCESSING_MAIN_ASSIST_logfile.txt"), sep = "", header = F)
  analysisName <- strsplit(as.character(logfile_f [grep("analysis", logfile_f$V2), "V2"]), "analysis: ")[[1]][2]
  posthoc_res <- read.xlsx(paste0(odir, f, "/POSTHOC_all_DV_",  analysisName, ".xlsx"))
  posthoc_res <- posthoc_res[posthoc_res$IVs == "factor(Condition)",]
  namesOld <- names(posthoc_res)[which(!names(posthoc_res) %in% c("DV", "Intervention_Phase"))]
  names(posthoc_res)[which(!names(posthoc_res) %in% c("DV", "Intervention_Phase"))] <- paste0(namesOld,"_",analysisName ) 
  main_res =  read.xlsx(paste0(odir, f, "/MAIN_all_DV_",  analysisName, ".xlsx"))
  main_res <- main_res[which(main_res$IVs == "factor(Condition):factor(Intervention_Phase)"),]
  namesOld <- names(main_res)[which(!names(main_res) %in% "DV")]
  names(main_res)[which(!names(main_res) %in% "DV")] <- paste0(namesOld,"_",analysisName ) 
  
  if (exists ("posthoc_res_merged") == F){
    posthoc_res_merged <- posthoc_res
  } else {
    posthoc_res_merged <- merge(posthoc_res_merged, posthoc_res, by = c("DV", "Intervention_Phase"), all = T, sort=F)
  }
  if (exists ("main_res_merged") == F){
    main_res_merged <- main_res
  } else {
    main_res_merged <- merge(main_res_merged, main_res, by = c("DV"), all = T, sort=F)
  }
} # end of f-loop

#posthoc_res_merged$Intervention_Phase <- factor(posthoc_res_merged$Intervention_Phase, levels = c("pre", "post1", "post2"))
#posthoc_res_merged<- posthoc_res_merged[with(posthoc_res_merged, order(Intervention_Phase)), ]


nameSelection <- names(posthoc_res_merged) [c(grep("Pr", names(posthoc_res_merged)), 
                                              grep("t-value", names(posthoc_res_merged)),
                                              grep("r_", names(posthoc_res_merged)),
                                              grep("CohensD", names(posthoc_res_merged)))]

posthoc_res_merged <- posthoc_res_merged[,which(names(posthoc_res_merged) %in% c("DV", "Intervention_Phase",nameSelection))]
main_res_merged <- main_res_merged[,which(names(main_res_merged) %in% c("DV",nameSelection))]


# save overview
save(posthoc_res_merged, file = paste0(odir,gsub("-", "_", Sys.Date()), "_posthoc_res_merged.RData"))
write.xlsx(x = posthoc_res_merged, file = paste0(odir,gsub("-", "_", Sys.Date()), "_posthoc_res_merged.xlsx"), row.names = F)
save(main_res_merged, file = paste0(odir,gsub("-", "_", Sys.Date()), "_main_res_merged.RData"))
write.xlsx(x = main_res_merged, file = paste0(odir,gsub("-", "_", Sys.Date()), "_main_res_merged.xlsx"), row.names = F)