## ---------------------------
##
## Script Purpose: Generate UpSet plot of main neuroimaging and cognitive associations with DNAm metrics 
##                
##
##
## ----------------------------#
## Uses results generated from the protein_vs_DNAm_global script
## https://github.com/EleanorSC/Inflammatory-DNAm_STRADL/tree/main/Results
## Specifiically, baseline model contains brain and cognitive associations with DNAm/protein signatures, where p < 0.05
## ----------------------------#
#baseline_model <- read.csv("baseline_model.csv")

# ----------------------------#
# Examine individual DNAms
# ----------------------------#

# ----------------------------#
# / END CREATE SUPPLEMENTARY TABLE
# ----------------------------#
Poor_brain_health_DNAm <- plot_global_methylation %>% 
  filter(brain_metric == "global grey matter" & estimate < 0 |
           brain_metric == "global white matter" & estimate < 0 |
           brain_metric == "total brain volume" & estimate < 0 |
           brain_metric == "global cortical thickness" & estimate < 0 |
           brain_metric == "global cortical volume" & estimate < 0 |
           brain_metric == "global cortical surface area" & estimate < 0 |
           brain_metric == "gFA" & estimate < 0 |
           brain_metric == "global subcortical volume" & estimate < 0 |
           brain_metric == "WMH" & estimate > 0 |
           brain_metric == "gMD" & estimate > 0 |
           brain_metric == "g" & estimate < 0 |
           brain_metric == "gf" & estimate < 0 |
           brain_metric == "relative brain age" & estimate > 0 |
           brain_metric == "APOE"  & estimate > 0 |
           brain_metric =="processing speed" & estimate < 0 |
           brain_metric =="executive function" & estimate < 0 |
           brain_metric =="vocabulary" & estimate < 0|
           brain_metric =="verbal declarative memory" & estimate < 0|
           brain_metric =="matrix reasoning" & estimate < 0) 

#MMP12_test <- plot_global_methylation %>% filter(DNAm == "MMP12" & brain_metric == "gMD")
CCL25_test <- plot_global_methylation %>% filter(DNAm == "CCL25" & FDR_significance == "Yes")
CCL25_test
### Which DNAm proxy has the most FDR significant hits?
Poor_brain_health_DNAm_pFDR <- Poor_brain_health_DNAm %>%
  filter(FDR_significance == "Yes") %>%
  group_by(DNAm, brain_metric) %>%
  mutate(number_significant =
           case_when(pFDR > 0.05 ~ 0,
                     TRUE ~ 1))

# ----------------------------#
# Examine individual DNAms
# ----------------------------#
sigs <- Poor_brain_health_DNAm_pFDR

### Examine the overlap for a eular plot
sigs_brainage <- sigs %>% filter(brain_metric == "relative brain age")  # n = 19 signatures
sigs_gm <- sigs %>% filter(brain_metric == "global grey matter")  # n = 15 signatures
sigs_cv <- sigs %>% filter(brain_metric == "global cortical volume")  # n = 14 signatures
sigs_ct <- sigs %>% filter(brain_metric == "global cortical thickness")  # n = 10 signatures
sigs_processingspeed <- sigs %>% filter(brain_metric == "processing speed")  # n = 11 signatures
sigs_wmh <- sigs %>% filter(brain_metric == "WMH")  # n = 9 signatures
sigs_gMD <- sigs %>% filter(brain_metric == "gMD")  # n = 6 signatures
sigs_gFA <- sigs %>% filter(brain_metric == "gFA")  # n = 7 signatures
sigs_mtr <- sigs %>% filter(brain_metric == "matrix reasoning")  # n = 7 signatures
sigs_g <- sigs %>% filter(brain_metric == "g")  # n = 1 signature
sigs_scv <- sigs %>% filter(brain_metric == "global subcortical volume")  # n = 1 signatures
sigs_wbv <- sigs %>% filter(brain_metric == "total brain volume")  # n = 2 signatures

sigs_gf<- sigs %>% filter(brain_metric == "gf") #0
sigs_verbal_declar<- sigs %>% filter(brain_metric == "verbal declarative memory") #0
sigs_exec<- sigs %>% filter(brain_metric == "executive function") #0

#gf<- c("","","","","","","","","","","","","","","","","")
#verbal_declar<- c("","","","","","","","","","","","","","","","","")
#gloal_CSA<- c("","","","","","","","","","","","","","","","","")
#vocab<- c("","","","","","","","","","","","","","","","","")
#global_WM<- c("","","","","","","","","","","","","","","","","")
#exec<- c("","","","","","","","","","","","","","","","","")
#APOE<- c("","","","","","","","","","","","","","","","","")
#global_CT<- c("","","","","","","","","","","","","","","","","")
#

relative_brain_age <-c(sigs_brainage$DNAm) # n = 19
total_brain_volume <- c(sigs_wbv$DNAm, "","","","","","","","","","","","","","","","","") # n =2 
global_GM <- c(sigs_gm$DNAm, "","","","") # n = 15
global_CV <-c(sigs_cv$DNAm, "","","","","") # n =14
global_CT<- c(sigs_ct$DNAm, "","","","","","","","","") # n 10
global_sCV<- c(sigs_scv$DNAm, "","","","","","","","","","","","","","","","","","") # n =1
processing_speed<- c(sigs_processingspeed$DNAm, "","","","","","","","") # n = 11
matrix_reasoning<- c(sigs_mtr$DNAm, "","","","","","","","","","","","") # n = 7
g<- c(sigs_g$DNAm,"","","","","","","","","","","","","","","","","","") # n =1
WMH<- c(sigs_wmh$DNAm,"","","","","","","","","","") # n = 9
gFA<- c(sigs_gFA$DNAm,"","","","","","","","","","","","") # n = 7
gMD<- c(sigs_gMD$DNAm,"","","","","","","","","","","","", "") # n =6 

dat <-  data.frame(relative_brain_age,
                   total_brain_volume,
                   global_GM,
                   global_CV,
                   global_CT,
                   global_sCV,
                   processing_speed,
                   matrix_reasoning,
                   g,WMH,gFA,gMD)

vennfun <- function(x) { 
  x$id <- seq(1, nrow(x))  #add a column of numbers (required for melt)
  xm <- reshape2::melt(x, id.vars="id", na.rm=TRUE)  #melt table into two columns (value & variable)
  xc <- reshape2::dcast(xm, value~variable, fun.aggregate=length)  #remove NA's, list presence/absence of each value for each variable (1 or 0)
  rownames(xc) <- xc$value  #value column = rownames (required for Venneuler)
  xc$value <- NULL  #remove redundent value column
  xc  #output the new dataframe
}

VennDat <- vennfun(dat)

VennDat = VennDat[-1,]

brain_metrics <- list(colnames(VennDat)) 
DNAm <- list(rownames(VennDat))

VennDat_t <- t(VennDat)

VennDat_x = VennDat[-1,]



upset(
  VennDat,
  brain_metrics,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=FALSE,
      width=0.3,
      text=list(size=7),
      mapping=aes(fill=rownames(VennDat)),
      geom_text(
        aes(label = rownames(VennDat),
            stat='count',
            position=position_fill(vjust = .5),
            na.rm=TRUE)
      )
    )
  ),
  # sort_sets=FALSE,
 # sort_sets = "descending",
 # sort_intersections = "descending" ,
  sort_sets = "ascending",
  sort_intersections = "ascending" ,
  width_ratio=0.1
) & scale_fill_manual(values = c(
  "#E5F5F9" ,"#1D91C0" ,"#67001F" ,"#F7FCFD" ,"#CB181D" ,"#78C679" ,"#F46D43" ,"#A6CEE3" ,"#FD8D3C" ,"#A6D854",
  "#D4B9DA" ,"#6A51A3" ,"#7F0000" ,"#D9D9D9" ,"#FFF7BC" ,"#000000" ,"#F0F0F0" ,"#C7EAE5" ,"#003C30" ,"#F16913",
  "#FFF7FB" ,"#8C6BB1" ,"#C7E9B4" ,"#762A83" ,"#FC9272" ,"#AE017E" ,"#F7F7F7" ,"#DF65B0" ,"#EF3B2C" ,"#74C476",
  #"indianred",
  "red", "yellow", "grey"))


# ----------------------------#
# / Upset Good brain health
# ----------------------------#
# ----------------------------#
# / END CREATE SUPPLEMENTARY TABLE
# ----------------------------#
Good_brain_health_DNAm <- plot_global_methylation %>% 
  filter(brain_metric == "global grey matter" & estimate > 0 |
           brain_metric == "global white matter" & estimate > 0 |
           brain_metric == "total brain volume" & estimate > 0 |
           brain_metric == "global cortical thickness" & estimate > 0 |
           brain_metric == "global cortical volume" & estimate > 0 |
           brain_metric == "global cortical surface area" & estimate > 0 |
           brain_metric == "gFA" & estimate > 0 |
           brain_metric == "global subcortical volume" & estimate > 0 |
           brain_metric == "WMH" & estimate < 0 |
           brain_metric == "gMD" & estimate < 0 |
           brain_metric == "g" & estimate > 0 |
           brain_metric == "gf" & estimate > 0 |
           brain_metric == "relative brain age" & estimate < 0 |
           brain_metric == "APOE"  & estimate < 0 |
           brain_metric =="processing speed" & estimate > 0 |
           brain_metric =="executive function" & estimate > 0 |
           brain_metric =="vocabulary" & estimate > 0|
           brain_metric =="verbal declarative memory" & estimate > 0|
           brain_metric =="matrix reasoning" & estimate > 0) 

### Which DNAm proxy has the most FDR significant hits?
Good_brain_health_DNAm_pFDR <- Good_brain_health_DNAm %>%
  filter(FDR_significance == "Yes") %>%
  group_by(DNAm, brain_metric) %>%
  mutate(number_significant =
           case_when(pFDR > 0.05 ~ 0,
                     TRUE ~ 1))

# ----------------------------#
# Examine individual DNAms
# ----------------------------#
sigs <- Good_brain_health_DNAm_pFDR %>% arrange(estimate)

### Examine the overlap for a eular plot
sigs_brainage <- sigs %>% filter(brain_metric == "relative brain age")  # n = 13 signatures
skimr::skim(sigs_brainage$DNAm)

sigs_gm <- sigs %>% filter(brain_metric == "global grey matter")  # n = 8 signatures
skimr::skim(sigs_gm$DNAm)

sigs_cv <- sigs %>% filter(brain_metric == "global cortical volume")  # n = 7 signatures
skimr::skim(sigs_cv$DNAm)

sigs_ct <- sigs %>% filter(brain_metric == "global cortical thickness")  # n = 12 signatures
skimr::skim(sigs_ct$DNAm)

sigs_processingspeed <- sigs %>% filter(brain_metric == "processing speed")  # n = 7 signatures
skimr::skim(sigs_processingspeed$DNAm)

sigs_wmh <- sigs %>% filter(brain_metric == "WMH")  # n = 7 signatures
skimr::skim(sigs_wmh$DNAm)

sigs_gMD <- sigs %>% filter(brain_metric == "gMD")  # n = 5 signatures
skimr::skim(sigs_gMD$DNAm)

sigs_gFA <- sigs %>% filter(brain_metric == "gFA")  # n = 2 signatures
skimr::skim(sigs_gFA$DNAm)

sigs_mtr <- sigs %>% filter(brain_metric == "matrix reasoning")  # n = 3 signatures
skimr::skim(sigs_mtr$DNAm)

sigs_g <- sigs %>% filter(brain_metric == "g")  # n = 0 signature
skimr::skim(sigs_g$DNAm)

sigs_scv <- sigs %>% filter(brain_metric == "global subcortical volume")  # n = 0 signatures
skimr::skim(sigs_scv$DNAm)

sigs_wbv <- sigs %>% filter(brain_metric == "total brain volume")  # n = 2 signatures
skimr::skim(sigs_wbv$DNAm)

sigs_gf<- sigs %>% filter(brain_metric == "gf") #0
skimr::skim(sigs_gf$DNAm)

sigs_verbal_declar<- sigs %>% filter(brain_metric == "verbal declarative memory") #0
skimr::skim(sigs_verbal_declar$DNAm)

sigs_exec<- sigs %>% filter(brain_metric == "executive function") #0
skimr::skim(sigs_exec$DNAm)

# ----------------------------#
# Put those with sig hits pFDR <0.05 into a dataframe
# ----------------------------#

relative_brain_age <-c(sigs_brainage$DNAm) # n = 13
total_brain_volume <- c(sigs_wbv$DNAm, "","","","","","","","","","", "") # n =2 
global_GM <- c(sigs_gm$DNAm, "","","","", "") # n = 8
global_CV <-c(sigs_cv$DNAm, "","","","","", "") # n = 7
global_CT<- c(sigs_ct$DNAm, "") # n 12
#global_sCV<- c(sigs_scv$DNAm, "","","","","","","","","","","","","","","","","","") # n = 0
processing_speed<- c(sigs_processingspeed$DNAm,"","","","","", "") # n = 7
matrix_reasoning<- c(sigs_mtr$DNAm, "","","","","","","","","","") # n = 3
#g<- c(sigs_g$DNAm,"","","","","","","","","","","","","","","","","","") # n =1
WMH<- c(sigs_wmh$DNAm,"","","","","", "") # n = 7
gFA<- c(sigs_gFA$DNAm,"","","","","","","","","","", "") # n = 2
gMD<- c(sigs_gMD$DNAm,"","","","","","","", "") # n =5 

dat <-  data.frame(relative_brain_age,
                   total_brain_volume,
                   global_GM,
                   global_CV,
                   global_CT,
                   #global_sCV,
                   processing_speed,
                   matrix_reasoning,
                   #g,
                   WMH,gFA,gMD)

vennfun <- function(x) { 
  x$id <- seq(1, nrow(x))  #add a column of numbers (required for melt)
  xm <- reshape2::melt(x, id.vars="id", na.rm=TRUE)  #melt table into two columns (value & variable)
  xc <- reshape2::dcast(xm, value~variable, fun.aggregate=length)  #remove NA's, list presence/absence of each value for each variable (1 or 0)
  rownames(xc) <- xc$value  #value column = rownames (required for Venneuler)
  xc$value <- NULL  #remove redundent value column
  xc  #output the new dataframe
}

VennDat <- vennfun(dat)

VennDat = VennDat[-1,]

brain_metrics <- list(colnames(VennDat)) 
DNAm <- list(rownames(VennDat))

VennDat_t <- t(VennDat)

VennDat_x = VennDat[-1,]



upset(
  VennDat,
  brain_metrics,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=FALSE,
      width=0.3,
      text=list(size=7),
      mapping=aes(fill=rownames(VennDat)),
      geom_text(
        aes(label = rownames(VennDat),
            stat='count',
            position=position_fill(vjust = .5),
            na.rm=TRUE)
      )
    )
  ),
  sort_sets = "ascending",
  sort_intersections = "ascending" ,
  width_ratio=0.1
) & scale_fill_manual(values = c(
  "#E5F5F9" ,"#1D91C0" ,"#67001F" ,"#F7FCFD" ,"#CB181D" ,"#78C679" ,"#F46D43" ,"#A6CEE3" ,"#FD8D3C" ,"#A6D854",
  "#D4B9DA" ,"#6A51A3" ,"#7F0000" ,"#D9D9D9" ,"#FFF7BC" ,"#000000" ,"#F0F0F0" ,"#C7EAE5" ,"#003C30" ,"#F16913",
  "#FFF7FB" ,"#8C6BB1" ,"#C7E9B4" ,"#762A83" ,"#FC9272" ,"#AE017E" ,"#F7F7F7" ,"#DF65B0" ,"#EF3B2C" ,"#74C476",
  #"indianred",
  "red", "yellow", "grey"))
