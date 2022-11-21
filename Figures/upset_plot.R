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
baseline_model <- read.csv("baseline_model.csv")

# ----------------------------#
# Examine individual DNAms
# ----------------------------#


#### Poor brain health
Poor_brain_health_DNAm <- baseline_model %>% filter(  omic_type == "DNAm" & brain_metric == "global grey matter" & estimate < 0 |
                                                        omic_type == "DNAm" & brain_metric == "global white matter" & estimate < 0 |
                                                        omic_type == "DNAm" & brain_metric == "total brain volume" & estimate < 0 |
                                                        omic_type == "DNAm" & brain_metric == "global cortical thickness" & estimate < 0 |
                                                        omic_type == "DNAm" & brain_metric == "global cortical volume" & estimate < 0 |
                                                        omic_type == "DNAm" & brain_metric == "global cortical surface area" & estimate < 0 |
                                                        omic_type == "DNAm" & brain_metric == "gFA" & estimate < 0 |
                                                        omic_type == "DNAm" & brain_metric == "global subcortical volume" & estimate < 0 |
                                                        omic_type == "DNAm" & brain_metric == "WMH" & estimate > 0 |
                                                        omic_type == "DNAm" & brain_metric == "gMD" & estimate > 0 |
                                                        omic_type == "DNAm" & brain_metric == "g" & estimate < 0 |
                                                        omic_type == "DNAm" & brain_metric == "gf" & estimate < 0 |
                                                        omic_type == "DNAm" & brain_metric == "relative brain age" & estimate > 0 |
                                                        omic_type == "DNAm" & brain_metric == "APOE" & estimate > 0 |
                                                        
                                                        omic_type == "DNAm" & brain_metric =="processing speed" & estimate < 0 |
                                                        omic_type == "DNAm" & brain_metric =="executive function" & estimate < 0 |
                                                        omic_type == "DNAm" & brain_metric =="vocabulary" & estimate < 0|
                                                        omic_type == "DNAm" & brain_metric =="verbal declarative memory" & estimate < 0|
                                                        omic_type == "DNAm" & brain_metric =="matrix reasoning" & estimate < 0)



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
sigs_brainage <- sigs %>% filter(brain_metric == "relative brain age")  # n = 17 signatures
sigs_gm <- sigs %>% filter(brain_metric == "global grey matter")  # n = 15 signatures
sigs_cv <- sigs %>% filter(brain_metric == "global cortical volume")  # n = 14 signatures
sigs_ct <- sigs %>% filter(brain_metric == "global cortical thickness")  # n = 10 signatures
sigs_processingspeed <- sigs %>% filter(brain_metric == "processing speed")  # n = 13 signatures
sigs_wmh <- sigs %>% filter(brain_metric == "WMH")  # n = 9 signatures
sigs_gMD <- sigs %>% filter(brain_metric == "gMD")  # n = 6 signatures
sigs_gFA <- sigs %>% filter(brain_metric == "gFA")  # n = 6 signatures
sigs_mtr <- sigs %>% filter(brain_metric == "matrix reasoning")  # n = 7 signatures
sigs_g <- sigs %>% filter(brain_metric == "g")  # n = 7 signatures
sigs_scv <- sigs %>% filter(brain_metric == "global subcortical volume")  # n = 1 signatures
sigs_wbv <- sigs %>% filter(brain_metric == "total brain volume")  # n = 2 signatures

sigs_gf<- sigs %>% filter(brain_metric == "gf") #0
sigs_verbal_declar<- sigs %>% filter(brain_metric == "verbal declarative memory") #0

gf<- c("","","","","","","","","","","","","","","","","")
verbal_declar<- c("","","","","","","","","","","","","","","","","")
gloal_CSA<- c("","","","","","","","","","","","","","","","","")
vocab<- c("","","","","","","","","","","","","","","","","")
global_WM<- c("","","","","","","","","","","","","","","","","")
exec<- c("","","","","","","","","","","","","","","","","")
APOE<- c("","","","","","","","","","","","","","","","","")
global_CT<- c("","","","","","","","","","","","","","","","","")


relative_brain_age <-c(sigs_brainage$DNAm) # n = 19
total_brain_volume <- c(sigs_wbv$DNAm, "","","","","","","","","","",
                        "","","","","","","") # n =2 
global_GM <- c(sigs_gm$DNAm, "", "","", "") # n = 15
global_CV <-c(sigs_cv$DNAm, "","","","","")
global_CT<- c(sigs_ct$DNAm, "","","","","","","", "", "")
global_sCV<- c(sigs_scv$DNAm, "","","","","","","","","","","","","","","","","", "")
processing_speed<- c(sigs_processingspeed$DNAm, "","","","","","","","")
matrix_reasoning<- c(sigs_mtr$DNAm, "","","","","","","","","","","","")
g<- c(sigs_g$DNAm,"","","","","","","","","","","","","","","","","", "")
WMH<- c(sigs_wmh$DNAm,"","","","","","","","","","")
gFA<- c(sigs_gFA$DNAm,"","","","","","","","","","","","")
gMD<- c(sigs_gFA$DNAm,"","","","","","","","","","","","")

dat <-  data.frame(relative_brain_age,total_brain_volume,global_GM,global_CV,
                   global_CT,global_sCV,processing_speed,matrix_reasoning,
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
brain_metrics <- list(colnames(VennDat)) 
DNAm <- list(rownames(VennDat))

VennDat_t <- t(VennDat)

VennDat_x = VennDat[-1,]

upset(
  VennDat,
  brain_metrics,
  
#  queries=list(
#    upset_query(set='g', fill='#6B0077'),
#    upset_query(set='processing_speed', fill='#73579B'),
#    upset_query(set='matrix_reasoning', fill='#828BBC')
#
#  ),
#  
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
  # to prevent connectors from getting the colorured
  # use `fill` instead of `color`, together with `shape='circle filled'`
#  matrix=intersection_matrix(
#    geom=geom_point(
#      shape='circle filled',
#      size=3.5,
#      stroke=0.45
#    )
#  ),
 # sort_sets=FALSE,
  sort_sets = "ascending",
  sort_intersections = "ascending" ,
  width_ratio=0.1
) & scale_fill_manual(values = c(
"#E5F5F9" ,"#1D91C0" ,"#67001F" ,"#F7FCFD" ,"#CB181D" ,"#78C679" ,"#F46D43" ,"#A6CEE3" ,"#FD8D3C" ,"#A6D854",
"#D4B9DA" ,"#6A51A3" ,"#7F0000" ,"#D9D9D9" ,"#FFF7BC" ,"#000000" ,"#F0F0F0" ,"#C7EAE5" ,"#003C30" ,"#F16913",
"#FFF7FB" ,"#8C6BB1" ,"#C7E9B4" ,"#762A83" ,"#FC9272" ,"#AE017E" ,"#F7F7F7" ,"#DF65B0" ,"#EF3B2C" ,"#74C476",
"red", "yellow", "grey"))

#& scale_fill_manual(values = c(
#     colorspace::sequential_hcl(33, 
#                                palette = "SunsetDark")))
##+
##  scale_fill_manual(values = c(
##    colorspace::sequential_hcl(101, 
##                               palette = "SunsetDark")))
#sigs_serp <- sigs %>% filter(DNAm == "SERPIND1") 
