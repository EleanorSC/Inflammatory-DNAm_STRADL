## ---------------------------
##
## Script Purpose: SUBCORTICAL volume associations with DNAm metrics 
##                
##                (nb.) Neuroimaging models must contain key covariates of (site / edits / batch / estimated ICV) alongside age and sex
## 
##                 This script does the following:
##
##                (1) Creates supplementary table for pFDR significant cortical volume associations;
##                (2) Creates barplot for pFDR significant cortical volume associations facetted by whether they associate with lower or increased volumes
##
##                packages needed are tidyverse & magrittr
## ---------------------------

## ---------------------------

setwd("/Users/eleanorc_worklaptop/repos/STRADL_inflammatory_DNAm")

## ----------------------------# 
# install relevant packages

install.packages("tidyverse")
install.packages("magrittr")

library(tidyverse)
library(magrittr)

## ----------------------------#

######Load up neuroimaging full dataset & DNAm data (that was cleaned in data prep)
STRADL_FreeSurfer <- read.csv("STRADL_Measures_FreeSurfer_Main.csv")
names(STRADL_FreeSurfer)[names(STRADL_FreeSurfer) == 'id'] <-
  'stradl_ID'
PROTEOMICS_DNAm_DATA <- read.csv("PROTEOMICS_DNAm_DATA.csv")

# First load up ID linkages
ID_link_attempt <- read.csv("STRADL_DNAm_target_REM_17April2020.csv")

# Next load up methylation and proteomics data
protein_data_full <- read.csv("prot_file_150621.csv", check.names=FALSE)
annotation_data <- read.csv("Annotation_data.csv", check.names=FALSE)

# clean protein data
newnames <- c()
for(colname in names(protein_data_full)){
  if(colname %in% annotation_data$SeqId){
    egname <- annotation_data[annotation_data$SeqId == colname,]$`Entrez Gene Name`
    newnames <- append(newnames, egname)
  } else {
    newnames <- append(newnames, colname)
  }
}

print(newnames)

names(protein_data_full) <- newnames

# Load main data
names(protein_data_full)[names(protein_data_full) == 'SampleId'] <- 'stradl_ID'
protein_data_full <- protein_data_full[ -c(1:6,8:32) ]
protein_data_only <- protein_data_full[c(1:4236)]
STRADL_main_data <- protein_data_full[c(1, 4237:4286)]

# n =709
Neuroimaging_DNAm <- merge(PROTEOMICS_DNAm_DATA,
                           STRADL_FreeSurfer,
                           by = "stradl_ID")


# add in some cognitive metrics

Cognitive_data_select <- STRADL_main_data %>% 
  select(stradl_ID, digit_symbol, verbal_total, vocabulary, LM, mr_correct) %>% 
  rename(processing_speed = digit_symbol,
         executive_function = verbal_total,
         verbal_declarative_memory = LM,
         matrix_reasoning = mr_correct)


Neuroimaging_DNAm <- merge(Neuroimaging_DNAm,
                           Cognitive_data_select,
                           by = "stradl_ID")

# Create list to loop through all DNAm signatures
FULL_DNAm_list <- c(
  "CRTAM"           ,
  "EZR"             ,
  "FcRL2"           ,
  "G.CSF"           ,
  "GDF.8"           ,
  "GZMA_olink"      ,
  "N.CDase"         ,
  "NEP"             ,
  "NMNAT1"          ,
  "NTRK3_olink"    ,
  "SIGLEC1"        ,
  "SKR3"           ,
  "SMPD1"          ,
  "CCL11"          ,
  "CD6"            ,
  "CXCL10_olink"   ,
  "CXCL11_olink"   ,
  "CXCL9"          ,
  "EN.RAGE"        ,
  "FGF.21"         ,
  "HGF"            ,
  "MMP.1_olink"    ,
  "OSM"            ,
  "TGF.alpha"      ,
  "VEGFA"          ,
  "CCL21"          ,
  "MMP9"           ,
  "MPO"            ,
  "NTRK3"          ,
  "TNFRSF17"       ,
  "MIA"            ,
  "CCL25"          ,
  "IGFBP1"         ,
  "LTF"            ,
  "BCAM"           ,
  "EDA"            ,
  "C5"             ,
  "GHR"            ,
  "IGFBP4"         ,
  "CLEC11A"        ,
  "VCAM1"          ,
  "LGALS4"         ,
  "CD209"          ,
  "IL19"           ,
  "CXCL11"         ,
  "MRC2"           ,
  "CCL18"          ,
  "RETN"           ,
  "C9"             ,
  "RARRES2"        ,
  "TNFRSF1B"       ,
  "IDUA"           ,
  "ADAMTS13"       ,
  "F7"             ,
  "GNLY"           ,
  "PIGR"           ,
  "WFIKKN2"        ,
  "FCER2"          ,
  "CD48"           ,
  "CD5L"           ,
  "CNTN4"          ,
  "FCGR3B"         ,
  "SERPIND1"       ,
  "LY9"            ,
  "THBS2"          ,
  "ACY1"           ,
  "BMP1"           ,
  "TPSB2"          ,
  "GZMA"           ,
  "INSR"           ,
  "SELE"           ,
  "MPL"            ,
  "B2M"            ,
  "LTA.LTB"      ,
  "CCL22"          ,
  "CCL17"          ,
  "ADIPOQ"         ,
  "CHIT1"          ,
  "HGFAC"          ,
  "ESM1"           ,
  "CXCL10"         ,
  "PAPPA"          ,
  "SERPINA3"       ,
  "MMP2"           ,
  "CRP"            ,
  "MST1"           ,
  "ENPP7"          ,
  "C4A.C4B"        ,
  "MMP12"          ,
  "NCAM1"          ,
  "CLEC11A.1"      ,
  "SLITRK5"        ,
  "AFM"            ,
  "SELL"           ,
  "LYZ"            ,
  "MMP1"           ,
  "SHBG"           ,
  "STC1"           ,
  "GP1BA"          ,
  "LGALS3BP"      ,
  "CD163"         ,
  "FAP"           ,
  "PRSS2"         ,
  "NOTCH1"        ,
  "ICAM5"         ,
  "S100A9"        ,
  "OMD"           ,
  "SEMA3E"        ,
  "SPOCK2"
)


FULL_neuroimaging_list <- c(
  
  ### Global measures
  "global.cerebral.wm",
  "global.total.gm",
  "global.wbv",
  
  "global_cortical_surface_area",
  "global_cortical_thickness",
  "global_cortical_volume",
  "global_subcortical_volume",
  
  "gFA",
  "gMD", 
  
  "Fazekas_Score_Total"
  
)

Neuroimaging_DNAm %<>% mutate(scv.cerebellum_GM =  vol.cerebellum.lh.gm + vol.cerebellum.rh.gm,
                              scv.cerebellum_WM =  vol.cerebellum.lh.wm + vol.cerebellum.rh.wm)


Neuroimaging_DNAm %<>% mutate(global_cortical_surface_area = 
                                hem.lh.csa + hem.rh.csa,
                              global_cortical_thickness = 
                                hem.lh.ct + hem.rh.ct,
                              global_cortical_volume = 
                                hem.lh.cv + hem.rh.cv,
                              global_subcortical_volume =
                                scv.bilat.accumbens + scv.bilat.amygdala + scv.bilat.caudate +
                                scv.bilat.hippocampus + scv.bilat.pallidum + scv.bilat.putamen +
                                scv.bilat.thalamus + scv.bilat.ventraldc + scv.cerebellum_GM + scv.cerebellum_WM + vol.brainstem
                              )

# Invert the polarity of these measures
Neuroimaging_DNAm %<>% mutate(gFA =
                                gFA*-1,
                              gMD = 
                                gMD*-1)

# Converting multiple varibles into a factor
Neuroimaging_DNAm %<>% mutate_at(c("sex", "site", "edited", "batch"),
                                 as.factor)

# Making a data frame with all combinations of variables
df <-
  as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1,
                metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
model_output <- function(list) {
  
  metric <- list$metric
  DNAm <- list$DNAm
  
  regression_summary <- summary(
    lm(
      scale(Neuroimaging_DNAm[[metric]]) ~
        scale(st_age)
      + sex
      + site
      + batch
      + edited
      + scale(est.icv.BAD)
      + scale(Neuroimaging_DNAm[[DNAm]]),
      data = Neuroimaging_DNAm
    )
  )
  
  regression_summary_tidy <- broom::tidy(regression_summary)
  regression_summary_tidy_complete <- regression_summary_tidy %>% mutate(r2 = regression_summary$r.squared)
  regression_summary_tidy_DNAm <- regression_summary_tidy_complete[8, c(2, 3, 5, 6)]
  
  return(regression_summary_tidy_DNAm)
  
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = model_output)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality = 
           case_when(
             metric == "global.cerebral.wm" ~ "global",
             metric == "global.total.gm" ~ "global",
             metric == "global.wbv" ~ "global",
             
             metric == "global_cortical_surface_area" ~ "cortical",
             metric == "global_cortical_thickness"~ "cortical",
             metric == "global_cortical_volume"~ "cortical",
             metric == "global_subcortical_volume"~ "subcortical",
             
             metric == "gFA"~ "WM_integrity",
             metric == "gMD"~ "WM_integrity",
             metric == "Fazekas_Score_Total" ~ "WM_integrity",
             
             TRUE ~ "misc")) %>%
  
  # Create a brain metric column
  mutate(
    brain_metric =
      case_when(
        
        ### Global
        metric == "global.cerebral.wm" ~ "global white matter",
        metric == "global.total.gm" ~ "global grey matter",
        metric == "global.wbv" ~ "total brain volume",
        
        metric == "global_cortical_surface_area" ~ "global cortical surface area",
        metric == "global_cortical_thickness"~ "global cortical thickness",
        metric == "global_cortical_volume"~ "global cortical volume",
        metric == "global_subcortical_volume"~ "global subcortical volume",
        
        metric == "gFA"~ "gFA",
        metric == "gMD"~ "gMD",
        
        metric == "Fazekas_Score_Total" ~ "WMH",
        
        TRUE ~ "misc"
      )
  ) %>%
  
  # Create pFDR column
  group_by(brain_metric, modality) %>%
  mutate(pFDR = p.adjust(p.value, method = "fdr")) %>%
  # Create a column to denote where there are significant hits
  mutate(significance =
           case_when(p.value < 0.05 ~ "Yes",
                     TRUE ~ "No")) %>%
  # Create a column to denote where there are FDR significant hits
  mutate(FDR_significance =
           case_when(pFDR < 0.05 ~ "Yes",
                     TRUE ~ "No")) %>%
  # Create a column to denote whether these are DNAm or proteins
  mutate(omic_type = "DNAm") %>%
  mutate(model = "Model 1 (n=709)")

###
plot_global_neuroimaging_methylation <- newdf



# ----------------------------#
# Run cognitive associations
# ----------------------------#

FULL_cognitive_list <- c("g", "gf", "Brain_age", "APOE",
                         "processing_speed",
                         "executive_function",
                         "vocabulary",
                         "verbal_declarative_memory",
                         "matrix_reasoning")

# Making a data frame with all combinations of variables
df <-
  as.data.frame(expand.grid(FULL_DNAm_list, FULL_cognitive_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1,
                metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
model_output <- function(list) {
  
  metric <- list$metric
  DNAm <- list$DNAm
  
  regression_summary <- summary(
    lm(
      scale(Neuroimaging_DNAm[[metric]]) ~
        scale(st_age)
      + sex
     # + site
     # + batch
     # + edited
     # + scale(est.icv.BAD)
      + scale(Neuroimaging_DNAm[[DNAm]]),
      data = Neuroimaging_DNAm
    )
  )
  
  regression_summary_tidy <- broom::tidy(regression_summary)
  regression_summary_tidy_complete <- regression_summary_tidy %>% mutate(r2 = regression_summary$r.squared)
  regression_summary_tidy_DNAm <- regression_summary_tidy_complete[4, c(2, 3, 5, 6)]
  
  return(regression_summary_tidy_DNAm)
  
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = model_output)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality = 
           case_when(
             metric == "g" ~ "cognitive",
             metric == "gf" ~ "cognitive",
             metric == "APOE" ~ "genetic",
             metric == "Brain_age" ~ "brain_ageing",
             
             metric =="processing_speed"~ "cognitive",
             metric =="executive_function"~ "cognitive",
             metric =="vocabulary"~ "cognitive",
             metric =="verbal_declarative_memory"~ "cognitive",
             metric =="matrix_reasoning"~ "cognitive",
             
             TRUE ~ "misc")) %>%
  
  # Create a brain metric column
  mutate(
    brain_metric =
      case_when(
        
        metric == "g" ~ "g",
        metric == "gf" ~ "gf",
        metric == "APOE" ~ "APOE",
        metric == "Brain_age" ~ "relative brain age",
        
        metric =="processing_speed"~ "processing speed",
        metric =="executive_function"~ "executive function",
        metric =="vocabulary"~ "vocabulary",
        metric =="verbal_declarative_memory"~ "verbal declarative memory",
        metric =="matrix_reasoning"~ "matrix reasoning",
        
        TRUE ~ "misc"
      )
  ) %>%
  
  # Create pFDR column
  group_by(brain_metric, modality) %>%
  mutate(pFDR = p.adjust(p.value, method = "fdr")) %>%
  # Create a column to denote where there are significant hits
  mutate(significance =
           case_when(p.value < 0.05 ~ "Yes",
                     TRUE ~ "No")) %>%
  # Create a column to denote where there are FDR significant hits
  mutate(FDR_significance =
           case_when(pFDR < 0.05 ~ "Yes",
                     TRUE ~ "No")) %>%
  # Create a column to denote whether these are DNAm or proteins
  mutate(omic_type = "DNAm") %>%
  mutate(model = "Model 1 (n=709)")

###
plot_cognitive_methylation <- newdf


# ----------------------------#
# Combine cognitive and brain analyses
# ----------------------------#

plot_global_methylation <- rbind(plot_global_neuroimaging_methylation,
                                 plot_cognitive_methylation)

# ----------------------------#
# Table of significant cortical volume regressions for supplementary document
# ----------------------------#

table_model1 <- plot_global_methylation %>%
  filter(FDR_significance == "Yes") %>%
  group_by(DNAm, brain_metric) %>%
  arrange(brain_metric,
          estimate,
          by_group = TRUE) %>%
  mutate(CI_lower =
           estimate - (1.96 * std.error),
         CI_upper =
           estimate + (1.96 * std.error)) %>%
  select(brain_metric, DNAm, estimate, CI_lower, CI_upper, r2, p.value, pFDR)

# ----------------------------#
# Examine which DNAm associate across multiple cortical regions
# ----------------------------#

### Which DNAm proxy has the most FDR significant hits?
sigs <- plot_global_methylation %>%
  #filter(estimate < 0 & FDR_significance == "Yes") %>%
  filter(
    brain_metric == "global grey matter" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "global white matter" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "total brain volume" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "global cortical thickness" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "global cortical volume" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "global cortical surface area" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "gFA" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "global subcortical volume" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "WMH" & estimate > 0 & FDR_significance == "Yes" |
      brain_metric == "gMD" & estimate > 0 & FDR_significance == "Yes" |
      brain_metric == "g" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "gf" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "relative brain age" & estimate > 0 & FDR_significance == "Yes" |
      brain_metric == "APOE" & estimate > 0 & FDR_significance == "Yes" |
      brain_metric == "processing speed" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "executive function" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "vocabulary" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "verbal declarative memory" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "matrix reasoning" & estimate < 0 & FDR_significance == "Yes"
  ) %>%
  group_by(DNAm, brain_metric) %>%
  mutate(number_significant =
           case_when(FDR_significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_sigs <-
  aggregate(sigs$number_significant,
            by = list(DNAm = sigs$DNAm),
            FUN = sum) %>%
  arrange(desc(x)) %>%
  rename(n_pFDR = x)


# ----------------------------#
# EULAR BEGIN
## ----------------------------#
#
#### Examine MMP12 as it associates across 9 different brain / cog metrics
#sigs_MMP12 <- sigs %>% filter(DNAm == "MMP12")
#
#
#### Examine the overlap for a eular plot
#sigs_brainage <- sigs %>% filter(brain_metric == "relative brain age")  # n = 17 signatures
#sigs_gm <- sigs %>% filter(brain_metric == "global grey matter")  # n = 15 signatures
#sigs_cv <- sigs %>% filter(brain_metric == "global cortical volume")  # n = 14 signatures
#sigs_ct <- sigs %>% filter(brain_metric == "global cortical thickness")  # n = 10 signatures
#sigs_processingspeed <- sigs %>% filter(brain_metric == "processing speed")  # n = 13 signatures
#sigs_wmh <- sigs %>% filter(brain_metric == "WMH")  # n = 9 signatures
#sigs_gMD <- sigs %>% filter(brain_metric == "gMD")  # n = 6 signatures
#sigs_gFA <- sigs %>% filter(brain_metric == "gFA")  # n = 6 signatures
#sigs_mtr <- sigs %>% filter(brain_metric == "matrix reasoning")  # n = 7 signatures
#sigs_g <- sigs %>% filter(brain_metric == "g")  # n = 7 signatures
#sigs_scv <- sigs %>% filter(brain_metric == "global subcortical volume")  # n = 1 signatures
#sigs_wbv <- sigs %>% filter(brain_metric == "total brain volume")  # n = 2 signatures
#
#
##overlap 
#
#
#
##overlap
#overlap <- list(
#  "relative_brain_age" = sigs_brainage$DNAm,
#  "total brain volume" = sigs_wbv$DNAm,
#  "global_GM" = sigs_gm$DNAm,
#  "global_CV" = sigs_cv$DNAm,
#  "global_CT" = sigs_ct$DNAm,
#  "global_sCV" = sigs_scv$DNAm,
#  "processing_speed" = sigs_processingspeed$DNAm,
#  "matrix_reasoning" = sigs_mtr$DNAm,
#  "g" = sigs_g$DNAm,
#  "WMH" = sigs_wmh$DNAm,
#  "gFA" = sigs_gFA$DNAm,
#  "gMD" = sigs_gFA$DNAm
#  ) %>% 
#  map2_dfr(., names(.), ~tibble::enframe(.x) %>% mutate(group=.y)) %>% 
#  add_count(value) %>%  
#  group_by(group) %>% 
#  rename(n_brain_metric = name) %>% 
#  mutate(overlaps = case_when(n > 1 ~ "multiple",
#                              TRUE ~ "single")) %>%
#  filter(overlaps == "multiple") %>%
#  arrange(value) 
#  #summarise(group2 = ifelse(n() > 1, "multiple", group)) %>% 
#  count(group2)
#  
####### overlap?
#  relative_brain_age <- sigs_brainage$DNAm
#  total_brain_volume <- sigs_wbv$DNAm
#  global_GM <- sigs_gm$DNAm
#  global_CV <-  sigs_cv$DNAm
#  global_CT <- sigs_ct$DNAm
#  global_sCV <-  sigs_scv$DNAm
#  processing_speed <-  sigs_processingspeed$DNAm
#  matrix_reasoning <-  sigs_mtr$DNAm
#  g <- sigs_g$DNAm
#  WMH <- sigs_wmh$DNAm
#  gFA <- sigs_gFA$DNAm
#  gMD <- sigs_gFA$DNAm  
#  
##baseoverlap <- intersect(global_GM, global_CV)
#  
#  install.packages("venneuler")
#  library(venneuler)
#  
#  
#  H <- list(
#    "relative_brain_age" = sigs_brainage$DNAm,
#    "total brain volume" = sigs_wbv$DNAm,
#    "global_GM" = sigs_gm$DNAm,
#    "global_CV" = sigs_cv$DNAm,
#    "global_CT" = sigs_ct$DNAm,
#    "global_sCV" = sigs_scv$DNAm,
#    "processing_speed" = sigs_processingspeed$DNAm,
#    "matrix_reasoning" = sigs_mtr$DNAm,
#    "g" = sigs_g$DNAm,
#    "WMH" = sigs_wmh$DNAm,
#    "gFA" = sigs_gFA$DNAm,
#    "gMD" = sigs_gFA$DNAm)
#  
#  Hmat <- lapply(seq_along(H), function (i) {
#    Set <- names(H)[i]
#    data.frame(Element = H[[i]], Set)
#  })
#  Hmat <- do.call(rbind, Hmat)
#  vd <- venneuler(Hmat)
#  plot(vd, legend = TRUE, quantities= TRUE)
#  
#  plot(venneuler::venneuler(dat))
#  
#  #########
#  
#  relative_brain_age <- c(sigs_brainage$DNAm)
#  total_brain_volume <- c(sigs_wbv$DNAm)
#  global_GM <- c(sigs_gm$DNAm)
#  global_CV<- c(sigs_cv$DNAm)
#  global_CT<- c(sigs_ct$DNAm)
#  global_sCV<- c(sigs_scv$DNAm)
#  processing_speed<- c(sigs_processingspeed$DNAm)
#  matrix_reasoning<- c(sigs_mtr$DNAm)
#  g<- c(sigs_g$DNAm)
#  WMH<- c(sigs_wmh$DNAm)
#  gFA<- c(sigs_gFA$DNAm)
#  gMD<- c(sigs_gFA$DNAm)
#  
#  
#  relative_brain_age <-
#    c(
#      "NEP",
#      "SKR3",
#      "FGF.21",
#      "MMP.1_olink",
#      "VEGFA",
#      "IGFBP4",
#      "CCL18",
#      "PIGR",
#      "THBS2",
#      "ACY1" ,
#      "SELE",
#      "CRP",
#      "MMP12",
#      "MMP1",
#      "STC1",
#      "LGALS3BP",
#      "PRSS2"
#    )
#  
#  total_brain_volume <- c("IGFBP4",
#                          "PIGR",
#                          "",
#                          "",
#                          "",
#                          "",
#                          "",
#                          "",
#                          "",
#                          "",
#                          "",
#                          "",
#                          "",
#                          "",
#                          "",
#                          "",
#                          "")
#  
#  global_GM <-
#    c(
#      "SKR3"   ,
#      "FGF.21" ,
#      "HGF"   ,
#      "VEGFA"  ,
#      "MMP9"   ,
#      "IGFBP4" ,
#      "RARRES2"  ,
#      "PIGR"  ,
#      "SERPIND1",
#      "ICAM5" ,
#      "THBS2"  ,
#      "CCL17"  ,
#      "CRP"   ,
#      "MMP12"  ,
#      "PRSS2",
#      "",
#      ""
#    )
#  
#  global_CV <-
#    c(
#      "SKR3"   ,
#      "HGF"    ,
#      "TGF.alpha" ,
#      "VEGFA"  ,
#      "MMP9"  ,
#      "THBS2"  ,
#      "CCL17"  ,
#      "MMP12"     ,
#      "PRSS2"  ,
#      "ICAM5",
#      "IGFBP4"  ,
#      "RARRES2"  ,
#      "PIGR"  ,
#      "SERPIND1" ,
#      "",
#      "",
#      ""
#    )
#  
#  global_CT<- c(sigs_ct$DNAm, "","","","","","","")
#  
#  global_sCV<- c(sigs_scv$DNAm, "","","","","","","","","","","","","","","","")
#  
#  processing_speed<- c(sigs_processingspeed$DNAm, "","","","")
#  
#  matrix_reasoning<- c(sigs_mtr$DNAm, "","","","","","","","","","")
#  
#  g<- c(sigs_g$DNAm,"","","","","","","","","","","","","","")
#  
#  WMH<- c(sigs_wmh$DNAm,"","","","","","","","")
#  
#  gFA<- c(sigs_gFA$DNAm,"","","","","","","","","","")
#  
#  gMD<- c(sigs_gFA$DNAm,"","","","","","","","","","")
#  
#  
#  dat <-  data.frame(relative_brain_age,total_brain_volume,global_GM,global_CV,
#                     global_CT,global_sCV,processing_speed,matrix_reasoning,
#                     g,WMH,gFA,gMD)
#  
#  dat2 <-  data.frame(relative_brain_age, total_brain_volume, global_GM, global_CV)
#  
#  dat3 <-  data.frame(processing_speed,matrix_reasoning,g,relative_brain_age)
#  
#  dat4 <-  data.frame(relative_brain_age, total_brain_volume, global_GM,g)
#  
#  dat5 <-  data.frame(relative_brain_age, total_brain_volume, global_GM,processing_speed)
#  
# 
#  dat6 <-  data.frame(global_CV, global_sCV, g, global_GM,total_brain_volume )
#  
#  vennfun <- function(x) { 
#    x$id <- seq(1, nrow(x))  #add a column of numbers (required for melt)
#    xm <- reshape2::melt(x, id.vars="id", na.rm=TRUE)  #melt table into two columns (value & variable)
#    xc <- reshape2::dcast(xm, value~variable, fun.aggregate=length)  #remove NA's, list presence/absence of each value for each variable (1 or 0)
#    rownames(xc) <- xc$value  #value column = rownames (required for Venneuler)
#    xc$value <- NULL  #remove redundent value column
#    xc  #output the new dataframe
#  }
#  
#  VennDat <- vennfun(dat6)
#  
#  genes.venn <- eulerr::euler(VennDat, shape = "ellipse")
#  
#  plot(genes.venn,
#       quantities = TRUE,
#      # labels = TRUE,
#       legend = list(labels = c("cortical volume", 
#                                "subcortical volume",
#                                "g", 
#                                "GM volume",
#                                "wbv"))
#    
#      )
#  
#  # ----------------------------#
#  # EULAR END
#  # ----------------------------#
#  
#  myV <- nVennR::plotVenn(as.list(dat2))
#
#  
#  
#  myV4 <-nVennR::plotVenn(list(a=c(1, 2, 3), b=c(3, 4, 5), c=c(3, 6, 1)), 
#                   nCycles = 2000, 
#                   setColors=c('red', 'green', 'blue'), 
#                   labelRegions=F, fontScale=2, opacity=0.2, borderWidth=2)
#  
 #install.packages("eulerr")
 #library(eulerr)
  
    

# 
# install.packages("VennDiagram")                           
# library(VennDiagram)  
#  
#overlap2 <- VennDiagram::calculate.overlap(
#  x = list(
#    "relative_brain_age" = sigs_brainage$DNAm,
#    "total brain volume" = sigs_wbv$DNAm,
#    "global_GM" = sigs_gm$DNAm,
#    "global_CV" = sigs_cv$DNAm,
#    "global_CT" = sigs_ct$DNAm,
#    "global_sCV" = sigs_scv$DNAm,
#    "processing_speed" = sigs_processingspeed$DNAm,
#    "matrix_reasoning" = sigs_mtr$DNAm,
#    "g" = sigs_g$DNAm,
#    "WMH" = sigs_wmh$DNAm,
#    "gFA" = sigs_gFA$DNAm,
#    "gMD" = sigs_gFA$DNAm
#  )
#)

### eular
#ULAR <-  c(relative_brain_age = 17, 
#           global_GM = 15,
#           global_CV = 14,
#           global_CT = 10,
#           global_sCV = 1,
#           
#           g = 3,
#           processing_speed = 13, 
#           matrix_reasoning = 7,
#           
#           WMH = 9, 
#           gMD = 6,
#           gMD = 7,
#           
#           "relative_brain_age&global_GM" = 9, #SKR3, FGF.21,VEGFA,IGFBP4,PIGR,THBS2,	CRP, MMP12,PRSS2
#           "global_GM&global_CV" = 4, #CCL17 HGF ICAM5 IGFBP4 MMP12
#           "relative_brain_age&matrix_reasoning" = 4, #CCL18 SKR3 CRP IGFBP4
#           "relative_brain_age&global_CT" = 3, # CRP FGF.21
#           

#
# ----------------------------#
# PLOT - barplot of significant and FDR significant for good brain health measures
# ----------------------------#
plot2 <- plot_global_methylation


## Separate into direction of effect where increases in DNAm associate with poor brain health outcomes
## Excepting gMD and WMH which will show increases in effect size corresponding to poor brain health metrics...

Good_brain_health <- plot2 %>% filter(brain_metric == "global grey matter" & estimate > 0 |
                                        brain_metric == "global white matter" & estimate > 0 |
                                        brain_metric == "total brain volume" & estimate > 0 |
                                        brain_metric == "global cortical thickness" & estimate > 0 |
                                        brain_metric == "global cortical volume" & estimate > 0 |
                                        brain_metric == "global cortical surface area" & estimate > 0 |
                                        brain_metric == "gFA" & estimate > 0 |
                                        brain_metric == "global subcortical volume" & estimate > 0 |
                                        brain_metric == "WMH" & estimate > 0 |
                                        brain_metric == "gMD" & estimate > 0 |
                                        brain_metric == "g" & estimate > 0 |
                                        brain_metric == "gf" & estimate > 0 |
                                        brain_metric == "relative brain age" & estimate < 0 |
                                        brain_metric == "APOE" & estimate < 0 |
                                      
                                      brain_metric =="processing speed" & estimate > 0 |
                                      brain_metric =="executive function" & estimate > 0 |
                                      brain_metric =="vocabulary" & estimate > 0|
                                      brain_metric =="verbal declarative memory" & estimate > 0|
                                      brain_metric =="matrix reasoning" & estimate > 0
                                      
                                      
                                        )

### Which DNAm proxy has the most FDR significant hits?
Good_brain_health %<>%
  mutate(number_significant =
           case_when(significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Good_brain_health <-
  aggregate(
    Good_brain_health$number_significant,
    by = list(brain_metric = Good_brain_health$brain_metric),
    FUN = sum
  )

Top_hits_Good_brain_health %<>% arrange(desc(x)) %>% rename(p = x)

# Run FDR
Good_brain_health %<>%
  mutate(number_significant_FDR =
           case_when(FDR_significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Good_brain_health_2 <-
  aggregate(
    Good_brain_health$number_significant_FDR,
    by = list(brain_metric = Good_brain_health$brain_metric),
    FUN = sum
  ) %>% rename(pFDR = x)


Good_brain <- merge(Top_hits_Good_brain_health,
                    Top_hits_Good_brain_health_2,
                    by = "brain_metric") %>%
  mutate(direction =
           "Favourable brain health outcome")

#### Poor brain health
Poor_brain_health <- plot2 %>% filter(brain_metric == "global grey matter" & estimate < 0 |
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
                                        brain_metric == "APOE" & estimate > 0 |
                                        
                                        brain_metric =="processing speed" & estimate < 0 |
                                        brain_metric =="executive function" & estimate < 0 |
                                        brain_metric =="vocabulary" & estimate < 0|
                                        brain_metric =="verbal declarative memory" & estimate < 0|
                                        brain_metric =="matrix reasoning" & estimate < 0)

### Which DNAm proxy has the most FDR significant hits?
Poor_brain_health %<>%
  mutate(number_significant =
           case_when(significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Poor_brain_health <-
  aggregate(
    Poor_brain_health$number_significant,
    by = list(brain_metric = Poor_brain_health$brain_metric),
    FUN = sum
  )

Top_hits_Poor_brain_health %<>% arrange(desc(x)) %>% rename(p = x)

# Run FDR
Poor_brain_health %<>%
  mutate(number_significant_FDR =
           case_when(FDR_significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Poor_brain_health_2 <-
  aggregate(
    Poor_brain_health$number_significant_FDR,
    by = list(brain_metric = Poor_brain_health$brain_metric),
    FUN = sum
  ) %>% rename(pFDR = x)


Poor_brain <- merge(Top_hits_Poor_brain_health,
                    Top_hits_Poor_brain_health_2,
                    by = "brain_metric") %>%
  mutate(direction =
           "Poor brain health outcome")

Significance_plot <- rbind(Poor_brain, Good_brain)

# ----------------------------#
# Code for Barplot of count of significant hits
# ----------------------------#

####
ggplot(Significance_plot,
       aes(
         y = reorder(brain_metric, p),
         x = p,
         fill = reorder(brain_metric, p)
       )) +
  
  geom_col(aes(y = reorder(brain_metric, p),
               x = p)) +
  
  geom_col(aes(
    y = reorder(brain_metric, p),
    x = pFDR,
    col = "grey"
  )) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",
    
    strip.text = element_text(
      size = 6,
      face = "bold",
      family = "sans",
      colour = "black"
    ),
    axis.text.x = element_text(
      vjust = 0.5,
      hjust = 1,
      size = 6,
      family = "sans"
    ),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(
      size = 7,
      face = "bold",
      family = "sans"
    ),
    axis.title.x = element_text(
      size = 7,
      face = "bold",
      family = "sans"
    )
  ) +
  
  labs(x = "number of significant assocations with DNAm signatures",
       y = "global brain metrics") +
  
  facet_wrap( ~ direction) +
  scale_fill_manual(values = c(
   # viridis::magma(n = 19)
    colorspace::sequential_hcl(19, palette = "SunsetDark")
    )) +
  scale_colour_manual(values = c("#808080"))
  
# scale_fill_manual(values = c(
#   "#666699",
#   # "#7F8D9C", # cerebellum
#   # accumbens
#   "#4C1D4BFF", #cerebellum GM
#   "#751F58FF", #putamen
#   "#A11A5BFF", #caudate
#   "#AA0144",
#   "#CB1B4FFF", #palli
#   "#E83F3FFF", #amy
#   "#F2704DFF", #thala
#   "#F69C73FF", #vdc
#   "#F7C5A5FF", #brainstem
#   "#FAEBDDFF" # hippo
#   
#   #viridis::rocket(n = 11)
# )) +
# 
# scale_colour_manual(values = c("#808080")) 
#+
#  xlim(0, 27)

## ----------------------------# 
# MODEL 3 - LIFESTYLE
## ----------------------------#

## ----------------------------# 
# MODEL 3 - LIFESTYLE
## ----------------------------#

phenotypes_STRADL <- read.csv("phenotypes.csv")
ID_link_attempt <- read.csv("STRADL_DNAm_target_REM_17April2020.csv")
GS_link <- ID_link_attempt %>% select("stradl_ID", "GS_id")
names(phenotypes_STRADL)[names(phenotypes_STRADL) == 'id'] <- 'GS_id'
phenotypes_STRADL <- merge(GS_link, phenotypes_STRADL, by = "GS_id")


sensitivity_analysis_lifestyle <- phenotypes_STRADL %>% select(stradl_ID, 
                                                               units, 
                                                               drink_status, 
                                                               hypertension_category, 
                                                               bmi, 
                                                               whr, 
                                                               body_fat)
#STRADL_main_data
### Clean up this lifestyle data
sensitivity_analysis_lifestyle %<>% 
  mutate(hypertension = case_when(
    hypertension_category == 1 ~ 1,
    hypertension_category == 2 ~ 1,
    hypertension_category == 3 ~ 1,
    TRUE ~ 0),
    CurrentDrinker = case_when(
      drink_status == 1 ~ 1,
      TRUE ~ 0))

# n = 778
Lifestyle_DNAm <- merge(PROTEOMICS_DNAm_DATA, 
                        sensitivity_analysis_lifestyle,
                        by = "stradl_ID")

#####
Neuroimaging_DNAm_lifestyle <- merge(Lifestyle_DNAm, 
                                     STRADL_FreeSurfer,
                                     by = "stradl_ID")


Neuroimaging_DNAm_lifestyle %<>% mutate(scv.cerebellum_GM =  vol.cerebellum.lh.gm + vol.cerebellum.rh.gm,
                                        scv.cerebellum_WM =  vol.cerebellum.lh.wm + vol.cerebellum.rh.wm)


Neuroimaging_DNAm_lifestyle %<>% mutate(global_cortical_surface_area = 
                                          hem.lh.csa + hem.rh.csa,
                                        global_cortical_thickness = 
                                          hem.lh.ct + hem.rh.ct,
                                        global_cortical_volume = 
                                          hem.lh.cv + hem.rh.cv,
                                        global_subcortical_volume =
                                          scv.bilat.accumbens + scv.bilat.amygdala + scv.bilat.caudate +
                                          scv.bilat.hippocampus + scv.bilat.pallidum + scv.bilat.putamen +
                                          scv.bilat.thalamus + scv.bilat.ventraldc + scv.cerebellum_GM + scv.cerebellum_WM + vol.brainstem
)

# Invert the polarity of these measures
Neuroimaging_DNAm_lifestyle %<>% mutate(gFA =
                                          gFA*-1,
                                        gMD = 
                                          gMD*-1)


# Making a data frame with all combinations of variables
df <-
  as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1,
                metric = Var2)

# Converting multiple varibles into a factor
Neuroimaging_DNAm_lifestyle %<>% mutate_at(c("sex", "site", "edited", "batch",
                                             "hypertension", "CurrentSmoker", "CurrentDrinker"),
                                           as.factor)

# Making a data frame with all combinations of variables
df <- as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1, metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
model_output <- function(list) {
  
  metric <- list$metric
  DNAm <- list$DNAm
  
  regression_summary <- summary(
    lm(
      scale(Neuroimaging_DNAm_lifestyle[[metric]]) ~ 
        scale(st_age)
      + sex
      + site
      + batch
      + edited
      
      + hypertension
      + CurrentSmoker
      + CurrentDrinker
      + scale(bmi)
      
      + scale(est.icv.BAD)
      + scale(Neuroimaging_DNAm_lifestyle[[DNAm]]),
      data = Neuroimaging_DNAm_lifestyle
    )
  )
  
  regression_summary_tidy <- broom::tidy(regression_summary)
  regression_summary_tidy_complete <- regression_summary_tidy %>% mutate(r2 = regression_summary$r.squared)
  regression_summary_tidy_DNAm <- regression_summary_tidy_complete[12, c(2, 3, 5, 6)]
  
  return(regression_summary_tidy_DNAm)
  
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = model_output)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality = 
           case_when(
             metric == "global.cerebral.wm" ~ "global",
             metric == "global.total.gm" ~ "global",
             metric == "global.wbv" ~ "global",
             
             metric == "global_cortical_surface_area" ~ "cortical",
             metric == "global_cortical_thickness"~ "cortical",
             metric == "global_cortical_volume"~ "cortical",
             metric == "global_subcortical_volume"~ "subcortical",
             
             metric == "gFA"~ "WM_integrity",
             metric == "gMD"~ "WM_integrity",
             metric == "Fazekas_Score_Total" ~ "WM_integrity",
             
             TRUE ~ "misc")) %>%
  
  # Create a brain metric column
  mutate(
    brain_metric =
      case_when(
        
        ### Global
        metric == "global.cerebral.wm" ~ "global white matter",
        metric == "global.total.gm" ~ "global grey matter",
        metric == "global.wbv" ~ "total brain volume",
        
        metric == "global_cortical_surface_area" ~ "global cortical surface area",
        metric == "global_cortical_thickness"~ "global cortical thickness",
        metric == "global_cortical_volume"~ "global cortical volume",
        metric == "global_subcortical_volume"~ "global subcortical volume",
        
        metric == "gFA"~ "gFA",
        metric == "gMD"~ "gMD",
        
        metric == "Fazekas_Score_Total" ~ "WMH",
        
        TRUE ~ "misc"
      )
  ) %>%
  
  # Create pFDR column
  group_by(brain_metric, modality) %>%
  mutate(pFDR = p.adjust(p.value, method = "fdr")) %>%
  # Create a column to denote where there are significant hits
  mutate(significance =
           case_when(p.value < 0.05 ~ "Yes",
                     TRUE ~ "No")) %>%
  # Create a column to denote where there are FDR significant hits
  mutate(FDR_significance =
           case_when(pFDR < 0.05 ~ "Yes",
                     TRUE ~ "No")) %>%
  # Create a column to denote whether these are DNAm or proteins
  mutate(omic_type = "DNAm") %>%
  mutate(model = "Model 3 Lifestyle (n=709)")

###
plot_global_neuroimaging_methylation_lifestyle <- newdf



# ----------------------------#
# Run cognitive associations
# ----------------------------#

FULL_cognitive_list <- c("g", "gf", "Brain_age", "APOE",
                         "processing_speed",
                         "executive_function",
                         "vocabulary",
                         "verbal_declarative_memory",
                         "matrix_reasoning")

# Making a data frame with all combinations of variables
df <-
  as.data.frame(expand.grid(FULL_DNAm_list, FULL_cognitive_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1,
                metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
model_output <- function(list) {
  
  metric <- list$metric
  DNAm <- list$DNAm
  
  regression_summary <- summary(
    lm(
      scale(Neuroimaging_DNAm_lifestyle[[metric]]) ~
        scale(st_age)
      + sex
      #+ site
      #+ batch
      #+ edited
      #+ scale(est.icv.BAD)
      + scale(Neuroimaging_DNAm_lifestyle[[DNAm]]),
      data = Neuroimaging_DNAm_lifestyle
    )
  )
  
  regression_summary_tidy <- broom::tidy(regression_summary)
  regression_summary_tidy_complete <- regression_summary_tidy %>% mutate(r2 = regression_summary$r.squared)
  regression_summary_tidy_DNAm <- regression_summary_tidy_complete[4, c(2, 3, 5, 6)]
  
  return(regression_summary_tidy_DNAm)
  
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = model_output)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality = 
           case_when(
             metric == "g" ~ "cognitive",
             metric == "gf" ~ "cognitive",
             metric == "APOE" ~ "genetic",
             metric == "Brain_age" ~ "brain_ageing",
             
             metric =="processing_speed"~ "cognitive",
             metric =="executive_function"~ "cognitive",
             metric =="vocabulary"~ "cognitive",
             metric =="verbal_declarative_memory"~ "cognitive",
             metric =="matrix_reasoning"~ "cognitive",
             
             TRUE ~ "misc")) %>%
  
  # Create a brain metric column
  mutate(
    brain_metric =
      case_when(
        
        metric == "g" ~ "g",
        metric == "gf" ~ "gf",
        metric == "APOE" ~ "APOE",
        metric == "Brain_age" ~ "relative brain age",
        
        metric =="processing_speed"~ "processing speed",
        metric =="executive_function"~ "executive function",
        metric =="vocabulary"~ "vocabulary",
        metric =="verbal_declarative_memory"~ "verbal declarative memory",
        metric =="matrix_reasoning"~ "matrix reasoning",
        
        TRUE ~ "misc"
      )
  ) %>%
  
  # Create pFDR column
  group_by(brain_metric, modality) %>%
  mutate(pFDR = p.adjust(p.value, method = "fdr")) %>%
  # Create a column to denote where there are significant hits
  mutate(significance =
           case_when(p.value < 0.05 ~ "Yes",
                     TRUE ~ "No")) %>%
  # Create a column to denote where there are FDR significant hits
  mutate(FDR_significance =
           case_when(pFDR < 0.05 ~ "Yes",
                     TRUE ~ "No")) %>%
  # Create a column to denote whether these are DNAm or proteins
  mutate(omic_type = "DNAm") %>%
  mutate(model = "Model 3 Lifestyle (n=709)")

###
plot_cognitive_methylation_lifestyle <- newdf

# ----------------------------#
# Combine cognitive and brain analyses
# ----------------------------#

plot_global_methylation <- rbind(plot_global_neuroimaging_methylation_lifestyle,
                                 plot_cognitive_methylation_lifestyle)

# ----------------------------#
# Table of significant cortical volume regressions for supplementary document
# ----------------------------#

table_model1 <- plot_global_methylation %>%
  filter(FDR_significance == "Yes") %>%
  group_by(DNAm, brain_metric) %>%
  arrange(brain_metric,
          estimate,
          by_group = TRUE) %>%
  mutate(CI_lower =
           estimate - (1.96 * std.error),
         CI_upper =
           estimate + (1.96 * std.error)) %>%
  select(brain_metric, DNAm, estimate, CI_lower, CI_upper, r2, p.value, pFDR)

# ----------------------------#
# Examine which DNAm associate across multiple cortical regions
# ----------------------------#


### Which DNAm proxy has the most FDR significant hits?
sigs <- plot_global_methylation %>%
  #filter(estimate < 0 & FDR_significance == "Yes") %>%
  filter(
    brain_metric == "global grey matter" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "global white matter" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "total brain volume" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "global cortical thickness" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "global cortical volume" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "global cortical surface area" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "gFA" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "global subcortical volume" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "WMH" & estimate > 0 & FDR_significance == "Yes" |
      brain_metric == "gMD" & estimate > 0 & FDR_significance == "Yes" |
      brain_metric == "g" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "gf" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "relative brain age" & estimate > 0 & FDR_significance == "Yes" |
      brain_metric == "APOE" & estimate > 0 & FDR_significance == "Yes" |
      brain_metric == "processing speed" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "executive function" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "vocabulary" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "verbal declarative memory" & estimate < 0 & FDR_significance == "Yes" |
      brain_metric == "matrix reasoning" & estimate < 0 & FDR_significance == "Yes"
  ) %>%
  group_by(DNAm, brain_metric) %>%
  mutate(number_significant =
           case_when(FDR_significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_sigs <-
  aggregate(sigs$number_significant,
            by = list(DNAm = sigs$DNAm),
            FUN = sum) %>%
  arrange(desc(x)) %>%
  rename(n_pFDR = x)
plot2 <- plot_global_methylation


## Separate into direction of effect where increases in DNAm associate with poor brain health outcomes
## Excepting gMD and WMH which will show increases in effect size corresponding to poor brain health metrics...

Good_brain_health <- plot2 %>% filter(brain_metric == "global grey matter" & estimate > 0 |
                                        brain_metric == "global white matter" & estimate > 0 |
                                        brain_metric == "total brain volume" & estimate > 0 |
                                        brain_metric == "global cortical thickness" & estimate > 0 |
                                        brain_metric == "global cortical volume" & estimate > 0 |
                                        brain_metric == "global cortical surface area" & estimate > 0 |
                                        brain_metric == "gFA" & estimate > 0 |
                                        brain_metric == "global subcortical volume" & estimate > 0 |
                                        brain_metric == "WMH" & estimate > 0 |
                                        brain_metric == "gMD" & estimate > 0 |
                                        brain_metric == "g" & estimate > 0 |
                                        brain_metric == "gf" & estimate > 0 |
                                        brain_metric == "relative brain age" & estimate < 0 |
                                        brain_metric == "APOE" & estimate < 0 |
                                        
                                        brain_metric =="processing speed" & estimate > 0 |
                                        brain_metric =="executive function" & estimate > 0 |
                                        brain_metric =="vocabulary" & estimate > 0|
                                        brain_metric =="verbal declarative memory" & estimate > 0|
                                        brain_metric =="matrix reasoning" & estimate > 0
                                      
                                      
)

### Which DNAm proxy has the most FDR significant hits?
Good_brain_health %<>%
  mutate(number_significant =
           case_when(significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Good_brain_health <-
  aggregate(
    Good_brain_health$number_significant,
    by = list(brain_metric = Good_brain_health$brain_metric),
    FUN = sum
  )

Top_hits_Good_brain_health %<>% arrange(desc(x)) %>% rename(p = x)

# Run FDR
Good_brain_health %<>%
  mutate(number_significant_FDR =
           case_when(FDR_significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Good_brain_health_2 <-
  aggregate(
    Good_brain_health$number_significant_FDR,
    by = list(brain_metric = Good_brain_health$brain_metric),
    FUN = sum
  ) %>% rename(pFDR = x)


Good_brain <- merge(Top_hits_Good_brain_health,
                    Top_hits_Good_brain_health_2,
                    by = "brain_metric") %>%
  mutate(direction =
           "Favourable brain health outcome")

#### Poor brain health
Poor_brain_health <- plot2 %>% filter(brain_metric == "global grey matter" & estimate < 0 |
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
                                        brain_metric == "APOE" & estimate > 0 |
                                        
                                        brain_metric =="processing speed" & estimate < 0 |
                                        brain_metric =="executive function" & estimate < 0 |
                                        brain_metric =="vocabulary" & estimate < 0|
                                        brain_metric =="verbal declarative memory" & estimate < 0|
                                        brain_metric =="matrix reasoning" & estimate < 0)

### Which DNAm proxy has the most FDR significant hits?
Poor_brain_health %<>%
  mutate(number_significant =
           case_when(significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Poor_brain_health <-
  aggregate(
    Poor_brain_health$number_significant,
    by = list(brain_metric = Poor_brain_health$brain_metric),
    FUN = sum
  )

Top_hits_Poor_brain_health %<>% arrange(desc(x)) %>% rename(p = x)

# Run FDR
Poor_brain_health %<>%
  mutate(number_significant_FDR =
           case_when(FDR_significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Poor_brain_health_2 <-
  aggregate(
    Poor_brain_health$number_significant_FDR,
    by = list(brain_metric = Poor_brain_health$brain_metric),
    FUN = sum
  ) %>% rename(pFDR = x)


Poor_brain <- merge(Top_hits_Poor_brain_health,
                    Top_hits_Poor_brain_health_2,
                    by = "brain_metric") %>%
  mutate(direction =
           "Poor brain health outcome")

Significance_plot <- rbind(Poor_brain, Good_brain)

# ----------------------------#
# Code for Barplot of count of significant hits
# ----------------------------#

####
ggplot(Significance_plot,
       aes(
         y = reorder(brain_metric, p),
         x = p,
         fill = reorder(brain_metric, p)
       )) +
  
  geom_col(aes(y = reorder(brain_metric, p),
               x = p)) +
  
  geom_col(aes(
    y = reorder(brain_metric, p),
    x = pFDR,
    col = "grey"
  )) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",
    
    strip.text = element_text(
      size = 6,
      face = "bold",
      family = "sans",
      colour = "black"
    ),
    axis.text.x = element_text(
      vjust = 0.5,
      hjust = 1,
      size = 6,
      family = "sans"
    ),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(
      size = 7,
      face = "bold",
      family = "sans"
    ),
    axis.title.x = element_text(
      size = 7,
      face = "bold",
      family = "sans"
    )
  ) +
  
  labs(x = "number of significant assocations with DNAm signatures",
       y = "global brain metrics") +
  
  facet_wrap( ~ direction) +
  scale_fill_manual(values = c(
    # viridis::magma(n = 19)
    colorspace::sequential_hcl(19, palette = "SunsetDark")
  )) +
  scale_colour_manual(values = c("#808080"))
