

#Load methylation regressions created previously
#write.csv(plot_global_methylation, "plot_global_methylation.csv")
#write.csv(plot_global_methylation_lifestyle, "plot_global_methylation_lifestyle.csv")
# can be accessed at: https://github.com/EleanorSC/Inflammatory-DNAm_STRADL/tree/main/results/gglobal_metrics

plot_global_methylation <- read.csv("plot_global_methylation.csv")
plot_global_methylation_lifestyle <- read.csv("plot_global_methylation_lifestyle.csv")


# ----------------------------#
# PROTEOMICS
# ----------------------------#

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
  "ACY1_proteomics"    ,
  "ADAMTS13_proteomics"  ,
  "ADIPOQ_proteomics",
  "AFM_proteomics"       ,
  "B2M_proteomics"      ,
  "BCAM_proteomics"     ,
  "BMP1_proteomics"   ,
  "C4A.C4B_proteomics" ,
  "C5_proteomics"      ,
  "C9_proteomics"       ,
  "CCL17_proteomics"     ,
  "CCL18_proteomics"   ,
  "CCL21_proteomics"    ,
  "CCL22_proteomics"  ,
  "CCL25_proteomics"  ,
  "CD163_proteomics" ,
  "CD209_proteomics"   ,
  "CD48_proteomics"      ,
  "CD5L_proteomics"      ,
  "CHIT1_proteomics"     ,
  "CLEC11A_proteomics"  ,
  "CLEC11A.1_proteomics",
  "CNTN4_proteomics"  ,
  "CRP_proteomics"     ,
  "CXCL10_proteomics"  ,
  "CXCL11_proteomics"    ,
  "EDA_proteomics"       ,
  "ENPP7_proteomics"     ,
  "ESM1_proteomics"     ,
  "F7_proteomics"       ,
  "FAP_proteomics"    ,
  "FCER2_proteomics"   ,
  "FCGR3B_proteomics"  ,
  "GHR_proteomics"       ,
  "GZMA_proteomics" ,
  "GNLY_proteomics"      ,
  "GP1BA_proteomics"     ,
  "HGFAC_proteomics"    ,
  "ICAM5_proteomics"  ,
  "IDUA_proteomics"    ,
  "IGFBP1_proteomics"  ,
  "IGFBP4_proteomics"    ,
  "IL19_proteomics"      ,
  "INSR_proteomics"      ,
  "LGALS3BP_proteomics" ,
  "LGALS4_proteomics"   ,
  "LTA.LTB_proteomics",
  "LTF_proteomics"   ,
  "LY9_proteomics"     ,
  "LYZ_proteomics"     ,
  "MIA_proteomics"       ,
  "MMP1_proteomics"     ,
  "MMP12_proteomics"    ,
  "MMP2_proteomics"   ,
  "MMP9_proteomics"    ,
  "MPL_proteomics"     ,
  "MPO_proteomics"       ,
  "MRC2_proteomics"      ,
  "MST1_proteomics"      ,
  "NCAM1_proteomics"    ,
  "NOTCH1_proteomics"   ,
  "NTRK3_proteomics"  ,
  "OMD_proteomics"     ,
  "PAPPA_proteomics"   ,
  "PIGR_proteomics"      ,
  "PRSS2_proteomics"     ,
  "RARRES2_proteomics"   ,
  "RETN_proteomics"     ,
  "S100A9_proteomics"   ,
  "SELE_proteomics"   ,
  "SELL_proteomics"   ,
  "SEMA3E_proteomics"  ,
  "SERPINA3_proteomics"  ,
  "SERPIND1_proteomics" ,
  "SHBG_proteomics"      ,
  "SLITRK5_proteomics"  ,
  "SPOCK2_proteomics"   ,
  "STC1_proteomics"   ,
  "THBS2_proteomics"   ,
  "TNFRSF17_proteomics"  ,
  "TNFRSF1B_proteomics" ,
  "TPSB2_proteomics"     ,
  "VCAM1_proteomics"    ,
  "WFIKKN2_proteomics",
  
  # LBC trained
  
  "CRTAM_proteomics"    ,
  "EZR_proteomics"      ,
  "NMNAT1_proteomics"   ,
  "SMPD1_proteomics"    ,
  "CCL11_proteomics"  ,
  "CXCL9_proteomics"    ,
  "HGF_proteomics"      ,
  "OSM_proteomics"      ,
  "VEGFA_proteomics" )

#  FcRL2    , #NO PROTEIN EQUIVALENT
#  G.CSF    , #NO PROTEIN EQUIVALENT
#  GDF.8    , #NO PROTEIN EQUIVALENT
#  N.CDase  , #NO PROTEIN EQUIVALENT
#  NEP      , #NO PROTEIN EQUIVALENT
#  SIGLEC1  , #NO PROTEIN EQUIVALENT
#  SKR3     , #NO PROTEIN EQUIVALENT
#  CD6      , #NO PROTEIN EQUIVALENT
#  EN.RAGE  , #NO PROTEIN EQUIVALENT
#  FGF.21   , #NO PROTEIN EQUIVALENT
#  TGF.alpha, #NO PROTEIN EQUIVALENT



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
  
  # rename DNAm column to remove _proteomics
  mutate(DNAm =  gsub('_proteomics', '', DNAm)) %>%
  
  # Create a column to denote whether these are DNAm or proteins
  mutate(omic_type = "protein") %>%
  mutate(model = "Model 1 (n=709)")

###
plot_global_neuroimaging_protein <- newdf

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
  
  # rename DNAm column to remove _proteomics
  mutate(DNAm =  gsub('_proteomics', '', DNAm)) %>%
  
  # Create a column to denote whether these are DNAm or proteins
  mutate(omic_type = "protein") %>%
  mutate(model = "Model 1 (n=709)")

###
plot_cognitive_protein <- newdf


# ----------------------------#
# Combine cognitive and brain analyses
# ----------------------------#

plot_global_protein <- rbind(plot_global_neuroimaging_protein,
                                 plot_cognitive_protein)


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

Neuroimaging_DNAm_lifestyle <- merge(Neuroimaging_DNAm_lifestyle,
                                     Cognitive_data_select,
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
  
  # rename DNAm column to remove _proteomics
  mutate(DNAm =  gsub('_proteomics', '', DNAm)) %>%
  
  # Create a column to denote whether these are DNAm or proteins
  mutate(omic_type = "protein") %>%
  mutate(model = "Model 3 Lifestyle (n=709)")

###
plot_global_neuroimaging_protein_lifestyle <- newdf



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
      
      + hypertension
      + CurrentSmoker
      + CurrentDrinker
      + scale(bmi)
      
      + scale(Neuroimaging_DNAm_lifestyle[[DNAm]]),
      data = Neuroimaging_DNAm_lifestyle
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
  
  # rename DNAm column to remove _proteomics
  mutate(DNAm =  gsub('_proteomics', '', DNAm)) %>%
  
  # Create a column to denote whether these are DNAm or proteins
  mutate(omic_type = "protein") %>%
  mutate(model = "Model 3 Lifestyle (n=709)")

###
plot_cognitive_protein_lifestyle <- newdf

# ----------------------------#
# Combine cognitive and brain analyses
# ----------------------------#

plot_global_protein_lifestyle <- rbind(plot_global_neuroimaging_protein_lifestyle,
                                           plot_cognitive_protein_lifestyle)


## ----------------------------#
# Combining both sets of models to see percentage attenuation
## ----------------------------#

add_on <- plot_global_protein_lifestyle %>%
  select(estimate, std.error, r2, p.value, pFDR) %>%
  rename(
    estimate_lifestyle = estimate,
    r2_lifestyle = r2,
    p.value_lifestyle = p.value,
    pFDR_lifestyle = pFDR,
    std.error_lifestyle = std.error
  )

add_on <- subset(add_on, select = -c(brain_metric))

table_new <- cbind(plot_global_protein, add_on)

table_new %<>% mutate(percentage_increase_decrease =
                        100 * (estimate - estimate_lifestyle)) %>%
  
  # filter(FDR_significance == "Yes" & estimate < 0) %>%
  # filter(significance == "Yes") %>%
  
  group_by(DNAm, brain_metric) %>%
  
  
  arrange(brain_metric,
          estimate,
          by_group = TRUE) %>%
  
  mutate(
    CI_lower =
      estimate - (1.96 * std.error),
    CI_upper =
      estimate + (1.96 * std.error),
    
    CI_lower_lifestyle =
      estimate_lifestyle - (1.96 * std.error_lifestyle),
    CI_upper_lifestyle =
      estimate_lifestyle + (1.96 * std.error_lifestyle),
    
    better_model =
      case_when(r2_lifestyle > r2 ~ "Yes",
                TRUE ~ "No")) %>%
  select(
    brain_metric,
    DNAm,
    # Model 1
    estimate,
    CI_lower,
    CI_upper,
    p.value,
    pFDR,
    # Model 2 
    estimate_lifestyle,
    CI_lower_lifestyle,
    CI_upper_lifestyle,
    p.value_lifestyle,
    pFDR_lifestyle,
    # R2
    r2,
    r2_lifestyle,
    better_model,
    percentage_increase_decrease
  ) 


# ----------------------------#
# write to .csv
# ----------------------------#

Ultra <- rbind(
  plot_global_protein,
  plot_global_protein_lifestyle,
  plot_global_methylation,
  plot_global_methylation_lifestyle
)

#################################

baseline_model <- rbind(plot_global_protein,
                        plot_global_methylation) %>%
  filter(significance == "Yes")


Good_brain_health_protein <- baseline_model %>% filter(
                                                 omic_type == "protein" & brain_metric == "global grey matter" & estimate > 0 |
                                                 omic_type == "protein" & brain_metric == "global white matter" & estimate > 0 |
                                                 omic_type == "protein" & brain_metric == "total brain volume" & estimate > 0 |
                                                 omic_type == "protein" & brain_metric == "global cortical thickness" & estimate > 0 |
                                                 omic_type == "protein" & brain_metric == "global cortical volume" & estimate > 0 |
                                                 omic_type == "protein" & brain_metric == "global cortical surface area" & estimate > 0 |
                                                 omic_type == "protein" & brain_metric == "gFA" & estimate > 0 |
                                                 omic_type == "protein" & brain_metric == "global subcortical volume" & estimate > 0 |
                                                 omic_type == "protein" & brain_metric == "WMH" & estimate > 0 |
                                                 omic_type == "protein" & brain_metric == "gMD" & estimate > 0 |
                                                 omic_type == "protein" & brain_metric == "g" & estimate > 0 |
                                                 omic_type == "protein" & brain_metric == "gf" & estimate > 0 |
                                                 omic_type == "protein" & brain_metric == "relative brain age" & estimate < 0 |
                                                 omic_type == "protein" & brain_metric == "APOE" & estimate < 0 |
                                                
                                                 omic_type == "protein" & brain_metric =="processing speed" & estimate > 0 |
                                                 omic_type == "protein" & brain_metric =="executive function" & estimate > 0 |
                                                 omic_type == "protein" & brain_metric =="vocabulary" & estimate > 0|
                                                 omic_type == "protein" & brain_metric =="verbal declarative memory" & estimate > 0|
                                                 omic_type == "protein" & brain_metric =="matrix reasoning" & estimate > 0
                                               
                                              )

### Which DNAm proxy has the most FDR significant hits?

Good_brain_health_protein %<>%
  mutate(number_significant =
           case_when(significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Good_brain_health_protein <-
  aggregate(
    Good_brain_health_protein$number_significant,
    by = list(DNAm = Good_brain_health_protein$DNAm),
    FUN = sum) %>%
  arrange(desc(x)) %>%
  rename(n_p = x)


# Run FDR
Good_brain_health_protein %<>%
  mutate(number_significant_FDR =
           case_when(FDR_significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Good_brain_health_protein_2 <-
  aggregate(
    Good_brain_health_protein$number_significant_FDR,
    by = list(DNAm = Good_brain_health_protein$DNAm),
    FUN = sum) %>%
  arrange(desc(x)) %>%
  rename(n_pFDR = x)



Good_brain_protein <- merge(Top_hits_Good_brain_health_protein,
                    Top_hits_Good_brain_health_protein_2,
                    by = "DNAm") %>%
  mutate(direction =
           "Favourable brain health outcome",
         omic_type = "protein")

#### Poor brain health
Poor_brain_health_protein <- baseline_model %>% filter(  omic_type == "protein" & brain_metric == "global grey matter" & estimate < 0 |
                                                 omic_type == "protein" & brain_metric == "global white matter" & estimate < 0 |
                                                 omic_type == "protein" & brain_metric == "total brain volume" & estimate < 0 |
                                                 omic_type == "protein" & brain_metric == "global cortical thickness" & estimate < 0 |
                                                 omic_type == "protein" & brain_metric == "global cortical volume" & estimate < 0 |
                                                 omic_type == "protein" & brain_metric == "global cortical surface area" & estimate < 0 |
                                                 omic_type == "protein" & brain_metric == "gFA" & estimate < 0 |
                                                 omic_type == "protein" & brain_metric == "global subcortical volume" & estimate < 0 |
                                                 omic_type == "protein" & brain_metric == "WMH" & estimate > 0 |
                                                 omic_type == "protein" & brain_metric == "gMD" & estimate > 0 |
                                                 omic_type == "protein" & brain_metric == "g" & estimate < 0 |
                                                 omic_type == "protein" & brain_metric == "gf" & estimate < 0 |
                                                 omic_type == "protein" & brain_metric == "relative brain age" & estimate > 0 |
                                                 omic_type == "protein" & brain_metric == "APOE" & estimate > 0 |
                                                 
                                                 omic_type == "protein" & brain_metric =="processing speed" & estimate < 0 |
                                                 omic_type == "protein" & brain_metric =="executive function" & estimate < 0 |
                                                 omic_type == "protein" & brain_metric =="vocabulary" & estimate < 0|
                                                 omic_type == "protein" & brain_metric =="verbal declarative memory" & estimate < 0|
                                                 omic_type == "protein" & brain_metric =="matrix reasoning" & estimate < 0)

### Which DNAm proxy has the most FDR significant hits?
### Which DNAm proxy has the most FDR significant hits?
Poor_brain_health_protein %<>%
  mutate(number_significant =
           case_when(significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Poor_brain_health_protein <-
  aggregate(
    Poor_brain_health_protein$number_significant,
    by = list(DNAm = Poor_brain_health_protein$DNAm),
    FUN = sum) %>%
  arrange(desc(x)) %>%
  rename(n_p = x)


# Run FDR
Poor_brain_health_protein %<>%
  mutate(number_significant_FDR =
           case_when(FDR_significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Poor_brain_health_protein_2 <-
  aggregate(
    Poor_brain_health_protein$number_significant_FDR,
    by = list(DNAm = Poor_brain_health_protein$DNAm),
    FUN = sum) %>%
  arrange(desc(x)) %>%
  rename(n_pFDR = x)



Poor_brain_protein <- merge(Top_hits_Poor_brain_health_protein,
                    Top_hits_Poor_brain_health_protein_2,
                    by = "DNAm") %>%
  mutate(direction =
           "Poor brain health outcome",
         omic_type = "protein")

#####################################
Good_brain_health_DNAm <- baseline_model %>% filter(
  omic_type == "DNAm" & brain_metric == "global grey matter" & estimate > 0 |
    omic_type == "DNAm" & brain_metric == "global white matter" & estimate > 0 |
    omic_type == "DNAm" & brain_metric == "total brain volume" & estimate > 0 |
    omic_type == "DNAm" & brain_metric == "global cortical thickness" & estimate > 0 |
    omic_type == "DNAm" & brain_metric == "global cortical volume" & estimate > 0 |
    omic_type == "DNAm" & brain_metric == "global cortical surface area" & estimate > 0 |
    omic_type == "DNAm" & brain_metric == "gFA" & estimate > 0 |
    omic_type == "DNAm" & brain_metric == "global subcortical volume" & estimate > 0 |
    omic_type == "DNAm" & brain_metric == "WMH" & estimate > 0 |
    omic_type == "DNAm" & brain_metric == "gMD" & estimate > 0 |
    omic_type == "DNAm" & brain_metric == "g" & estimate > 0 |
    omic_type == "DNAm" & brain_metric == "gf" & estimate > 0 |
    omic_type == "DNAm" & brain_metric == "relative brain age" & estimate < 0 |
    omic_type == "DNAm" & brain_metric == "APOE" & estimate < 0 |
    
    omic_type == "DNAm" & brain_metric =="processing speed" & estimate > 0 |
    omic_type == "DNAm" & brain_metric =="executive function" & estimate > 0 |
    omic_type == "DNAm" & brain_metric =="vocabulary" & estimate > 0|
    omic_type == "DNAm" & brain_metric =="verbal declarative memory" & estimate > 0|
    omic_type == "DNAm" & brain_metric =="matrix reasoning" & estimate > 0
  
)

### Which DNAm proxy has the most FDR significant hits?

Good_brain_health_DNAm %<>%
  mutate(number_significant =
           case_when(significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Good_brain_health_DNAm <-
  aggregate(
    Good_brain_health_DNAm$number_significant,
    by = list(DNAm = Good_brain_health_DNAm$DNAm),
    FUN = sum) %>%
  arrange(desc(x)) %>%
  rename(n_p = x)


# Run FDR
Good_brain_health_DNAm %<>%
  mutate(number_significant_FDR =
           case_when(FDR_significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Good_brain_health_DNAm_2 <-
  aggregate(
    Good_brain_health_DNAm$number_significant_FDR,
    by = list(DNAm = Good_brain_health_DNAm$DNAm),
    FUN = sum) %>%
  arrange(desc(x)) %>%
  rename(n_pFDR = x)



Good_brain_DNAm <- merge(Top_hits_Good_brain_health_DNAm,
                         Top_hits_Good_brain_health_DNAm_2,
                         by = "DNAm") %>%
  mutate(direction =
           "Favourable brain health outcome",
         omic_type = "DNAm")

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
### Which DNAm proxy has the most FDR significant hits?
Poor_brain_health_DNAm %<>%
  mutate(number_significant =
           case_when(significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Poor_brain_health_DNAm <-
  aggregate(
    Poor_brain_health_DNAm$number_significant,
    by = list(DNAm = Poor_brain_health_DNAm$DNAm),
    FUN = sum) %>%
  arrange(desc(x)) %>%
  rename(n_p = x)


# Run FDR
Poor_brain_health_DNAm %<>%
  mutate(number_significant_FDR =
           case_when(FDR_significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Poor_brain_health_DNAm_2 <-
  aggregate(
    Poor_brain_health_DNAm$number_significant_FDR,
    by = list(DNAm = Poor_brain_health_DNAm$DNAm),
    FUN = sum) %>%
  arrange(desc(x)) %>%
  rename(n_pFDR = x)



Poor_brain_DNAm <- merge(Top_hits_Poor_brain_health_DNAm,
                         Top_hits_Poor_brain_health_DNAm_2,
                         by = "DNAm") %>%
  mutate(direction =
           "Poor brain health outcome",
         omic_type = "DNAm")

#####################################
Significance_plot <- rbind(Poor_brain_protein, Good_brain_protein,
                           Poor_brain_DNAm, Good_brain_DNAm)

# ----------------------------#
# Code for Barplot of count of significant hits
# ----------------------------#

Significance_plot %>% group_by(direction)

Significance_plot_poor <- Significance_plot %>% filter(direction == "Poor brain health outcome")

Significance_plot_poor_2 <- Significance_plot_poor %>% filter(n_p > 0 & omic_type == "DNAm")

#test <- list(Significance_plot_poor_2$DNAm)
# ----------------------------#
# Only plot instances where DNAm have significant hits to compare with protein equivalents
# ----------------------------#

Significance_plot_poor_3 <- Significance_plot_poor %>% filter( DNAm == "ACY1"       |DNAm ==  "AFM"    | DNAm==      "B2M"          | DNAm== "BCAM"      | DNAm ==   "BMP1"         | DNAm =="C4A.C4B" |  DNAm ==   "CCL11"     |  
                                                               DNAm == "CCL17"      |DNAm ==  "CCL18"  | DNAm==      "CCL22"        | DNAm== "CCL25"     | DNAm ==   "CD5L"         | DNAm =="CHIT1"   |  DNAm ==   "CLEC11A.1" |  
                                                               DNAm == "CRP"        |DNAm ==  "CXCL10" | DNAm==      "CXCL10_olink" | DNAm== "CXCL11"    | DNAm ==   "CXCL9"        | DNAm =="EN.RAGE" |  DNAm ==   "ENPP7"     |  
                                                               DNAm == "FCGR3B"     |DNAm ==  "FGF.21" | DNAm==      "G.CSF"        | DNAm== "GHR"       | DNAm ==   "HGF"          | DNAm =="HGFAC"   |  DNAm ==   "ICAM5"     |  
                                                               DNAm == "IDUA"       |DNAm ==  "IGFBP4" | DNAm==      "LGALS3BP"     | DNAm== "LYZ"       | DNAm ==   "MMP.1_olink"  | DNAm =="MMP1"    |  DNAm ==   "MMP12"     |  
                                                               DNAm == "MMP9"       |DNAm ==  "MST1"   | DNAm==      "NEP"          | DNAm== "OSM"       | DNAm ==   "PIGR"         | DNAm =="PRSS2"   |  DNAm ==   "RARRES2"   |  
                                                               DNAm == "RETN"       |DNAm ==  "S100A9" | DNAm==      "SELE"         | DNAm== "SERPIND1"  | DNAm ==   "SIGLEC1"      | DNAm =="SKR3"    |  DNAm ==   "STC1"      |  
                                                               DNAm == "TGF.alpha"  |DNAm ==  "THBS2"  | DNAm==      "TNFRSF17"     | DNAm== "TPSB2"     | DNAm ==   "VEGFA")

####
ggplot(Significance_plot_poor_3,
       aes(
         y = reorder(DNAm, n_p),
         x = n_p
         
       )) +
  
  geom_col(aes(y = reorder(DNAm, n_p),
               x = n_p,
               fill = reorder(DNAm, n_p))) +
  
  geom_col(aes(
    y = reorder(DNAm, n_p),
    x = n_pFDR,
    colour = c("grey"),
    alpha = 0.95
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
  
  labs(x = "number of significant assocations with brain and cognitive phenotypes",
       y = "") +
  
  facet_wrap( ~ omic_type
              #scales = "free"
              ) +
  scale_fill_manual(values = c(
    # viridis::magma(n = 19)
    colorspace::sequential_hcl(74, palette = "SunsetDark")
  )) +
  scale_colour_manual(values = c(
    "#BEBEBE"
    # "#808080"
  ))


# ----------------------------#
# Examine instances where protein has more hits
# ----------------------------#
DNAm_RARRES2 <- Poor_brain_health_DNAm %>% filter(DNAm == "RARRES2")
protein_RARRES2 <- Poor_brain_health_protein %>% filter(DNAm == "RARRES2")


Poor_brain_health_DNAm2 <- Poor_brain_health_DNAm%>% filter(pFDR < 0.05)
skimr::skim(abs(Poor_brain_health_DNAm2$estimate))

Poor_brain_health_DNAm2 <- Poor_brain_health_DNAm %>% filter(estimate < 0.053 & pFDR < 0.05)