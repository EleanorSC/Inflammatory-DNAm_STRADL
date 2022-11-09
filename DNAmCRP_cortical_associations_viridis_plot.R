global_neuroimaging_dataset <- merge(STRADL_lifestyle_covariates, 
                                     dataset_n655, 
                                     by = "stradl_ID") 


### Examine global measures, add L and R hemispheres

global_neuroimaging_dataset %<>% mutate(global_cortical_surface_area = 
                                          hem.lh.csa + hem.rh.csa) 

global_neuroimaging_dataset %<>% mutate(global_cortical_thickness = 
                                          hem.lh.ct + hem.rh.ct) 

global_neuroimaging_dataset %<>% mutate(global_cortical_volume = 
                                          hem.lh.cv + hem.rh.cv) 


#skimr::skim(global_neuroimaging_dataset)

FULL_neuroimaging_list <- c(
  ### Global measures
  "global.cerebral.wm",
  "global.total.gm",
  "global.wbv",
  # "global.wbv.vent",
  "global_cortical_volume",
  "global_cortical_thickness",
  "global_cortical_surface_area"
  
)

FULL_DNAm_list <- c("CRP", 
                    "IGFBP4",
                    "PIGR",
                    "MMP12")

# Converting multiple varibles into a factor
global_neuroimaging_dataset %<>% mutate_at(c("sex", "site", "edited", "batch"),
                                           as.factor)

# Making a data frame with all combinations of variables
df <- as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1, metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
get_summary_values <- function(list) {
  metric <- list$metric
  DNAm <- list$DNAm
  tib <-
    broom::tidy(summary(lm(
      scale(global_neuroimaging_dataset[[metric]]) ~ 
        scale(st_age)
      + sex
      + site
      + batch
      + edited
      + scale(Standardised_ICV)
      + scale(global_neuroimaging_dataset[[DNAm]]),
      data = global_neuroimaging_dataset
    )))[8, c(2, 3, 5)]
  return(tib)
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = get_summary_values)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality =
           case_when(
             str_detect(metric, "global") ~ "global",
             TRUE ~ "misc")
  ) %>%
  
  # Create a hemisphere column
  mutate(Hemisphere =
           case_when(
             str_detect(metric, "\\.l") ~ "Left",
             str_detect(metric, "\\.h") ~ "Right",
             TRUE ~ "Bilateral"
           )) %>%
  
  
  # Create a brain metric column
  mutate(brain_metric =
           case_when(
             
             ### Global
             metric == "global.cerebral.wm" ~ "global white matter",
             metric == "global.total.gm" ~ "global grey matter",
             metric == "global.wbv" ~ "total brain volume",
             #metric == "global.wbv.vent" ~ "global ventricles volume",
             metric == "global_cortical_volume" ~ "global cortical volume",
             metric == "global_cortical_thickness" ~ "global cortical thickness",
             metric == "global_cortical_surface_area" ~ "global cortical surface area",
             
             TRUE ~ "misc")
  ) %>%
  
  # Create pFDR column
  # group_by(metric, Hemisphere) %>%
  group_by(DNAm) %>%
  mutate(pFDR = p.adjust(p.value, method = "fdr")) %>%
  # Create a column to denote where there are significant hits 
  mutate(significance =
           case_when(p.value < 0.05 ~ "Yes",
                     TRUE ~ "No")) %>%
  # Create a column to denote where there are FDR significant hits 
  mutate(FDR_significance =
           case_when(pFDR < 0.05 ~ "Yes",
                     TRUE ~ "No"))%>%
  # Create a column to denote where there are FDR significant hits, significant, or not
  mutate(SIG =
           case_when(FDR_significance == "Yes" ~ "pFDR",
                     significance == "Yes" ~ "p",
                     TRUE ~ "not_significant"))

###
#plot %<>% filter(DNAm == "CRP")

plot_global <- newdf


#### SUBCORTICVAL
global_neuroimaging_dataset %<>% mutate(scv.accumbens =  scv.lh.accumbens + scv.rh.accumbens,
                                        scv.amygdala =  scv.lh.amygdala + scv.rh.amygdala,
                                        scv.caudate =  scv.lh.caudate + scv.rh.caudate,
                                        scv.hippocampus =  scv.lh.hippocampus + scv.rh.hippocampus,
                                        scv.pallidum =  scv.lh.pallidum + scv.rh.pallidum,
                                        scv.putamen =  scv.lh.putamen + scv.rh.putamen,
                                        scv.thalamus =  scv.lh.thalamus + scv.rh.thalamus,
                                        scv.ventraldc =  scv.lh.ventraldc + scv.rh.amygdala,
                                        scv.cerebellum_GM =  vol.cerebellum.lh.gm + vol.cerebellum.rh.gm,
                                        scv.cerebellum_WM =  vol.cerebellum.lh.wm + vol.cerebellum.rh.wm)


FULL_neuroimaging_list <- c(
  
  ### subcortical volumes
  "scv.accumbens",
  "scv.amygdala",
  "scv.caudate",
  "scv.hippocampus",
  "scv.pallidum",
  "scv.putamen",
  "scv.thalamus",
  "scv.ventraldc",
  "vol.brainstem",
  "scv.cerebellum_GM",
  "scv.cerebellum_WM"
)
# Making a data frame with all combinations of variables
df <- as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1, metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
# e.g. try
# info <- list(metric = "ACR", DNAm = "ACY1")
# get_summary_values(info)
# It will return you estimate, std.error and p.value in a tibble
get_summary_values <- function(list) {
  metric <- list$metric
  DNAm <- list$DNAm
  tib <-
    broom::tidy(summary(lm(
      scale(global_neuroimaging_dataset[[metric]]) ~ 
        scale(st_age)
      + sex
      + site
      + batch
      + edited
      + scale(Standardised_ICV)
      + scale(global_neuroimaging_dataset[[DNAm]]),
      data = global_neuroimaging_dataset
    )))[8, c(2, 3, 5)]
  return(tib)
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = get_summary_values)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality =
           case_when(
             str_detect(metric, "scv.") ~ "subcortical volumes",
             str_detect(metric, "vol") ~ "subcortical volumes",
             TRUE ~ "misc")
  ) %>%
  
  # Create a hemisphere column
  mutate(Hemisphere =
           case_when(
             str_detect(metric, "\\.lh.") ~ "Left",
             str_detect(metric, "\\.rh.") ~ "Right",
             TRUE ~ "Bilateral"
           )) %>%
  
  # Create a brain metric column
  mutate(brain_metric =
           case_when(
             
             ### Global
             metric == "scv.accumbens" ~ "accumbens",
             metric == "scv.amygdala" ~ "amygdala",
             metric == "scv.caudate"~ "caudate",
             metric == "scv.hippocampus"~ "hippocampus",
             metric == "scv.pallidum"~ "pallidum",
             metric == "scv.putamen"~ "putamen",
             metric == "scv.thalamus"~ "thalamus",
             metric == "scv.ventraldc"~ "ventral diencephalon",
             metric == "vol.brainstem"~ "brainstem",
             metric == "scv.cerebellum_GM"~ "cerebellum (GM)",
             metric == "scv.cerebellum_WM"~ "cerebellum (WM)",
             
             TRUE ~ "misc")
  ) %>%
  
  
  
  # Create pFDR column
  # group_by(metric, Hemisphere) %>%
  group_by(DNAm, metric) %>%
  mutate(pFDR = p.adjust(p.value, method = "fdr")) %>%
  # Create a column to denote where there are significant hits 
  mutate(significance =
           case_when(p.value < 0.05 ~ "Yes",
                     TRUE ~ "No")) %>%
  # Create a column to denote where there are FDR significant hits 
  mutate(FDR_significance =
           case_when(pFDR < 0.05 ~ "Yes",
                     TRUE ~ "No"))%>%
  # Create a column to denote where there are FDR significant hits, significant, or not
  mutate(SIG =
           case_when(FDR_significance == "Yes" ~ "pFDR",
                     significance == "Yes" ~ "p",
                     TRUE ~ "not_significant"))

###
plot_subcortical <- newdf

####### dMRI
FULL_neuroimaging_list <- c("ACR",
                            "ACR.L",
                            "ACR.R",
                            "ALIC",
                            "ALIC.L",
                            "ALIC.R",
                            "BCC",
                            "CC",
                            "CGC",
                            "CGC.L",
                            "CGC.R",
                            "CGH",
                            "CGH.L",
                            "CGH.R",
                            "CR",
                            "CR.L",
                            "CR.R",
                            "CST",
                            "CST.L",
                            "CST.R",
                            "EC",
                            "EC.L",
                            "EC.R",
                            "FX",
                            "FX.ST.L",
                            "FX.ST.R",
                            "FX.ST",
                            "GCC",
                            "IC",
                            "IC.L",
                            "IC.R",
                            "IFO",
                            "IFO.L",
                            "IFO.R",
                            "PCR",
                            "PCR.L",
                            "PCR.R",
                            "PLIC",
                            "PLIC.L",
                            "PLIC.R",
                            "PTR",
                            "PTR.L",
                            "PTR.R",
                            "RLIC",
                            "RLIC.L",
                            "RLIC.R",
                            "SCC",
                            "SCR",
                            "SCR.L",
                            "SCR.R",
                            "SFO",
                            "SFO.L",
                            "SFO.R",
                            "SLF",
                            "SLF.L",
                            "SLF.R",
                            "SS",
                            "SS.L",
                            "SS.R",
                            "UNC",
                            "UNC.L",
                            "UNC.R",
                            "gFA",
                            "gMD")

DTI_2 <-merge(STRADL_lifestyle_covariates,STRADL_DTI, by = "stradl_ID")
DTI_dataset <- merge(KORA_LBC_DNAm, DTI_2, by = "stradl_ID")
# n = 683

# Invert polarity
DTI_dataset  %<>% mutate(gFA = -1*gFA,
                         gMD = -1*gMD)

# Converting multiple varibles into a factor
DTI_dataset %<>% mutate_at(c("sex", "Site"),
                           as.factor)

#note change FXST -> FX.ST
names(DTI_dataset)[names(DTI_dataset) == "FXST"] <- "FX.ST"

#DTI_dataset %<>% drop_na(Site)

skimr::skim(DTI_dataset$st_age)

# Making a data frame with all combinations of variables
df <- as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1, metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
get_summary_values <- function(list) {
  metric <- list$metric
  DNAm <- list$DNAm
  tib <-
    broom::tidy(summary(lm(
      scale(DTI_dataset[[metric]]) ~ 
        scale(st_age)
      + sex
      + Site
      + scale(DTI_dataset[[DNAm]]),
      data = DTI_dataset
    )))[5, c(2, 3, 5)]
  return(tib)
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = get_summary_values)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality = 
           case_when(
             metric == "gFA" ~ "global",
             metric == "gMD" ~ "global",
             TRUE ~ "White matter integrity"
           )) %>%
  
  # Create a hemisphere column
  mutate(Hemisphere =
           case_when(
             str_detect(metric, "\\.L") ~ "Left",
             str_detect(metric, "\\.R") ~ "Right",
             TRUE ~ "Bilateral"
           )) %>%
  
  # Create a brain metric column
  mutate(
    brain_metric =
      case_when(
        ### Global
        metric == "gFA" ~ "gFA",
        metric == "gMD" ~ "gMD",
        
        metric == "ACR" ~ "ACR (FA)",
        #Anterior corona radiata
        metric == "ACR.L" ~ "ACR (FA)",
        #Anterior corona radiata
        metric == "ACR.R" ~ "ACR (FA)",
        #Anterior corona radiata
        metric == "ALIC" ~ "ALIC (FA)",
        #Anterior Limb of Internal Capsule
        metric == "ALIC.L" ~ "ALIC (FA)",
        #Anterior Limb of Internal Capsule
        metric == "ALIC.R" ~ "ALIC (FA)",
        #Anterior Limb of Internal Capsule
        metric == "CR" ~ "CR (FA)",
        #Corona Radiata
        metric == "CR.L" ~ "CR (FA)",
        #Corona Radiata
        metric == "CR.R" ~ "CR (FA)",
        #Corona Radiata
        metric == "IC" ~ "IC (FA)",
        #Internal Capsule
        metric == "IC.L" ~ "IC (FA)",
        #Internal Capsule
        metric == "IC.R" ~ "IC (FA)",
        #Internal Capsule
        metric == "RLIC" ~ "RLIC (FA)",
        #Retrolenticular part of internal capsule
        metric == "RLIC.L" ~ "RLIC (FA)",
        #Retrolenticular part of internal capsule
        metric == "RLIC.R" ~ "RLIC (FA)",
        #Retrolenticular part of internal capsule
        
        # str_detect(metric, "BCC") ~ "BCC (FA)",
        metric == "BCC" ~ "BCC (FA)",
        metric == "BCC.L" ~ "BCC (FA)",
        metric == "BCC.R" ~ "BCC (FA)",
        # #Body of corpus callosum
        metric == "CC" ~ "CC (FA)",
        metric == "CC.L" ~ "CC (FA)",
        metric == "CC.R" ~ "CC (FA)",
        #str_detect(metric, "CC") ~ "CC (FA)",
        # # Genu of corpus callosum
        metric == "CGC" ~ "CGC (FA)",
        metric == "CGC.L" ~ "CGC (FA)",
        metric == "CGC.R" ~ "CGC (FA)",
        #Cingulum-Cingulate Gyrus
        metric ==  "CGH" ~ "CGH (FA)",
        metric ==  "CGH.L" ~ "CGH (FA)",
        metric ==  "CGH.R" ~ "CGH (FA)",
        #Cingulum-Hippocampus
        
        metric == "CST" ~ "CST (FA)",
        metric == "CST.L" ~ "CST (FA)",
        metric == "CST.R" ~ "CST (FA)",
        #Corticospinal tract
        metric == "EC" ~ "EC (FA)",
        metric == "EC.L" ~ "EC (FA)",
        metric == "EC.R" ~ "EC (FA)",
        #External Capsule
        metric == "FX" ~ "FX (FA)",
        metric == "FX.L" ~ "FX (FA)",
        metric == "FX.R" ~ "FX (FA)",
        #Fornix-Column&Body
        metric == "FX.ST" ~ "FX.ST (FA)",
        metric == "FX.ST.L" ~ "FX.ST (FA)",
        metric == "FX.ST.R" ~ "FX.ST (FA)",
        #Fornix-Cres/Stria Terminalis #note I changed FXST -> FX.ST
        metric == "GCC" ~ "GCC (FA)",
        metric == "GCC.L" ~ "GCC (FA)",
        metric == "GCC.R" ~ "GCC (FA)",
        #Genu of corpus callosum
        
        metric == "IFO" ~ "IFO (FA)",
        metric == "IFO.L" ~ "IFO (FA)",
        metric == "IFO.R" ~ "IFO (FA)",
        #Inferior fronto-occipital fasciculus
        metric == "PCR" ~ "PCR (FA)",
        metric == "PCR.L" ~ "PCR (FA)",
        metric == "PCR.R" ~ "PCR (FA)",
        #Posterior corona radiata
        metric == "PLIC" ~ "PLIC (FA)",
        metric == "PLIC.L" ~ "PLIC (FA)",
        metric == "PLIC.R" ~ "PLIC (FA)",
        #Posterior limb of internal capsule
        metric == "PTR" ~ "PTR (FA)",
        metric == "PTR.L" ~ "PTR (FA)",
        metric == "PTR.R" ~ "PTR (FA)",
        #Posterior Thalamic Radiation
        
        metric == "SCC" ~ "SCC (FA)",
        metric == "SCC.L" ~ "SCC (FA)",
        metric == "SCC.R" ~ "SCC (FA)",
        # Splenium of corpus callosum
        metric == "SCR" ~ "SCR (FA)",
        metric == "SCR.L" ~ "SCR (FA)",
        metric == "SCR.R" ~ "SCR (FA)",
        #Superior corona radiata
        metric == "SFO" ~ "SFO (FA)",
        metric == "SFO.L" ~ "SFO (FA)",
        metric == "SFO.R" ~ "SFO (FA)",
        #Superior fronto-occipital fasciculus
        metric == "SLF" ~ "SLF (FA)",
        metric == "SLF.L" ~ "SLF (FA)",
        metric == "SLF.R" ~ "SLF (FA)",
        #Superior longitudinal fasciculus
        metric ==  "SS" ~ "SS (FA)",
        metric ==  "SS.L" ~ "SS (FA)",
        metric ==  "SS.R" ~ "SS (FA)",
        #Sagittal Striatum (includes inf. longitudinal fasciculus and inferior-fronto-occipital fasciculus)
        metric == "UNC" ~ "UNC (FA)",
        metric == "UNC.L" ~ "UNC (FA)",
        metric == "UNC.R" ~ "UNC (FA)",
        
        #Uncinate Fasciculus
        TRUE ~ "misc"
      )
  ) %>% 
  
  # Create pFDR column
  # group_by(metric, Hemisphere) %>%
  group_by(DNAm) %>%
  mutate(pFDR = p.adjust(p.value, method = "fdr")) %>%
  # Create a column to denote where there are significant hits 
  mutate(significance =
           case_when(p.value < 0.05 ~ "Yes",
                     TRUE ~ "No")) %>%
  # Create a column to denote where there are FDR significant hits 
  mutate(FDR_significance =
           case_when(pFDR < 0.05 ~ "Yes",
                     TRUE ~ "No"))%>%
  # Create a column to denote where there are FDR significant hits, significant, or not
  mutate(SIG =
           case_when(FDR_significance == "Yes" ~ "pFDR",
                     significance == "Yes" ~ "p",
                     TRUE ~ "not_significant"))

###


plot_dMRI <- newdf

# Cortical volume regressions

FULL_neuroimaging_list <- c(
  
  ###Bilateral ones
  "cv.bilat.bankssts",
  "cv.bilat.caudalanteriorcingulate",
  "cv.bilat.caudalmiddlefrontal",
  "cv.bilat.cuneus",
  "cv.bilat.entorhinal",
  "cv.bilat.frontalpole",
  "cv.bilat.fusiform",
  "cv.bilat.inferiorparietal",
  "cv.bilat.inferiortemporal",
  "cv.bilat.insula",
  "cv.bilat.isthmuscingulate",
  "cv.bilat.lateraloccipital",
  "cv.bilat.lateralorbitofrontal",
  "cv.bilat.lingual",
  "cv.bilat.medialorbitofrontal",
  "cv.bilat.middletemporal",
  "cv.bilat.paracentral",
  "cv.bilat.parahippocampal",
  "cv.bilat.parsopercularis",
  "cv.bilat.parsorbitalis",
  "cv.bilat.parstriangularis",
  "cv.bilat.pericalcarine",
  "cv.bilat.postcentral",
  "cv.bilat.posteriorcingulate",
  "cv.bilat.precentral",
  "cv.bilat.precuneus",
  "cv.bilat.rostralanteriorcingulate",
  "cv.bilat.rostralmiddlefrontal",
  "cv.bilat.superiorfrontal",
  "cv.bilat.superiorparietal",
  "cv.bilat.superiortemporal",
  "cv.bilat.supramarginal",
  "cv.bilat.temporalpole",
  "cv.bilat.transversetemporal",
  
  ### Cortical volume regressions 
  "cv.lh.bankssts",
  "cv.lh.caudalanteriorcingulate",
  "cv.lh.caudalmiddlefrontal",
  "cv.lh.cuneus",
  "cv.lh.entorhinal",
  "cv.lh.frontalpole",
  "cv.lh.fusiform",
  "cv.lh.inferiorparietal",
  "cv.lh.inferiortemporal",
  "cv.lh.insula",
  "cv.lh.isthmuscingulate",
  "cv.lh.lateraloccipital",
  "cv.lh.lateralorbitofrontal",
  "cv.lh.lingual",
  "cv.lh.medialorbitofrontal",
  "cv.lh.middletemporal",
  "cv.lh.paracentral",
  "cv.lh.parahippocampal",
  "cv.lh.parsopercularis",
  "cv.lh.parsorbitalis",
  "cv.lh.parstriangularis",
  "cv.lh.pericalcarine",
  "cv.lh.postcentral",
  "cv.lh.posteriorcingulate",
  "cv.lh.precentral",
  "cv.lh.precuneus",
  "cv.lh.rostralanteriorcingulate",
  "cv.lh.rostralmiddlefrontal",
  "cv.lh.superiorfrontal",
  "cv.lh.superiorparietal",
  "cv.lh.superiortemporal",
  "cv.lh.supramarginal",
  "cv.lh.temporalpole",
  "cv.lh.transversetemporal",
  "cv.rh.bankssts",
  "cv.rh.caudalanteriorcingulate",
  "cv.rh.caudalmiddlefrontal",
  "cv.rh.cuneus",
  "cv.rh.entorhinal",
  "cv.rh.frontalpole",
  "cv.rh.fusiform",
  "cv.rh.inferiorparietal",
  "cv.rh.inferiortemporal",
  "cv.rh.insula",
  "cv.rh.isthmuscingulate",
  "cv.rh.lateraloccipital",
  "cv.rh.lateralorbitofrontal",
  "cv.rh.lingual",
  "cv.rh.medialorbitofrontal",
  "cv.rh.middletemporal",
  "cv.rh.paracentral",
  "cv.rh.parahippocampal",
  "cv.rh.parsopercularis",
  "cv.rh.parsorbitalis",
  "cv.rh.parstriangularis",
  "cv.rh.pericalcarine",
  "cv.rh.postcentral",
  "cv.rh.posteriorcingulate",
  "cv.rh.precentral",
  "cv.rh.precuneus",
  "cv.rh.rostralanteriorcingulate",
  "cv.rh.rostralmiddlefrontal",
  "cv.rh.superiorfrontal",
  "cv.rh.superiorparietal",
  "cv.rh.superiortemporal",
  "cv.rh.supramarginal",
  "cv.rh.temporalpole",
  "cv.rh.transversetemporal"
  
)

# Converting multiple varibles into a factor
global_neuroimaging_dataset %<>% mutate_at(c("sex", "site", "edited", "batch"),
                                           as.factor)

# Making a data frame with all combinations of variables
df <- as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1, metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
get_summary_values <- function(list) {
  metric <- list$metric
  DNAm <- list$DNAm
  tib <-
    broom::tidy(summary(lm(
      scale(global_neuroimaging_dataset[[metric]]) ~ 
        scale(st_age)
      + sex
      + site
      + batch
      + edited
      + scale(Standardised_ICV)
      + scale(global_neuroimaging_dataset[[DNAm]]),
      data = global_neuroimaging_dataset
    )))[8, c(2, 3, 5)]
  return(tib)
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = get_summary_values)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality =
           case_when(
             
             ### Cortical volumes 
             str_detect(metric, "bankssts") ~ "cortical volumes",
             str_detect(metric, "caudalanteriorcingulate") ~ "cortical volumes",
             str_detect(metric, "caudalmiddlefrontal") ~ "cortical volumes",
             str_detect(metric, "cuneus") ~ "cortical volumes",
             str_detect(metric, "entorhinal") ~ "cortical volumes",
             str_detect(metric, "frontalpole")~ "cortical volumes",
             str_detect(metric, "fusiform")~ "cortical volumes",
             str_detect(metric, "inferiorparietal") ~ "cortical volumes",
             str_detect(metric, "inferiortemporal") ~ "cortical volumes",
             str_detect(metric, "insula") ~ "cortical volumes",
             str_detect(metric, "isthmuscingulate") ~ "cortical volumes",
             str_detect(metric, "lateraloccipital") ~ "cortical volumes",
             str_detect(metric, "lateralorbitofrontal") ~ "cortical volumes",
             str_detect(metric, "lingual") ~ "cortical volumes",
             str_detect(metric, "medialorbitofrontal") ~ "cortical volumes",
             str_detect(metric, "middletemporal") ~ "cortical volumes",
             str_detect(metric, "paracentral") ~ "cortical volumes",
             str_detect(metric, "parahippocampal") ~ "cortical volumes",
             str_detect(metric, "parsopercularis") ~ "cortical volumes",
             str_detect(metric, "parsorbitalis") ~ "cortical volumes",
             str_detect(metric, "parstriangularis") ~ "cortical volumes",
             str_detect(metric, "pericalcarine") ~ "cortical volumes",
             str_detect(metric, "postcentral") ~ "cortical volumes",
             str_detect(metric, "posteriorcingulate") ~ "cortical volumes",
             str_detect(metric, "precentral") ~ "cortical volumes",
             str_detect(metric, "precuneus") ~ "cortical volumes",
             str_detect(metric, "rostralanteriorcingulate") ~ "cortical volumes",
             str_detect(metric, "rostralmiddlefrontal") ~ "cortical volumes",
             str_detect(metric, "superiorfrontal") ~ "cortical volumes",
             str_detect(metric, "superiorparietal") ~ "cortical volumes",
             str_detect(metric, "superiortemporal") ~ "cortical volumes",
             str_detect(metric, "supramarginal") ~ "cortical volumes",
             str_detect(metric, "temporalpole") ~ "cortical volumes",
             str_detect(metric, "transversetemporal") ~ "cortical volumes",
             TRUE ~ "misc")
  ) %>%
  
  # Create a hemisphere column
  mutate(Hemisphere =
           case_when(
             str_detect(metric, "\\.lh.") ~ "Left",
             str_detect(metric, "\\.rh.") ~ "Right",
             TRUE ~ "Bilateral"
           )) %>%
  
  
  # Create a brain metric column
  # Create a brain metric column
  mutate(brain_metric =
           case_when(
             
             ### cortical
             str_detect(metric, "bankssts") ~ "banks of superior temporal sulcus",
             str_detect(metric, "caudalanteriorcingulate") ~ "caudal anterior cingulate",
             str_detect(metric, "caudalmiddlefrontal") ~ "caudal middle frontal",
             str_detect(metric, "cv.rh.cuneus") ~ "cuneus",
             str_detect(metric, "cv.lh.cuneus") ~ "cuneus",
             str_detect(metric, "cv.bilat.cuneus") ~ "cuneus",
             str_detect(metric, "entorhinal") ~ "entorhinal",
             str_detect(metric, "frontalpole") ~ "frontal pole",
             str_detect(metric, "fusiform") ~ "fusiform",
             str_detect(metric, "inferiorparietal") ~ "inferior parietal",
             str_detect(metric, "inferiortemporal") ~ "inferior temporal",
             str_detect(metric, "temporalpole") ~ "temporal pole",
             str_detect(metric, "insula") ~ "insula",
             str_detect(metric, "isthmuscingulate") ~ "isthmus cingulate",
             str_detect(metric, "lateraloccipital") ~ "lateral occipital",
             str_detect(metric, "lateralorbitofrontal") ~ "lateral orbitofrontal",
             str_detect(metric, "lingual") ~ "lingual",
             str_detect(metric, "medialorbitofrontal") ~ "medial orbitofrontal",
             str_detect(metric, "middletemporal") ~ "middle temporal",
             str_detect(metric, "paracentral") ~ "paracentral",
             str_detect(metric, "precuneus") ~ "precuneus",
             str_detect(metric, "parahippocampal") ~ "parahippocampal",
             str_detect(metric, "parsopercularis") ~ "pars opercularis",
             str_detect(metric, "parsorbitalis") ~ "pars orbitalis",
             str_detect(metric, "parstriangularis") ~ "pars triangularis",
             str_detect(metric, "pericalcarine") ~ "pericalcarine",
             str_detect(metric, "postcentral") ~ "postcentral",
             str_detect(metric, "posteriorcingulate") ~ "posterior cingulate",
             str_detect(metric, "precentral") ~ "precentral",
             str_detect(metric, "rostralanteriorcingulate") ~ "rostral anterior cingulate",
             str_detect(metric, "rostralmiddlefrontal") ~ "rostral middle frontal",
             str_detect(metric, "superiorfrontal") ~ "superior frontal",
             str_detect(metric, "superiorparietal") ~ "superior parietal",
             str_detect(metric, "superiortemporal") ~ "superior temporal",
             str_detect(metric, "supramarginal") ~ "supra marginal",
             str_detect(metric, "temporal pole") ~ "temporal pole",
             str_detect(metric, "transversetemporal") ~ "transverse temporal",
             TRUE ~ "misc"
           )) %>%
  
  
  
  # Create pFDR column
  # group_by(metric, Hemisphere) %>%
  group_by(DNAm, modality) %>%
  mutate(pFDR = p.adjust(p.value, method = "fdr")) %>%
  # Create a column to denote where there are significant hits 
  mutate(significance =
           case_when(p.value < 0.05 ~ "Yes",
                     TRUE ~ "No")) %>%
  # Create a column to denote where there are FDR significant hits 
  mutate(FDR_significance =
           case_when(pFDR < 0.05 ~ "Yes",
                     TRUE ~ "No"))%>%
  # Create a column to denote where there are FDR significant hits, significant, or not
  mutate(SIG =
           case_when(FDR_significance == "Yes" ~ "pFDR",
                     significance == "Yes" ~ "p",
                     TRUE ~ "not_significant"))

###
#plot %<>% filter(DNAm == "CRP")

plot_cortical_volumes <- newdf

#### Cortical thickness regrssions
FULL_neuroimaging_list <- c(
  
  ###Bilateral ones
  "ct.bilat.bankssts",
  "ct.bilat.caudalanteriorcingulate",
  "ct.bilat.caudalmiddlefrontal",
  "ct.bilat.cuneus",
  "ct.bilat.entorhinal",
  "ct.bilat.frontalpole",
  "ct.bilat.fusiform",
  "ct.bilat.inferiorparietal",
  "ct.bilat.inferiortemporal",
  "ct.bilat.insula",
  "ct.bilat.isthmuscingulate",
  "ct.bilat.lateraloccipital",
  "ct.bilat.lateralorbitofrontal",
  "ct.bilat.lingual",
  "ct.bilat.medialorbitofrontal",
  "ct.bilat.middletemporal",
  "ct.bilat.paracentral",
  "ct.bilat.parahippocampal",
  "ct.bilat.parsopercularis",
  "ct.bilat.parsorbitalis",
  "ct.bilat.parstriangularis",
  "ct.bilat.pericalcarine",
  "ct.bilat.postcentral",
  "ct.bilat.posteriorcingulate",
  "ct.bilat.precentral",
  "ct.bilat.precuneus",
  "ct.bilat.rostralanteriorcingulate",
  "ct.bilat.rostralmiddlefrontal",
  "ct.bilat.superiorfrontal",
  "ct.bilat.superiorparietal",
  "ct.bilat.superiortemporal",
  "ct.bilat.supramarginal",
  "ct.bilat.temporalpole",
  "ct.bilat.transversetemporal",
  
  ### Cortical volume regressions 
  "ct.lh.bankssts",
  "ct.lh.caudalanteriorcingulate",
  "ct.lh.caudalmiddlefrontal",
  "ct.lh.cuneus",
  "ct.lh.entorhinal",
  "ct.lh.frontalpole",
  "ct.lh.fusiform",
  "ct.lh.inferiorparietal",
  "ct.lh.inferiortemporal",
  "ct.lh.insula",
  "ct.lh.isthmuscingulate",
  "ct.lh.lateraloccipital",
  "ct.lh.lateralorbitofrontal",
  "ct.lh.lingual",
  "ct.lh.medialorbitofrontal",
  "ct.lh.middletemporal",
  "ct.lh.paracentral",
  "ct.lh.parahippocampal",
  "ct.lh.parsopercularis",
  "ct.lh.parsorbitalis",
  "ct.lh.parstriangularis",
  "ct.lh.pericalcarine",
  "ct.lh.postcentral",
  "ct.lh.posteriorcingulate",
  "ct.lh.precentral",
  "ct.lh.precuneus",
  "ct.lh.rostralanteriorcingulate",
  "ct.lh.rostralmiddlefrontal",
  "ct.lh.superiorfrontal",
  "ct.lh.superiorparietal",
  "ct.lh.superiortemporal",
  "ct.lh.supramarginal",
  "ct.lh.temporalpole",
  "ct.lh.transversetemporal",
  "ct.rh.bankssts",
  "ct.rh.caudalanteriorcingulate",
  "ct.rh.caudalmiddlefrontal",
  "ct.rh.cuneus",
  "ct.rh.entorhinal",
  "ct.rh.frontalpole",
  "ct.rh.fusiform",
  "ct.rh.inferiorparietal",
  "ct.rh.inferiortemporal",
  "ct.rh.insula",
  "ct.rh.isthmuscingulate",
  "ct.rh.lateraloccipital",
  "ct.rh.lateralorbitofrontal",
  "ct.rh.lingual",
  "ct.rh.medialorbitofrontal",
  "ct.rh.middletemporal",
  "ct.rh.paracentral",
  "ct.rh.parahippocampal",
  "ct.rh.parsopercularis",
  "ct.rh.parsorbitalis",
  "ct.rh.parstriangularis",
  "ct.rh.pericalcarine",
  "ct.rh.postcentral",
  "ct.rh.posteriorcingulate",
  "ct.rh.precentral",
  "ct.rh.precuneus",
  "ct.rh.rostralanteriorcingulate",
  "ct.rh.rostralmiddlefrontal",
  "ct.rh.superiorfrontal",
  "ct.rh.superiorparietal",
  "ct.rh.superiortemporal",
  "ct.rh.supramarginal",
  "ct.rh.temporalpole",
  "ct.rh.transversetemporal"
  
)

# Converting multiple varibles into a factor
global_neuroimaging_dataset %<>% mutate_at(c("sex", "site", "edited", "batch"),
                                           as.factor)

# Making a data frame with all combinations of variables
df <- as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1, metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
# e.g. try
# info <- list(metric = "ACR", DNAm = "ACY1")
# get_summary_values(info)
# It will return you estimate, std.error and p.value in a tibble
get_summary_values <- function(list) {
  metric <- list$metric
  DNAm <- list$DNAm
  tib <-
    broom::tidy(summary(lm(
      scale(global_neuroimaging_dataset[[metric]]) ~ 
        scale(st_age)
      + sex
      + site
      + batch
      + edited
      + scale(Standardised_ICV)
      + scale(global_neuroimaging_dataset[[DNAm]]),
      data = global_neuroimaging_dataset
    )))[8, c(2, 3, 5)]
  return(tib)
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = get_summary_values)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality =
           case_when(
             
             ### Cortical thickness 
             str_detect(metric, "bankssts") ~ "cortical thickness",
             str_detect(metric, "caudalanteriorcingulate") ~ "cortical thickness",
             str_detect(metric, "caudalmiddlefrontal") ~ "cortical thickness",
             str_detect(metric, "cuneus") ~ "cortical thickness",
             str_detect(metric, "entorhinal") ~ "cortical thickness",
             str_detect(metric, "frontalpole")~ "cortical thickness",
             str_detect(metric, "fusiform")~ "cortical thickness",
             str_detect(metric, "inferiorparietal") ~ "cortical thickness",
             str_detect(metric, "inferiortemporal") ~ "cortical thickness",
             str_detect(metric, "insula") ~ "cortical thickness",
             str_detect(metric, "isthmuscingulate") ~ "cortical thickness",
             str_detect(metric, "lateraloccipital") ~ "cortical thickness",
             str_detect(metric, "lateralorbitofrontal") ~ "cortical thickness",
             str_detect(metric, "lingual") ~ "cortical thickness",
             str_detect(metric, "medialorbitofrontal") ~ "cortical thickness",
             str_detect(metric, "middletemporal") ~ "cortical thickness",
             str_detect(metric, "paracentral") ~ "cortical thickness",
             str_detect(metric, "parahippocampal") ~ "cortical thickness",
             str_detect(metric, "parsopercularis") ~ "cortical thickness",
             str_detect(metric, "parsorbitalis") ~ "cortical thickness",
             str_detect(metric, "parstriangularis") ~ "cortical thickness",
             str_detect(metric, "pericalcarine") ~ "cortical thickness",
             str_detect(metric, "postcentral") ~ "cortical thickness",
             str_detect(metric, "posteriorcingulate") ~ "cortical thickness",
             str_detect(metric, "precentral") ~ "cortical thickness",
             str_detect(metric, "precuneus") ~ "cortical thickness",
             str_detect(metric, "rostralanteriorcingulate") ~ "cortical thickness",
             str_detect(metric, "rostralmiddlefrontal") ~ "cortical thickness",
             str_detect(metric, "superiorfrontal") ~ "cortical thickness",
             str_detect(metric, "superiorparietal") ~ "cortical thickness",
             str_detect(metric, "superiortemporal") ~ "cortical thickness",
             str_detect(metric, "supramarginal") ~ "cortical thickness",
             str_detect(metric, "temporalpole") ~ "cortical thickness",
             str_detect(metric, "transversetemporal") ~ "cortical thickness",
             TRUE ~ "misc")
  ) %>%
  
  # Create a hemisphere column
  mutate(Hemisphere =
           case_when(
             str_detect(metric, "\\.lh.") ~ "Left",
             str_detect(metric, "\\.rh.") ~ "Right",
             TRUE ~ "Bilateral"
           )) %>%
  
  
  # Create a brain metric column
  # Create a brain metric column
  mutate(brain_metric =
           case_when(
             
             ### cortical
             str_detect(metric, "bankssts") ~ "banks of superior temporal sulcus",
             str_detect(metric, "caudalanteriorcingulate") ~ "caudal anterior cingulate",
             str_detect(metric, "caudalmiddlefrontal") ~ "caudal middle frontal",
             str_detect(metric, "ct.rh.cuneus") ~ "cuneus",
             str_detect(metric, "ct.lh.cuneus") ~ "cuneus",
             str_detect(metric, "ct.bilat.cuneus") ~ "cuneus",
             str_detect(metric, "entorhinal") ~ "entorhinal",
             str_detect(metric, "frontalpole") ~ "frontal pole",
             str_detect(metric, "fusiform") ~ "fusiform",
             str_detect(metric, "inferiorparietal") ~ "inferior parietal",
             str_detect(metric, "inferiortemporal") ~ "inferior temporal",
             str_detect(metric, "temporalpole") ~ "temporal pole",
             str_detect(metric, "insula") ~ "insula",
             str_detect(metric, "isthmuscingulate") ~ "isthmus cingulate",
             str_detect(metric, "lateraloccipital") ~ "lateral occipital",
             str_detect(metric, "lateralorbitofrontal") ~ "lateral orbitofrontal",
             str_detect(metric, "lingual") ~ "lingual",
             str_detect(metric, "medialorbitofrontal") ~ "medial orbitofrontal",
             str_detect(metric, "middletemporal") ~ "middle temporal",
             str_detect(metric, "paracentral") ~ "paracentral",
             str_detect(metric, "precuneus") ~ "precuneus",
             str_detect(metric, "parahippocampal") ~ "parahippocampal",
             str_detect(metric, "parsopercularis") ~ "pars opercularis",
             str_detect(metric, "parsorbitalis") ~ "pars orbitalis",
             str_detect(metric, "parstriangularis") ~ "pars triangularis",
             str_detect(metric, "pericalcarine") ~ "pericalcarine",
             str_detect(metric, "postcentral") ~ "postcentral",
             str_detect(metric, "posteriorcingulate") ~ "posterior cingulate",
             str_detect(metric, "precentral") ~ "precentral",
             str_detect(metric, "rostralanteriorcingulate") ~ "rostral anterior cingulate",
             str_detect(metric, "rostralmiddlefrontal") ~ "rostral middle frontal",
             str_detect(metric, "superiorfrontal") ~ "superior frontal",
             str_detect(metric, "superiorparietal") ~ "superior parietal",
             str_detect(metric, "superiortemporal") ~ "superior temporal",
             str_detect(metric, "supramarginal") ~ "supra marginal",
             str_detect(metric, "temporal pole") ~ "temporal pole",
             str_detect(metric, "transversetemporal") ~ "transverse temporal",
             TRUE ~ "misc"
           )) %>%
  
  
  
  # Create pFDR column
  # group_by(metric, Hemisphere) %>%
  group_by(DNAm, modality) %>%
  mutate(pFDR = p.adjust(p.value, method = "fdr")) %>%
  # Create a column to denote where there are significant hits 
  mutate(significance =
           case_when(p.value < 0.05 ~ "Yes",
                     TRUE ~ "No")) %>%
  # Create a column to denote where there are FDR significant hits 
  mutate(FDR_significance =
           case_when(pFDR < 0.05 ~ "Yes",
                     TRUE ~ "No"))%>%
  
  # Create a column to denote where there are FDR significant hits, significant, or not
  mutate(SIG =
           case_when(FDR_significance == "Yes" ~ "pFDR",
                     significance == "Yes" ~ "p",
                     TRUE ~ "not_significant"))

###
#plot %<>% filter(DNAm == "CRP")

plot_cortical_thickness <- newdf

### Cortical surface area

FULL_neuroimaging_list <- c(
  
  ###Bilateral ones
  "csa.bilat.bankssts",
  "csa.bilat.caudalanteriorcingulate",
  "csa.bilat.caudalmiddlefrontal",
  "csa.bilat.cuneus",
  "csa.bilat.entorhinal",
  "csa.bilat.frontalpole",
  "csa.bilat.fusiform",
  "csa.bilat.inferiorparietal",
  "csa.bilat.inferiortemporal",
  "csa.bilat.insula",
  "csa.bilat.isthmuscingulate",
  "csa.bilat.lateraloccipital",
  "csa.bilat.lateralorbitofrontal",
  "csa.bilat.lingual",
  "csa.bilat.medialorbitofrontal",
  "csa.bilat.middletemporal",
  "csa.bilat.paracentral",
  "csa.bilat.parahippocampal",
  "csa.bilat.parsopercularis",
  "csa.bilat.parsorbitalis",
  "csa.bilat.parstriangularis",
  "csa.bilat.pericalcarine",
  "csa.bilat.postcentral",
  "csa.bilat.posteriorcingulate",
  "csa.bilat.precentral",
  "csa.bilat.precuneus",
  "csa.bilat.rostralanteriorcingulate",
  "csa.bilat.rostralmiddlefrontal",
  "csa.bilat.superiorfrontal",
  "csa.bilat.superiorparietal",
  "csa.bilat.superiortemporal",
  "csa.bilat.supramarginal",
  "csa.bilat.temporalpole",
  "csa.bilat.transversetemporal",
  
  ### Cortical volume regressions 
  "csa.lh.bankssts",
  "csa.lh.caudalanteriorcingulate",
  "csa.lh.caudalmiddlefrontal",
  "csa.lh.cuneus",
  "csa.lh.entorhinal",
  "csa.lh.frontalpole",
  "csa.lh.fusiform",
  "csa.lh.inferiorparietal",
  "csa.lh.inferiortemporal",
  "csa.lh.insula",
  "csa.lh.isthmuscingulate",
  "csa.lh.lateraloccipital",
  "csa.lh.lateralorbitofrontal",
  "csa.lh.lingual",
  "csa.lh.medialorbitofrontal",
  "csa.lh.middletemporal",
  "csa.lh.paracentral",
  "csa.lh.parahippocampal",
  "csa.lh.parsopercularis",
  "csa.lh.parsorbitalis",
  "csa.lh.parstriangularis",
  "csa.lh.pericalcarine",
  "csa.lh.postcentral",
  "csa.lh.posteriorcingulate",
  "csa.lh.precentral",
  "csa.lh.precuneus",
  "csa.lh.rostralanteriorcingulate",
  "csa.lh.rostralmiddlefrontal",
  "csa.lh.superiorfrontal",
  "csa.lh.superiorparietal",
  "csa.lh.superiortemporal",
  "csa.lh.supramarginal",
  "csa.lh.temporalpole",
  "csa.lh.transversetemporal",
  "csa.rh.bankssts",
  "csa.rh.caudalanteriorcingulate",
  "csa.rh.caudalmiddlefrontal",
  "csa.rh.cuneus",
  "csa.rh.entorhinal",
  "csa.rh.frontalpole",
  "csa.rh.fusiform",
  "csa.rh.inferiorparietal",
  "csa.rh.inferiortemporal",
  "csa.rh.insula",
  "csa.rh.isthmuscingulate",
  "csa.rh.lateraloccipital",
  "csa.rh.lateralorbitofrontal",
  "csa.rh.lingual",
  "csa.rh.medialorbitofrontal",
  "csa.rh.middletemporal",
  "csa.rh.paracentral",
  "csa.rh.parahippocampal",
  "csa.rh.parsopercularis",
  "csa.rh.parsorbitalis",
  "csa.rh.parstriangularis",
  "csa.rh.pericalcarine",
  "csa.rh.postcentral",
  "csa.rh.posteriorcingulate",
  "csa.rh.precentral",
  "csa.rh.precuneus",
  "csa.rh.rostralanteriorcingulate",
  "csa.rh.rostralmiddlefrontal",
  "csa.rh.superiorfrontal",
  "csa.rh.superiorparietal",
  "csa.rh.superiortemporal",
  "csa.rh.supramarginal",
  "csa.rh.temporalpole",
  "csa.rh.transversetemporal"
  
)

# Converting multiple varibles into a factor
global_neuroimaging_dataset %<>% mutate_at(c("sex", "site", "edited", "batch"),
                                           as.factor)

# Making a data frame with all combinations of variables
df <- as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1, metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
get_summary_values <- function(list) {
  metric <- list$metric
  DNAm <- list$DNAm
  tib <-
    broom::tidy(summary(lm(
      scale(global_neuroimaging_dataset[[metric]]) ~ 
        scale(st_age)
      + sex
      + site
      + batch
      + edited
      + scale(Standardised_ICV)
      + scale(global_neuroimaging_dataset[[DNAm]]),
      data = global_neuroimaging_dataset
    )))[8, c(2, 3, 5)]
  return(tib)
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = get_summary_values)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality =
           case_when(
             
             ### Cortical surface area 
             str_detect(metric, "bankssts") ~ "cortical surface area",
             str_detect(metric, "caudalanteriorcingulate") ~ "cortical surface area",
             str_detect(metric, "caudalmiddlefrontal") ~ "cortical surface area",
             str_detect(metric, "cuneus") ~ "cortical surface area",
             str_detect(metric, "entorhinal") ~ "cortical surface area",
             str_detect(metric, "frontalpole")~ "cortical surface area",
             str_detect(metric, "fusiform")~ "cortical surface area",
             str_detect(metric, "inferiorparietal") ~ "cortical surface area",
             str_detect(metric, "inferiortemporal") ~ "cortical surface area",
             str_detect(metric, "insula") ~ "cortical surface area",
             str_detect(metric, "isthmuscingulate") ~ "cortical surface area",
             str_detect(metric, "lateraloccipital") ~ "cortical surface area",
             str_detect(metric, "lateralorbitofrontal") ~ "cortical surface area",
             str_detect(metric, "lingual") ~ "cortical surface area",
             str_detect(metric, "medialorbitofrontal") ~ "cortical surface area",
             str_detect(metric, "middletemporal") ~ "cortical surface area",
             str_detect(metric, "paracentral") ~ "cortical surface area",
             str_detect(metric, "parahippocampal") ~ "cortical surface area",
             str_detect(metric, "parsopercularis") ~ "cortical surface area",
             str_detect(metric, "parsorbitalis") ~ "cortical surface area",
             str_detect(metric, "parstriangularis") ~ "cortical surface area",
             str_detect(metric, "pericalcarine") ~ "cortical surface area",
             str_detect(metric, "postcentral") ~ "cortical surface area",
             str_detect(metric, "posteriorcingulate") ~ "cortical surface area",
             str_detect(metric, "precentral") ~ "cortical surface area",
             str_detect(metric, "precuneus") ~ "cortical surface area",
             str_detect(metric, "rostralanteriorcingulate") ~ "cortical surface area",
             str_detect(metric, "rostralmiddlefrontal") ~ "cortical surface area",
             str_detect(metric, "superiorfrontal") ~ "cortical surface area",
             str_detect(metric, "superiorparietal") ~ "cortical surface area",
             str_detect(metric, "superiortemporal") ~ "cortical surface area",
             str_detect(metric, "supramarginal") ~ "cortical surface area",
             str_detect(metric, "temporalpole") ~ "cortical surface area",
             str_detect(metric, "transversetemporal") ~ "cortical surface area",
             TRUE ~ "misc")
  ) %>%
  
  # Create a hemisphere column
  mutate(Hemisphere =
           case_when(
             str_detect(metric, "\\.lh.") ~ "Left",
             str_detect(metric, "\\.rh.") ~ "Right",
             TRUE ~ "Bilateral"
           )) %>%
  
  
  # Create a brain metric column
  # Create a brain metric column
  mutate(brain_metric =
           case_when(
             
             ### cortical
             str_detect(metric, "bankssts") ~ "banks of superior temporal sulcus",
             str_detect(metric, "caudalanteriorcingulate") ~ "caudal anterior cingulate",
             str_detect(metric, "caudalmiddlefrontal") ~ "caudal middle frontal",
             str_detect(metric, "csa.rh.cuneus") ~ "cuneus",
             str_detect(metric, "csa.lh.cuneus") ~ "cuneus",
             str_detect(metric, "csa.bilat.cuneus") ~ "cuneus",
             str_detect(metric, "entorhinal") ~ "entorhinal",
             str_detect(metric, "frontalpole") ~ "frontal pole",
             str_detect(metric, "fusiform") ~ "fusiform",
             str_detect(metric, "inferiorparietal") ~ "inferior parietal",
             str_detect(metric, "inferiortemporal") ~ "inferior temporal",
             str_detect(metric, "temporalpole") ~ "temporal pole",
             str_detect(metric, "insula") ~ "insula",
             str_detect(metric, "isthmuscingulate") ~ "isthmus cingulate",
             str_detect(metric, "lateraloccipital") ~ "lateral occipital",
             str_detect(metric, "lateralorbitofrontal") ~ "lateral orbitofrontal",
             str_detect(metric, "lingual") ~ "lingual",
             str_detect(metric, "medialorbitofrontal") ~ "medial orbitofrontal",
             str_detect(metric, "middletemporal") ~ "middle temporal",
             str_detect(metric, "paracentral") ~ "paracentral",
             str_detect(metric, "precuneus") ~ "precuneus",
             str_detect(metric, "parahippocampal") ~ "parahippocampal",
             str_detect(metric, "parsopercularis") ~ "pars opercularis",
             str_detect(metric, "parsorbitalis") ~ "pars orbitalis",
             str_detect(metric, "parstriangularis") ~ "pars triangularis",
             str_detect(metric, "pericalcarine") ~ "pericalcarine",
             str_detect(metric, "postcentral") ~ "postcentral",
             str_detect(metric, "posteriorcingulate") ~ "posterior cingulate",
             str_detect(metric, "precentral") ~ "precentral",
             str_detect(metric, "rostralanteriorcingulate") ~ "rostral anterior cingulate",
             str_detect(metric, "rostralmiddlefrontal") ~ "rostral middle frontal",
             str_detect(metric, "superiorfrontal") ~ "superior frontal",
             str_detect(metric, "superiorparietal") ~ "superior parietal",
             str_detect(metric, "superiortemporal") ~ "superior temporal",
             str_detect(metric, "supramarginal") ~ "supra marginal",
             str_detect(metric, "temporal pole") ~ "temporal pole",
             str_detect(metric, "transversetemporal") ~ "transverse temporal",
             TRUE ~ "misc"
           )) %>%
  
  
  
  # Create pFDR column
  # group_by(metric, Hemisphere) %>%
  group_by(DNAm, modality) %>%
  mutate(pFDR = p.adjust(p.value, method = "fdr")) %>%
  # Create a column to denote where there are significant hits 
  mutate(significance =
           case_when(p.value < 0.05 ~ "Yes",
                     TRUE ~ "No")) %>%
  # Create a column to denote where there are FDR significant hits 
  mutate(FDR_significance =
           case_when(pFDR < 0.05 ~ "Yes",
                     TRUE ~ "No"))%>%
  
  # Create a column to denote where there are FDR significant hits, significant, or not
  mutate(SIG =
           case_when(FDR_significance == "Yes" ~ "pFDR",
                     significance == "Yes" ~ "p",
                     TRUE ~ "not_significant"))

###
#plot %<>% filter(DNAm == "CRP")
plot_cortical_surface_area <- newdf

#####



plot <- rbind(plot_global, 
              plot_subcortical, 
              plot_cortical_volumes,
              plot_cortical_thickness,
              plot_cortical_surface_area,
              plot_dMRI)

#######




# Cortical structures, by hemisphere plot

plot2 <- plot

plot2 %<>% filter(DNAm == "CRP")

plot2 %<>% filter(Hemisphere == "Bilateral")

#plot2 %<>% filter(Hemisphere == "Left" |
#                    Hemisphere == "Right")
#
plot2 %<>% filter(modality == "cortical volumes" |
                    modality ==  "cortical thickness" |
                    modality == "cortical surface area")


#plot2 <- plot2[with(plot2, order(brain_metric, modality, -estimate)),]

plot2 %<>%  arrange(metric, desc(estimate)) %>%
  group_by(modality)

plot2$brain_metric <- factor(plot2$brain_metric, levels = unique(plot2$brain_metric))


###### REORDER metric by effect size
plot2 %<>% mutate(modality = factor(
  modality,
  levels = c(
    #"global",
    #        "subcortical volumes",
    "cortical volumes",
    "cortical thickness",
    "cortical surface area"
    # "White matter integrity"
  )
))



#plot2 %<>%  arrange(metric, estimate) %>%
#  mutate(brain_metric = factor(brain_metric, levels = modality)) %>%


############
############

##### for effect sizes plot
facetSettings <-
  theme(strip.background = element_rect(
    fill = "#edf2fb",
    colour = "black",
    size = 1
  ))

ggplot(plot2,
       
       aes(
         #  x = reorder(metric, -(estimate)),
         x = reorder(brain_metric, -(estimate)),
         y = estimate,
          group = metric
         #group = Hemisphere
       )) +
  
  geom_point(
    aes(
      #col = modality,
      col = reorder(brain_metric, -(estimate)),
      #alpha = Hemisphere,
     # alpha = reorder(metric, -(estimate)),
      shape = SIG
    ),
    size = 1.3,
    position = position_dodge(width = 0.9),
    stroke = 0.9
  ) +
  
  
  coord_flip() +
  geom_errorbar(
    aes(
      ymin = estimate - (1.96 * std.error),
      ymax = estimate + (1.96 * std.error)
    ),
    position = position_dodge(0.9),
    width = 0.4,
    colour = "darkgrey",
    alpha = 0.6,
    size = 0.8
  ) +
  
  theme_classic() +
  theme(legend.position = "none") +
  
  theme(
    #axis.text.x = element_blank(),
    # axis.ticks = element_blank(),
    axis.text.x = element_text(
      #  angle = 90,
     # vjust = 0.5,
      #hjust = 1,
      size = 7,
      #face = "bold",
      family = "sans"
    ),
    strip.text = element_text(
      size = 8,
      face = "bold",
      family = "sans",
      colour = "black"
    ),
    axis.text.y = element_text(size = 7)
    #panel.grid.minor = element_line(colour = "grey", linetype = "dashed"),
    
  ) +
  
  labs(x = "",
       y = "") +
  
  geom_hline(
    yintercept = 0,
    color = "darkgrey",
    linetype = "dashed",
    size = 0.2,
    alpha = 0.5
  ) +
  viridis::scale_color_viridis(discrete = TRUE,
                               option = "D",
                               direction = -1) +
  scale_shape_manual(values = c(16,1,5
  )) +
  #  scale_colour_manual(values = c(#"#440154FF",
  #    #                             "#404788FF",
  #    "#2D708EFF",
  #    "#238A8DFF",
  #    "#20A387FF"
  #    # "darkgrey" 
  #    )) +
  facet_wrap( ~ modality,
              #scales="free"
              #scales = "free_x",
              nrow = 1) +
  facetSettings +
  scale_alpha_manual(values=c(0.6, 1)) 