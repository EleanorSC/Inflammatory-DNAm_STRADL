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

######Load up neuroimaging full dataset & DNAm data (that was cleaned in data prep)
STRADL_FreeSurfer <- read.csv("STRADL_Measures_FreeSurfer_Main.csv")
names(STRADL_FreeSurfer)[names(STRADL_FreeSurfer) == 'id'] <-
  'stradl_ID'
PROTEOMICS_DNAm_DATA <- read.csv("PROTEOMICS_DNAm_DATA.csv")

# n =709
Neuroimaging_DNAm <- merge(PROTEOMICS_DNAm_DATA,
                           STRADL_FreeSurfer,
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
  "scv.bilat.accumbens",  
  "scv.bilat.amygdala",
  "scv.bilat.caudate",  
  "scv.bilat.hippocampus",  
  "scv.bilat.pallidum",  
  "scv.bilat.putamen",  
  "scv.bilat.thalamus",  
  "scv.bilat.ventraldc", 
  "scv.cerebellum_GM",  
  "scv.cerebellum_WM",
  "vol.brainstem")

Neuroimaging_DNAm %<>% mutate(scv.cerebellum_GM =  vol.cerebellum.lh.gm + vol.cerebellum.rh.gm,
                              scv.cerebellum_WM =  vol.cerebellum.lh.wm + vol.cerebellum.rh.wm)

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
  mutate(modality = "subcortical") %>%
  
  # Create a brain metric column
  mutate(
    brain_metric =
      case_when(
        ### cortical
        metric == "scv.bilat.accumbens" ~ "accumbens",
        metric == "scv.bilat.amygdala" ~ "amygdala",
        metric == "scv.bilat.caudate"~ "caudate",
        metric == "scv.bilat.hippocampus"~ "hippocampus",
        metric == "scv.bilat.pallidum"~ "pallidum",
        metric == "scv.bilat.putamen"~ "putamen",
        metric == "scv.bilat.thalamus"~ "thalamus",
        metric == "scv.bilat.ventraldc"~ "ventral diencephalon",
        metric == "vol.bilat.brainstem"~ "brainstem",
        metric == "scv.cerebellum_GM"~ "cerebellum (GM)",
        metric == "scv.cerebellum_WM"~ "cerebellum (WM)",
        metric == "vol.brainstem" ~ "brainstem",
        
        TRUE ~ "misc"
      )
  ) %>%
  
  # Create pFDR column
  # group_by(metric, Hemisphere) %>%
  group_by(brain_metric) %>%
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
plot_subcortical_methylation <- newdf

# ----------------------------#
# Table of significant cortical volume regressions for supplementary document
# ----------------------------#
table_model1 <- plot_subcortical_methylation %>%
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
sigs <- plot_subcortical_methylation %>%
  filter(estimate < 0 & FDR_significance == "Yes") %>%
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
# PLOT - barplot of significant and FDR significant for good brain health measures
# ----------------------------#
plot2 <- plot_subcortical_methylation

Good_brain_health <- plot2 %>% filter(estimate > 0)

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
           "Greater subcortical volume")

#### Poor brain health
Poor_brain_health <- plot2 %>% filter(estimate < 0)

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
           "Lower subcortical volume")

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
       y = "subcortical volumes") +
  
  facet_wrap( ~ direction) +
  
  scale_fill_manual(values = c(
    "#666699",
   # "#7F8D9C", # cerebellum
     # accumbens
    "#4C1D4BFF", #cerebellum GM
    "#751F58FF", #putamen
    "#A11A5BFF", #caudate
   "#AA0144",
    "#CB1B4FFF", #palli
    "#E83F3FFF", #amy
    "#F2704DFF", #thala
    "#F69C73FF", #vdc
    "#F7C5A5FF", #brainstem
    "#FAEBDDFF" # hippo
    
    #viridis::rocket(n = 11)
                               )) +
  scale_colour_manual(values = c("#808080")) 
#+
#  xlim(0, 27)

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
  mutate(modality = "subcortical") %>%
  
  # Create a brain metric column
  mutate(
    brain_metric =
      case_when(
        ### cortical
        metric == "scv.bilat.accumbens" ~ "accumbens",
        metric == "scv.bilat.amygdala" ~ "amygdala",
        metric == "scv.bilat.caudate"~ "caudate",
        metric == "scv.bilat.hippocampus"~ "hippocampus",
        metric == "scv.bilat.pallidum"~ "pallidum",
        metric == "scv.bilat.putamen"~ "putamen",
        metric == "scv.bilat.thalamus"~ "thalamus",
        metric == "scv.bilat.ventraldc"~ "ventral diencephalon",
        metric == "vol.bilat.brainstem"~ "brainstem",
        metric == "scv.cerebellum_GM"~ "cerebellum (GM)",
        metric == "scv.cerebellum_WM"~ "cerebellum (WM)",
        metric == "vol.brainstem" ~ "brainstem",
        
        TRUE ~ "misc"
      )
  ) %>%
  # Create pFDR column
  # group_by(metric, Hemisphere) %>%
  group_by(brain_metric) %>%
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
plot_subcortical_lifestyle_methylation <- newdf


## ----------------------------#
# Combining both sets of models to see percentage attenuation
## ----------------------------#

add_on <- plot_subcortical_lifestyle_methylation %>%
  select(estimate, std.error, r2, p.value, pFDR) %>%
  rename(
    estimate_lifestyle = estimate,
    r2_lifestyle = r2,
    p.value_lifestyle = p.value,
    pFDR_lifestyle = pFDR,
    std.error_lifestyle = std.error
  )

add_on <- subset(add_on, select = -c(brain_metric))

table_new <- cbind(plot_subcortical_methylation, add_on)

table_new %<>% mutate(percentage_increase_decrease =
                        100 * (estimate - estimate_lifestyle)) %>%
  
  # filter(FDR_significance == "Yes" & estimate < 0) %>%
  filter(significance == "Yes") %>%
  
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
                TRUE ~ "No")
  ) %>%
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

write.csv(table_new, "subcortical_volume_regressions_both_models.csv")


table_pFDR <- table_new %>% filter(pFDR < 0.05)

table_pFDR_lyf <- table_pFDR %>% filter(pFDR_lifestyle < 0.05)

table_model1 <- plot_subcortical_methylation %>%
  filter(FDR_significance == "Yes" & estimate > 0)

range(table_model1$estimate)
unique(table_model1$DNAm)