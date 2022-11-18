## ---------------------------
##
## Script Purpose: Cortical volume associations with DNAm metrics 
##                
##                (nb.) Neuroimaging models must contain key covariates of (site / edits / batch / estimated ICV) alongside age and sex
## 
##                 This script does the following:
##
##                (1) Creates supplementary table for pFDR significant cortical volume associations;
##                (2) Creates barplot for pFDR significant cortical volume associations facetted by whether they associate with lower or increased volumes
##
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
  "cv.bilat.transversetemporal"
  
)


# Converting multiple varibles into a factor
Neuroimaging_DNAm %<>% mutate_at(c("sex", "site", "edited", "batch"),
                                 as.factor)

# Making a data frame with all combinations of variables
df <-
  as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1,
                metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
get_summary_values <- function(list) {
  metric <- list$metric
  DNAm <- list$DNAm
  tib <-
    broom::tidy(summary(
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
    ))[8, c(2, 3, 5)]
  return(tib)
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = get_summary_values)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality = "cortical") %>%
  
  # Create a brain metric column
  mutate(
    brain_metric =
      case_when(
        ### cortical
        str_detect(metric, "bankssts") ~ "banks of superior temporal sulcus",
        str_detect(metric, "caudalanteriorcingulate") ~ "caudal anterior cingulate",
        str_detect(metric, "caudalmiddlefrontal") ~ "caudal middle frontal",
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
plot_neuroimaging_methylation <- newdf

# ----------------------------#
# Table of significant cortical volume regressions for supplementary document
# ----------------------------#
table <- plot_neuroimaging_methylation %>%
  filter(FDR_significance == "Yes") %>%
  group_by(DNAm, brain_metric) %>%
  arrange(brain_metric,
          estimate,
          by_group = TRUE) %>%
  mutate(CI_lower =
           estimate - (1.96 * std.error),
         CI_upper =
           estimate + (1.96 * std.error)) %>%
  select(brain_metric, DNAm, estimate, CI_lower, CI_upper, p.value, pFDR)

# ----------------------------#
# Examine n number for this analysis
# ----------------------------#
Neuroimaging_DNAm %<>% mutate(global_cortical_volume =
                                hem.lh.cv + hem.rh.cv)

skimr::skim(Neuroimaging_DNAm$global_cortical_volume)

write.csv(table, "cortical_volume_regressions_n709.csv")

# ----------------------------#
# Table of beta values
# ----------------------------#
table <- plot_neuroimaging_methylation %>%
  filter(estimate < 0 & FDR_significance == "Yes") %>%
  group_by(DNAm, brain_metric) %>%
  select(brain_metric, DNAm, estimate, std.error, p.value, pFDR)
#filter(DNAm == "MMP12")
#arrange(estimate, groups = TRUE) %>%
summarize(mean_estimate = mean(estimate, na.rm = TRUE))


### Which DNAm proxy has the most FDR significant hits?
sigs <- plot_neuroimaging_methylation %>%
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
plot2 <- plot_neuroimaging_methylation

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
           "Increased cortical volume")

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
           "Decreased cortical volume")


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
       y = "cortical volumes") +
  
  facet_wrap( ~ direction) +
  
  scale_fill_manual(values = c(viridis::mako(n = 34))) +
  scale_colour_manual(values = c("#808080"))
  

