## ---------------------------
##
## Script Purpose: Looking at different global brain and cognitive associations with DNAm 
##
##                Things to note: 
##                (1) models have to be run separately for cognitive vs neuroimaging metrics owing to
##                Different covariates needing to be controlled for: 
##                (A) Cognitive models control for age and sex only
##                (B) Neuroimaging models must contain key covariates of (site / edits / batch / estimated ICV) alongside age and sex
##                (C) Relative brain age is calculated by taking brain age and regressing it on chronological age, hence only sex is controlled here: 
##
##                prot$brain_accel <- resid(lm(Brain_age ~ st_age, na.action = na.exclude, data = prot))                
##                
##                
##                
##                
##                 
##                 
##                 
##
##                
##
##  
## ----------------------------#

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
  "LTA|LTB"      ,
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
  "C4A|C4B"        ,
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


######Load up neuroimaging full dataset
STRADL_FreeSurfer <- read.csv("STRADL_Measures_FreeSurfer_Main.csv")
# Now link these to STRADL ID
names(STRADL_FreeSurfer)[names(STRADL_FreeSurfer) == 'id'] <- 'stradl_ID'


# n =709
Neuroimaging_DNAm <- merge(PROTEOMICS_DNAm_DATA, 
                           STRADL_FreeSurfer,
                           by = "stradl_ID")

# Invert the polarity of these measures
Neuroimaging_DNAm %<>% mutate(gFA =
                                gFA*-1,
                              gMD = 
                                gMD*-1)


Neuroimaging_DNAm %<>% mutate(global_cortical_surface_area = 
                                hem.lh.csa + hem.rh.csa,
                              global_cortical_thickness = 
                                hem.lh.ct + hem.rh.ct,
                              global_cortical_volume = 
                                hem.lh.cv + hem.rh.cv,
                              global_subcortical_volume =
                                scv.bilat.accumbens + scv.bilat.amygdala + scv.bilat.caudate +
                                scv.bilat.hippocampus + scv.bilat.pallidum + scv.bilat.putamen +
                                scv.bilat.thalamus + scv.bilat.ventraldc
                              ) 


skimr::skim(Neuroimaging_DNAm)

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

# Converting multiple varibles into a factor
Neuroimaging_DNAm %<>% mutate_at(c("sex", "site", "edited", "batch"),
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
      scale(Neuroimaging_DNAm[[metric]]) ~ 
        scale(st_age)
      + sex
      + site
      + batch
      + edited
      + scale(est.icv.BAD)
      + scale(Neuroimaging_DNAm[[DNAm]]),
      data = Neuroimaging_DNAm
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
             metric == "global.cerebral.wm" ~ "global",
             metric == "global.total.gm" ~ "global",
             metric == "global.wbv" ~ "global",
             
             metric == "global_cortical_surface_area" ~ "cortical",
             metric == "global_cortical_thickness"~ "cortical",
             metric == "global_cortical_volume"~ "cortical",
             metric == "global_subcortical_volume"~ "subcortical",
             
             metric == "gFA"~ "global",
             metric == "gMD"~ "global",
             metric == "Fazekas_Score_Total" ~ "global",
             
             TRUE ~ "misc")
  ) %>%
  
  # Create a brain metric column
  mutate(brain_metric =
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
             
             TRUE ~ "misc")
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
  mutate(omic_type = "DNAm")

###
plot_neuroimaging_methylation <- newdf

## ----------------------------# 
# COMBINE ANALYSES TABLES
## ----------------------------#

##### REORDER cognitive by global then domains
plot2 <- plot_neuroimaging_methylation

## ----------------------------# 
# WRITE TABLE OF RESULTS FOR SUPPLEMENTARY
## ----------------------------#
table2 <- plot2 %>% 
  group_by(brain_metric) %>% 
  arrange(brain_metric, 
          estimate,
          by_group = TRUE) %>% 
  mutate(CI_lower =
           estimate - (1.96 * std.error),
         CI_upper =
           estimate + (1.96 * std.error)) %>%
  select(DNAm, brain_metric, estimate,std.error,CI_lower, CI_upper, p.value, pFDR)


#write.csv(table, "supplementary_table_global_neuroimaging.csv")

plot2 %<>% filter(significance == "Yes" &
                    brain_metric == "total brain volume" |
                    brain_metric == "global grey matter" |
                    brain_metric == "global white matter"| 
                    brain_metric == "global cortical volume" |
                    brain_metric == "global subcortical volume" |
                    brain_metric == "gFA" |
                    brain_metric == "gMD" |
                    brain_metric == "WMH")

##### REORDER cognitive by global then domains
plot2 %<>% mutate(brain_metric =
                    factor(
                      brain_metric,
                      levels = c(
                        
                        "total brain volume",
                        "global grey matter",
                        "global white matter",
                        "global cortical volume",
                        "global subcortical volume",
                        "gFA",
                        "gMD",
                        "WMH"
                      )
                    ))


##### for effect sizes plot
facetSettings <-
  theme(strip.background = element_rect(
    fill = "#EAE3F2", #purple
    colour = "black",
    size = 1
  ))

ggplot(plot2,
       
       aes(
         x = reorder(DNAm,(-estimate)),
         y = estimate,
         # alpha = reorder(DNAm,
         #                 (-estimate)),
         shape = FDR_significance,
         #col = DNAm,
         col = reorder(DNAm,
                       (estimate)),
         group = omic_type
         # alpha = significance
       )
) +
  
  geom_point(position = position_dodge(width = 0.9),
             size = 1.6,
             stroke = 0.9) +
  
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
  coord_flip() +
  theme(legend.position = "none") +
  
  theme(
    axis.text.x = element_text(
      size = 6),
    strip.text = element_text(
      size = 6,
      face = "bold",
      family = "sans",
      colour = "black"
    ),
    axis.text.y = element_text(size = 7),
    axis.title.x =element_text(
      size = 8,
      face = "bold",
      family = "sans",
      colour = "black"),
    axis.title.y =element_text(
      size = 8,
      face = "bold",
      family = "sans",
      colour = "black")
  ) +
  
  labs(y = "Standardized effect size",
       x = "DNAm signature") +
  
  geom_hline(
    yintercept = 0,
    color = "lightgrey",
    linetype = "dashed",
    size = 0.3,
    alpha = 0.5
  ) +
  
   viridis::scale_color_viridis(discrete = TRUE,
                                option = "F") +
   
  scale_shape_manual(values = c(1,
                                16)
  ) +
  facet_wrap(~ brain_metric,
             nrow = 1) +
  facetSettings 