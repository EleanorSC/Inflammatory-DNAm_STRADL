## ---------------------------
##
## Script Purpose: Looking at different global brain and cognitive associations with DNAm 

##                
##                (A) Neuroimaging models must contain key covariates of (site / edits / batch / estimated ICV) alongside age and sex
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


#write.csv(PROTEOMICS_DNAm_DATA, "PROTEOMICS_DNAm_DATA.csv")
PROTEOMICS_DNAm_DATA <- read.csv("PROTEOMICS_DNAm_DATA.csv")

# n =709
Neuroimaging_DNAm <- merge(PROTEOMICS_DNAm_DATA, 
                           STRADL_FreeSurfer,
                           by = "stradl_ID")


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
  mutate(modality = "cortical") %>%
  
  # Create a brain metric column
  mutate(brain_metric =
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
#%>%
#  mutate(model = "Model 1 (n=709)")

###
plot_neuroimaging_methylation <- newdf

plot2 <- plot_neuroimaging_methylation

#plot2 %<>% filter(significance == "Yes")

####


##### for effect sizes plot
facetSettings <-
  theme(strip.background = element_rect(
    fill = "#EAE3F2", #purple
    colour = "black",
    size = 1
  ))


ggplot(plot2, aes(x = brain_metric, 
                 y = -log(p.value),
                 shape = FDR_significance,
                 alpha = significance,
                 colour = DNAm
)
) +
  
  geom_point(
    #aes(col = DNAm),
             #size = -(estimate)
             size = 1.2, 
             alpha = 0.8
             
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_blank(),
    panel.grid.minor = element_line(colour = "grey", linetype = "dashed"),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) +
  
  labs(color = "Category",
       #size = "Effect size",
       x = "neuroimaging metric",
       y = "log(p-value)") +
  
  ggrepel::geom_text_repel(
    data = . %>% mutate(label = ifelse(p.value < 0.01, as.character(brain_metric), "")),
    aes(label = label),
    size = 3.1,
    box.padding = unit(0.7, "lines"),
    max.overlaps = Inf
  ) +
  
  geom_hline(
    yintercept = -log(0.01),
    color = "red",
    size = 1,
    alpha = 0.5
  ) +
  
  geom_hline(
    yintercept = -log(0.05),
    color = "darkgrey",
    size = 1,
    alpha = 0.5
  ) +  
  facet_wrap(~ DNAm) +
  
  viridis::scale_color_viridis(discrete = TRUE,
                               option = "F")












#
#ggplot(plot2,
#       
#       aes(
#         x = reorder(DNAm,(-estimate)),
#         y = estimate,
#         # alpha = reorder(DNAm,
#         #                 (-estimate)),
#         shape = FDR_significance,
#         #col = DNAm,
#         col = reorder(DNAm,
#                       (estimate)),
#         group = omic_type
#         # alpha = significance
#       )
#) +
#  
#  geom_point(position = position_dodge(width = 0.9),
#             size = 1.6,
#             stroke = 0.9) +
#  
#  geom_errorbar(
#    aes(
#      ymin = estimate - (1.96 * std.error),
#      ymax = estimate + (1.96 * std.error)
#    ),
#    position = position_dodge(0.9),
#    width = 0.4,
#    colour = "darkgrey",
#    alpha = 0.6,
#    size = 0.8
#  ) +
#  
#  theme_classic() +
#  coord_flip() +
#  theme(legend.position = "none") +
#  
#  theme(
#    axis.text.x = element_text(
#      size = 6),
#    strip.text = element_text(
#      size = 6,
#      face = "bold",
#      family = "sans",
#      colour = "black"
#    ),
#    axis.text.y = element_text(size = 7),
#    axis.title.x =element_text(
#      size = 8,
#      face = "bold",
#      family = "sans",
#      colour = "black"),
#    axis.title.y =element_text(
#      size = 8,
#      face = "bold",
#      family = "sans",
#      colour = "black")
#  ) +
#  
#  labs(y = "Standardized effect size",
#       x = "DNAm signature") +
#  
#  geom_hline(
#    yintercept = 0,
#    color = "lightgrey",
#    linetype = "dashed",
#    size = 0.3,
#    alpha = 0.5
#  ) +
#  
#  viridis::scale_color_viridis(discrete = TRUE,
#                               option = "F") +
#  
#  scale_shape_manual(values = c(1,
#                                16)
#  ) +
#  facet_wrap(~ brain_metric,
#             nrow = 1) +
#  facetSettings 
#
#
#