## ---------------------------
##
## Script Purpose: Regressions with individual proteins
##                 Iterates through two lists simultaneously - 
##                 (1) that contains all DNAm proxies &
##                 (2) that contains all neuroimaging metrics of interest
##                 Runs a regression for every Inflammatory-DNAm signature & every neuroimaging metric
##
##                In this script, we for significant associations for white matter tract FA only.
##
##  General notes for neuroimaging data
----------------------------
  #dataset709 <- merge(STRADL_lifestyle_covariates, dataset_n709, by = "stradl_ID")
  #dataset655 <- merge(STRADL_lifestyle_covariates, dataset_n655, by = "stradl_ID") 
  # DTI_2 <-merge(STRADL_lifestyle_covariates,STRADL_DTI, by = "stradl_ID")
  # master_neuroimaging_data <-merge(KORA_LBC_DNAm,DTI_2, by = "stradl_ID")
  # master_neuroimaging_data <- read.csv("master_neuroimaging_data.csv")
  
  #write.csv(dataset709, "STRADL_neuroimaging_n709.csv")  
  #write.csv(dataset655, "STRADL_neuroimaging_n655.csv")
  #write.csv(master_neuroimaging_data, "master_neuroimaging_data.csv") 
  
  #####

#library(dplyr)
#library(stats)
#ICV_corrected <- merge(STRADL_lifestyle_covariates, dataset_n655, by = "stradl_ID") 
#STRADL_DTI <- merge(gFA_gMD, 
#                    STRADL_WM_tracts_FA,
#                    by = "stradl_ID")

master_neuroimaging_data <- merge(ICV_corrected, STRADL_DTI , by = "stradl_ID") 

# Converting multiple varibles into a factor
master_neuroimaging_data %<>% mutate_at(c("sex", "site", "edited", "batch"),
                     as.factor)

skimr::skim(master_neuroimaging_data)

#glimpse(master_neuroimaging_data)

# 6000 entries
FULL_DNAm_list <- c(
                   # "ACY1",
                   # "ADAMTS13",
                   # "ADIPOQ",
                   # "AFM",
                   # "B2M",
                   # "BCAM",
                   # "BMP1",
                   # "C4A|C4B",
                   # "C5",
                   # "C9",
                   # "CCL17",
                   # "CCL18",
                   # "CCL21",
                   # "CCL22",
                   # "CCL25",
                   # "CD163",
                   # "CD209",
                   # "CD48",
                   # "CD5L",
                   # "CHIT1",
                   # "CLEC11A",
                   # "CLEC11A.1",
                   # "CNTN4",
                    "CRP",
                   # "CXCL10.y",
                   # "CXCL11.y",
                   # "EDA",
                   # "ENPP7",
                   # "ESM1",
                   # "F7",
                   # "FAP",
                   # "FCER2",
                   # "FCGR3B",
                   # "GHR",
                   # "GNLY",
                   # "GP1BA",
                   # "GZMA.y",
                   # "HGFAC",
                   # "ICAM5",
                   # "IDUA",
                   # "IGFBP1",
                   # "IGFBP4",
                   # "IL19",
                   # "INSR",
                   # "LGALS3BP", 
                   # "LGALS4",
                   # "LTA|LTB",
                   # "LTF",
                   # "LY9",
                   # "LYZ",
                   # "MIA",
                   # "MMP1",
                   # "MMP12",
                   # "MMP2",
                   # "MMP9",
                   # "MPL",
                   # "MPO",
                   # "MRC2",
                   # "MST1",
                   # "NCAM1",
                   # "NOTCH1",
                   # "NTRK3.y",
                   # "OMD",
                   # "PAPPA",
                   # "PIGR",
                   # "PRSS2",
                   # "RARRES2",
                   # "RETN",
                   # "S100A9",
                   # "SELE",
                   # "SELL",
                   # "SEMA3E",
                   # "SERPINA3",
                   # "SERPIND1",
                   # "SHBG",
                   # "SLITRK5",
                   # "SPOCK2",
                   # "STC1",
                   # "THBS2",
                   # "TNFRSF17",
                   # "TNFRSF1B",
                   # "TPSB2",
                   # "VCAM1",
                   # "WFIKKN2",
                   # 
                   # # LBC trained proxies
                   # ##### WHICH STRADL KORA AND STRADL LBC match?
                   # # GZMA r = 0.71, MMP.1 r = 0.46, CXCL10 r = 0.35, NTRK3 r = 0.26, and CXCL11 r = 0.09
                   # # x = LBC, y = KORA
                   # 
                   # "CCL11",
                   # "CD6",
                   # "CRTAM",
                   # "CXCL10.x",
                   # "CXCL11.x",
                   # "CXCL9",
                   # "EN.RAGE",
                   # "EZR",
                   # "FcRL2",
                   # "FGF.21",
                   # "G.CSF",
                   # "GDF.8",
                   # "GZMA.x",
                   # "HGF",
                   # "MMP.1",
                   # "N.CDase",
                   # "NEP",
                   # "NMNAT1",
                   # "NTRK3.x",
                   # "OSM",
                   # "SIGLEC1",
                   # "SKR3",
                   # "SMPD1",
                   # "TGF.alpha",
                    "VEGFA"
                    )

FULL_neuroimaging_list <- c(### Global measures
                            "global.cerebral.wm",
                            "global.total.gm",
                            "global.wbv",
                            "global.wbv.vent",
                            "gFA",
                            "gMD",
  
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
                            "cv.rh.transversetemporal",
                            
                            ### subcortical volumes
                            
                            "scv.lh.accumbens",
                            "scv.lh.amygdala",
                            "scv.lh.caudate",
                            "scv.lh.hippocampus",
                            "scv.lh.pallidum",
                            "scv.lh.putamen",
                            "scv.lh.thalamus",
                            "scv.lh.ventraldc",
                            "scv.rh.accumbens",
                            "scv.rh.amygdala",
                            "scv.rh.caudate",
                            "scv.rh.hippocampus",
                            "scv.rh.pallidum",
                            "scv.rh.putamen",
                            "scv.rh.thalamus",
                            "scv.rh.ventraldc",
                            "vol.brainstem",
                            "vol.cerebellum.lh.gm",
                            "vol.cerebellum.lh.wm",
                            "vol.cerebellum.rh.gm",
                            "vol.cerebellum.rh.wm"
                            
                            ### White matter tract FA
                            
                        #    "ACR",  
                        #    "ACR.L",
                        #    "ACR.R",
                        #    "ALIC",
                        #    "ALIC.L",
                        #    "ALIC.R",
                        #    "BCC",
                        #    "CC",
                        #    "CGC",
                        #    "CGC.L",
                        #    "CGC.R",
                        #    "CGH",
                        #    "CGH.L",
                        #    "CGH.R",
                        #    "CR",
                        #    "CR.L",
                        #    "CR.R",
                        #    "CST",
                        #    "CST.L",
                        #    "CST.R",
                        #    "EC",
                        #    "EC.L",
                        #    "EC.R",
                        #    "FX",
                        #    "FX.ST.L",
                        #    "FX.ST.R",
                        #    "FX.ST",
                        #    "GCC",
                        #    "IC",
                        #    "IC.L",
                        #    "IC.R",
                        #    "IFO",
                        #    "IFO.L",
                        #    "IFO.R",
                        #    "PCR",
                        #    "PCR.L",
                        #    "PCR.R",
                        #    "PLIC",
                        #    "PLIC.L",
                        #    "PLIC.R",
                        #    "PTR",
                        #    "PTR.L",
                        #    "PTR.R",
                        #    "RLIC",
                        #    "RLIC.L",
                        #    "RLIC.R",
                        #    "SCC",
                        #    "SCR",
                        #    "SCR.L",
                        #    "SCR.R",
                        #    "SFO",
                        #    "SFO.L",
                        #    "SFO.R",
                        #    "SLF",
                        #    "SLF.L",
                        #    "SLF.R",
                        #    "SS",
                        #    "SS.L",
                        #    "SS.R",
                        #    "UNC",
                        #    "UNC.L",
                        #    "UNC.R"
                        #    
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
      scale(master_neuroimaging_data[[metric]]) ~ scale(st_age)
      + sex
      + site
      + batch
      + edited
      + scale(Standardised_ICV)
      + scale(master_neuroimaging_data[[DNAm]]),
      data = master_neuroimaging_data
    )))[8, c(2, 3, 5)]
  return(tib)
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  # This splits input data frame into nrow(df) named lists
  # of the form list(metric = {metric}, DNAm = {DNAm})
  split(1:nrow(.)) %>%
  # map_dfr applies the get_summary_values function to each of these lists in turn,
  # and saves the output as a new row with 3 columns
  # in the newcols dataframe (NB dfr stands for data frame row)
  purrr::map_dfr(.f = get_summary_values)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%

  # Create a neurimaging modality column
  mutate(modality =
           case_when(
             str_detect(metric, "global") ~ "global",
             str_detect(metric, "gFA") ~ "global",
             str_detect(metric, "gMD") ~ "global",
             
             ### Sub-cortical volumes 
             str_detect(metric, "scv.") ~ "subcortical volumes",
             str_detect(metric, "vol.") ~ "subcortical volumes",
             
             ### Cortical volumes 
             metric == "cv.lh.bankssts" ~ "cortical volumes",
             metric == "cv.lh.caudalanteriorcingulate" ~ "cortical volumes",
             metric == "cv.lh.caudalmiddlefrontal"~ "cortical volumes",
             metric == "cv.lh.cuneus" ~ "cortical volumes",
             metric == "cv.lh.entorhinal" ~ "cortical volumes",
             metric == "cv.lh.frontalpole"~ "cortical volumes",
             metric == "cv.lh.fusiform"~ "cortical volumes",
             metric == "cv.lh.inferiorparietal" ~ "cortical volumes",
             metric == "cv.lh.inferiortemporal" ~ "cortical volumes",
             metric == "cv.lh.insula" ~ "cortical volumes",
             metric == "cv.lh.isthmuscingulate" ~ "cortical volumes",
             metric == "cv.lh.lateraloccipital" ~ "cortical volumes",
             metric == "cv.lh.lateralorbitofrontal" ~ "cortical volumes",
             metric == "cv.lh.lingual" ~ "cortical volumes",
             metric == "cv.lh.medialorbitofrontal" ~ "cortical volumes",
             metric == "cv.lh.middletemporal" ~ "cortical volumes",
             metric == "cv.lh.paracentral" ~ "cortical volumes",
             metric == "cv.lh.parahippocampal" ~ "cortical volumes",
             metric == "cv.lh.parsopercularis" ~ "cortical volumes",
             metric == "cv.lh.parsorbitalis" ~ "cortical volumes",
             metric == "cv.lh.parstriangularis" ~ "cortical volumes",
             metric == "cv.lh.pericalcarine" ~ "cortical volumes",
             metric == "cv.lh.postcentral" ~ "cortical volumes",
             metric == "cv.lh.posteriorcingulate" ~ "cortical volumes",
             metric == "cv.lh.precentral" ~ "cortical volumes",
             metric == "cv.lh.precuneus" ~ "cortical volumes",
             metric == "cv.lh.rostralanteriorcingulate" ~ "cortical volumes",
             metric == "cv.lh.rostralmiddlefrontal" ~ "cortical volumes",
             metric == "cv.lh.superiorfrontal" ~ "cortical volumes",
             metric == "cv.lh.superiorparietal" ~ "cortical volumes",
             metric == "cv.lh.superiortemporal" ~ "cortical volumes",
             metric == "cv.lh.supramarginal" ~ "cortical volumes",
             metric == "cv.lh.temporalpole" ~ "cortical volumes",
             metric == "cv.lh.transversetemporal" ~ "cortical volumes",
             metric == "cv.rh.bankssts" ~ "cortical volumes",
             metric == "cv.rh.caudalanteriorcingulate" ~ "cortical volumes",
             metric == "cv.rh.caudalmiddlefrontal" ~ "cortical volumes",
             metric == "cv.rh.cuneus" ~ "cortical volumes",
             metric == "cv.rh.entorhinal" ~ "cortical volumes",
             metric == "cv.rh.frontalpole" ~ "cortical volumes",
             metric == "cv.rh.fusiform" ~ "cortical volumes",
             metric == "cv.rh.inferiorparietal" ~ "cortical volumes",
             metric == "cv.rh.inferiortemporal" ~ "cortical volumes",
             metric == "cv.rh.insula" ~ "cortical volumes",
             metric == "cv.rh.isthmuscingulate" ~ "cortical volumes",
             metric == "cv.rh.lateraloccipital" ~ "cortical volumes",
             metric == "cv.rh.lateralorbitofrontal" ~ "cortical volumes",
             metric == "cv.rh.lingual" ~ "cortical volumes",
             metric == "cv.rh.medialorbitofrontal" ~ "cortical volumes",
             metric == "cv.rh.middletemporal" ~ "cortical volumes",
             metric == "cv.rh.paracentral" ~ "cortical volumes",
             metric == "cv.rh.parahippocampal" ~ "cortical volumes",
             metric == "cv.rh.parsopercularis" ~ "cortical volumes",
             metric == "cv.rh.parsorbitalis" ~ "cortical volumes",
             metric == "cv.rh.parstriangularis" ~ "cortical volumes",
             metric == "cv.rh.pericalcarine" ~ "cortical volumes",
             metric == "cv.rh.postcentral" ~ "cortical volumes",
             metric == "cv.rh.posteriorcingulate" ~ "cortical volumes",
             metric == "cv.rh.precentral" ~ "cortical volumes",
             metric == "cv.rh.precuneus" ~ "cortical volumes",
             metric == "cv.rh.rostralanteriorcingulate" ~ "cortical volumes",
             metric == "cv.rh.rostralmiddlefrontal" ~ "cortical volumes",
             metric == "cv.rh.superiorfrontal" ~ "cortical volumes",
             metric == "cv.rh.superiorparietal" ~ "cortical volumes",
             metric == "cv.rh.superiortemporal" ~ "cortical volumes",
             metric == "cv.rh.supramarginal" ~ "cortical volumes",
             metric == "cv.rh.temporalpole" ~ "cortical volumes",
             metric == "cv.rh.transversetemporal" ~ "cortical volumes",
             
             TRUE ~ "misc"
           )) %>%
  
  # Create a hemisphere column
  mutate(Hemisphere =
           case_when(
             str_detect(metric, "\\.L") ~ "Left",
             str_detect(metric, "\\.R") ~ "Right",
             str_detect(metric, "\\.lh") ~ "Left",
             str_detect(metric, "\\.rh") ~ "Right",
             TRUE ~ "Bilateral"
           )) %>%
  # Create pFDR column
  # group_by(metric, Hemisphere) %>%
  group_by(modality, Hemisphere, DNAm) %>%
  mutate(pFDR = p.adjust(p.value, method = "fdr")) %>%
  # Create a column to denote where there are significant hits 
  mutate(significance =
           case_when(p.value < 0.05 ~ "Yes",
                     TRUE ~ "No"))

###

plot <- newdf


ggplot(plot, aes(x = metric, 
                y = -log(p.value)
                )
       ) +
  
  geom_point(aes(col = modality),
                 #size = -(estimate)
                 size = 1.2, 
                 alpha = 0.8
                 
             ) +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_blank(),
    panel.grid.minor = element_line(colour = "grey", linetype = "dashed"),
    axis.ticks = element_blank()
  ) +
  
  labs(color = "Category",
       #size = "Effect size",
       x = "neuroimaging metric",
       y = "log(p-value)") +
  
  ggrepel::geom_text_repel(
    data = . %>% mutate(label = ifelse(p.value < 0.01, as.character(metric), "")),
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
  facet_wrap(~ DNAm, nrow = 2)