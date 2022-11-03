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
  # DTI_DNAm <-merge(KORA_LBC_DNAm,DTI_2, by = "stradl_ID")
  # DTI_DNAm <- read.csv("DTI_DNAm.csv")
  
  #write.csv(dataset709, "STRADL_neuroimaging_n709.csv")  
  #write.csv(dataset655, "STRADL_neuroimaging_n655.csv")
  #write.csv(DTI_DNAm, "DTI_DNAm.csv") 
  
  #####

#library(dplyr)
#library(stats)

# 6000 entries
FULL_DNAm_list <- c("ACY1",
                    "ADAMTS13",
                    "ADIPOQ",
                    "AFM",
                    "B2M",
                    "BCAM",
                    "BMP1",
                    "C4A|C4B",
                    "C5",
                    "C9",
                    "CCL17",
                    "CCL18",
                    "CCL21",
                    "CCL22",
                    "CCL25",
                    "CD163",
                    "CD209",
                    "CD48",
                    "CD5L",
                    "CHIT1",
                    "CLEC11A",
                    "CLEC11A.1",
                    "CNTN4",
                    "CRP",
                    "CXCL10.y",
                    "CXCL11.y",
                    "EDA",
                    "ENPP7",
                    "ESM1",
                    "F7",
                    "FAP",
                    "FCER2",
                    "FCGR3B",
                    "GHR",
                    "GNLY",
                    "GP1BA",
                    "GZMA.y",
                    "HGFAC",
                    "ICAM5",
                    "IDUA",
                    "IGFBP1",
                    "IGFBP4",
                    "IL19",
                    "INSR",
                    "LGALS3BP", 
                    "LGALS4",
                    "LTA|LTB",
                    "LTF",
                    "LY9",
                    "LYZ",
                    "MIA",
                    "MMP1",
                    "MMP12",
                    "MMP2",
                    "MMP9",
                    "MPL",
                    "MPO",
                    "MRC2",
                    "MST1",
                    "NCAM1",
                    "NOTCH1",
                    "NTRK3.y",
                    "OMD",
                    "PAPPA",
                    "PIGR",
                    "PRSS2",
                    "RARRES2",
                    "RETN",
                    "S100A9",
                    "SELE",
                    "SELL",
                    "SEMA3E",
                    "SERPINA3",
                    "SERPIND1",
                    "SHBG",
                    "SLITRK5",
                    "SPOCK2",
                    "STC1",
                    "THBS2",
                    "TNFRSF17",
                    "TNFRSF1B",
                    "TPSB2",
                    "VCAM1",
                    "WFIKKN2",
                    
                    # LBC trained proxies
                    ##### WHICH STRADL KORA AND STRADL LBC match?
                    # GZMA r = 0.71, MMP.1 r = 0.46, CXCL10 r = 0.35, NTRK3 r = 0.26, and CXCL11 r = 0.09
                    # x = LBC, y = KORA
                    
                    "CCL11",
                    "CD6",
                    "CRTAM",
                    "CXCL10.x",
                    "CXCL11.x",
                    "CXCL9",
                    "EN.RAGE",
                    "EZR",
                    "FcRL2",
                    "FGF.21",
                    "G.CSF",
                    "GDF.8",
                    "GZMA.x",
                    "HGF",
                    "MMP.1",
                    "N.CDase",
                    "NEP",
                    "NMNAT1",
                    "NTRK3.x",
                    "OSM",
                    "SIGLEC1",
                    "SKR3",
                    "SMPD1",
                    "TGF.alpha",
                    "VEGFA")
                       
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
                            "UNC.R")


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
      scale(DTI_DNAm[[metric]]) ~ scale(st_age)
      + sex
      + Site
      + scale(DTI_DNAm[[DNAm]]),
      data = DTI_DNAm
    )))[5, c(2, 3, 5)]
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
  # Create a hemisphere column
  mutate(Hemisphere =
           case_when(
             str_detect(metric, "\\.L") ~ "Left",
             str_detect(metric, "\\.R") ~ "Right",
             TRUE ~ "Bilateral"
           )) %>%
  # Create pFDR column
 # group_by(metric, Hemisphere) %>%
  group_by(metric, Hemisphere, DNAm) %>%
  mutate(pFDR = p.adjust(p.value, method = "fdr")) %>%
  # Create a column to denote where there are significant hits 
  mutate(significance =
           case_when(p.value < 0.05 ~ "Yes",
                     TRUE ~ "No"))

###

plot <- newdf

# Create a column "brain_var" where WM are labelled

plot %<>%
  mutate(
    brain_var =
      case_when(
        metric == "ACR" ~ "ACR",
        #Anterior corona radiata
        metric == "ACR.L" ~ "ACR",
        #Anterior corona radiata
        metric == "ACR.R" ~ "ACR",
        #Anterior corona radiata
        metric == "ALIC" ~ "ALIC",
        #Anterior Limb of Internal Capsule
        metric == "ALIC.L" ~ "ALIC",
        #Anterior Limb of Internal Capsule
        metric == "ALIC.R" ~ "ALIC",
        #Anterior Limb of Internal Capsule
        metric == "CR" ~ "CR",
        #Corona Radiata
        metric == "CR.L" ~ "CR",
        #Corona Radiata
        metric == "CR.R" ~ "CR",
        #Corona Radiata
        metric == "IC" ~ "IC",
        #Internal Capsule
        metric == "IC.L" ~ "IC",
        #Internal Capsule
        metric == "IC.R" ~ "IC",
        #Internal Capsule
        metric == "RLIC" ~ "RLIC",
        #Retrolenticular part of internal capsule
        metric == "RLIC.L" ~ "RLIC",
        #Retrolenticular part of internal capsule
        metric == "RLIC.R" ~ "RLIC",
        #Retrolenticular part of internal capsule
        
        str_detect(metric, "BCC") ~ "BCC",
        # #Body of corpus callosum
        str_detect(metric, "CC") ~ "CC",
        # # Genu of corpus callosum
        str_detect(metric, "CGC") ~ "CGC",
        #Cingulum-Cingulate Gyrus
        str_detect(metric, "CGH") ~ "CGH",
        #Cingulum-Hippocampus
        
        str_detect(metric, "CST") ~ "CST",
        #Corticospinal tract
        str_detect(metric, "EC") ~ "EC",
        #External Capsule
        str_detect(metric, "FX") ~ "FX",
        #Fornix-Column&Body
        str_detect(metric, "FX.ST") ~ "FX.ST",
        #Fornix-Cres/Stria Terminalis #note I changed FXST -> FX.ST
        str_detect(metric, "GCC") ~ "GCC",
        #Genu of corpus callosum
        
        str_detect(metric, "IFO") ~ "IFO",
        #Inferior fronto-occipital fasciculus
        str_detect(metric, "PCR") ~ "PCR",
        #Posterior corona radiata
        str_detect(metric, "PLIC") ~ "PLIC",
        #Posterior limb of internal capsule
        str_detect(metric, "PTR") ~ "PTR",
        #Posterior Thalamic Radiation
        
        str_detect(metric, "SCC") ~ "SCC",
        # Splenium of corpus callosum
        str_detect(metric, "SCR") ~ "SCR",
        #Superior corona radiata
        str_detect(metric, "SFO") ~ "SFO",
        #Superior fronto-occipital fasciculus
        str_detect(metric, "SLF") ~ "SLF",
        #Superior longitudinal fasciculus
        str_detect(metric, "SS") ~ "SS",
        #Sagittal Striatum (includes inf. longitudinal fasciculus and inferior-fronto-occipital fasciculus)
        str_detect(metric, "UNC") ~ "UNC",
        #Uncinate Fasciculus
        TRUE ~ "misc"
      )
  )

### Plot based on significance so make a new column like is p > 0.05 then concat with protein

plot %<>%
  mutate(significance =
           case_when(p.value < 0.05 ~ "Yes",
                     TRUE ~ "No"))


### Which DNAm proxy has the most significant hits?
plot %<>%
  mutate(number_significant =
           case_when(significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
top_significant_hits <- aggregate(plot$number_significant,
                                  by = list(DNAm = plot$DNAm),
                                  FUN = sum)

top_significant_hits %<>% arrange(desc(x))


plot$Hemisphere <- as.factor(plot$Hemisphere)
plot2 <- plot %>% filter(Hemisphere == "Right" | Hemisphere == "Left")


###### REORDER DNAm by no. sig
plot2 %<>% mutate(DNAm = factor(
  DNAm,
  levels = c(
    "SEMA3E",    "PIGR",      "SERPIND1",  "VEGFA",     "MMP12",     "MMP.1",     "NTRK3.x",   "SKR3",      "CCL11",    
     "CRP",       "NCAM1",     "NTRK3.y",   "CCL17",     "CCL18",     "FGF.21",    "ICAM5",     "PRSS2",     "SLITRK5",  
     "WFIKKN2",   "HGF",       "BCAM",      "IGFBP4",    "CNTN4",     "NOTCH1",    "OMD",       "TGF.alpha", "INSR",     
     "SPOCK2",    "AFM",       "CXCL11.x",  "GDF.8",     "GZMA.y",    "SIGLEC1",   "ESM1",      "MMP9",      "EDA",      
     "ACY1",      "C4A|C4B",   "CCL22",     "ENPP7",     "GZMA.x",    "IL19",      "NEP",       "SELL",      "LYZ",      
     "MMP1",      "CD209",     "CD5L",      "CD6",       "SHBG",      "THBS2",     "CD48",      "EN.RAGE",   "F7",       
     "LTA|LTB",   "MMP2",      "SELE",      "STC1",      "VCAM1",     "LTF",       "OSM",       "SERPINA3",  "ADAMTS13", 
     "CRTAM",     "GNLY",      "LGALS3BP",  "LGALS4",    "MIA" ,      "NMNAT1",    "RARRES2",   "ADIPOQ",    "BMP1",     
     "C9",        "CCL25",     "CLEC11A",   "CLEC11A.1", "CXCL9",     "FcRL2",     "G.CSF",     "MPO",       "TPSB2",   
     "B2M",       "C5",        "CCL21",     "CD163",     "CHIT1",     "CXCL10.x",  "CXCL10.y",  "CXCL11.y",  "EZR",      
     "FAP",       "FCER2",     "FCGR3B",    "GHR",       "GP1BA",     "HGFAC",     "IDUA",      "IGFBP1",    "LY9",      
     "MPL",       "MRC2",      "MST1",      "N.CDase",   "PAPPA",     "RETN",      "S100A9",    "SMPD1",     "TNFRSF17", 
     "TNFRSF1B")))



##### Now as a point plot instead of barplot

x <- ggplot(
  plot2,
  aes(
    x = reorder(brain_var,-estimate),
    y = estimate,
    colour = significance,
    group = Hemisphere,
    shape = Hemisphere
  )
) +
  
  geom_point(position = position_dodge(width = 0.9),
             size = 2.0,
             stroke = 0.9) +
  
  geom_errorbar(
    aes(
      ymin = estimate - (1.96 * std.error),
      ymax = estimate + (1.96 * std.error)
    ),
    position = position_dodge(0.9),
    width = 0.2,
    colour = "darkgrey",
    alpha = 0.9,
    size = 0.8
  ) +
  
  geom_hline(yintercept = 0, linetype = "dotted") +
  
  coord_flip() +
  #theme_bw() +
  theme_classic() +
  xlab("White matter tract FA") +
  ylab("Standardised effect size") +
  #scale_y_continuous(limits = c(-0.2, 0.2)) +
  facet_wrap( ~ DNAm, nrow = 5)


x +
  theme(
    axis.title.x = element_text(
      size = 11,
      face = "bold",
      colour = "black"
    ),
    axis.text.y = element_text(size = 9,
                               colour = "black"),
    axis.text.x = element_text(size = 6,
                               colour = "black"),
    axis.title.y = element_text(
      size = 11,
      face = "bold",
      colour = "black"
    )) +
      scale_shape_manual(values = c(16,
                                    1)) +
      scale_colour_manual(values = c("darkgrey",
                                    # "cornflowerblue"
                                    # "#404788FF"
                                    "#CC0A4E"))
    