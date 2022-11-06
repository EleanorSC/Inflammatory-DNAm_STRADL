## ---------------------------
##
## Script Purpose: Code for plot of significant DNAm associations with global grey matter 
##                 Uses analyses from 'global_neuroimaging_metrics'
##                 
##                 plot <- read.csv("plot.csv")
##                 
##
##                
##
##  
## ---------------------------
## ---------------------------
##
## Script Purpose: Regressions with individual DNAms for GLOBAL METRICS SIGNIFICANT ONLY
##                 Iterates through two lists simultaneously - 
##                 (1) that contains all DNAm proxies &
##                 (2) that contains all neuroimaging metrics of interest
##                 Runs a regression for every DNAm signature & every neuroimaging metric specified 
##
##                
##
##  Global_neuroimaging_metrics analysis script:
## ----------------------------
##                  global_neuroimaging_dataset <- merge(STRADL_lifestyle_covariates, 
##                                                       dataset_n655, 
##                                                       by = "stradl_ID") 
##                  
##                  
##                  ### Examine global measures, add L and R hemispheres
##                  
##                  global_neuroimaging_dataset %<>% mutate(global_cortical_surface_area = 
##                                                            hem.lh.csa + hem.rh.csa) 
##                  
##                  global_neuroimaging_dataset %<>% mutate(global_cortical_thickness = 
##                                                            hem.lh.ct + hem.rh.ct) 
##                  
##                  global_neuroimaging_dataset %<>% mutate(global_cortical_volume = 
##                                                            hem.lh.cv + hem.rh.cv) 
##                  
##                  # Converting multiple varibles into a factor
##                  global_neuroimaging_dataset %<>% mutate_at(c("sex", "site", "edited", "batch"),
##                                                             as.factor)
##                  
##                  skimr::skim(global_neuroimaging_dataset)
##                  
##                  FULL_neuroimaging_list <- c(
##                    ### Global measures
##                    "global.cerebral.wm",
##                    "global.total.gm",
##                    "global.wbv",
##                    # "global.wbv.vent",
##                    "global_cortical_volume",
##                    "global_cortical_thickness",
##                    "global_cortical_surface_area"
##                  )
##                  
##                  
##                  # 6000 entries
##                  FULL_DNAm_list <- c("ACY1",
##                                      "ADAMTS13",
##                                      "ADIPOQ",
##                                      "AFM",
##                                      "B2M",
##                                      "BCAM",
##                                      "BMP1",
##                                      "C4A|C4B",
##                                      "C5",
##                                      "C9",
##                                      "CCL17",
##                                      "CCL18",
##                                      "CCL21",
##                                      "CCL22",
##                                      "CCL25",
##                                      "CD163",
##                                      "CD209",
##                                      "CD48",
##                                      "CD5L",
##                                      "CHIT1",
##                                      "CLEC11A",
##                                      "CLEC11A.1",
##                                      "CNTN4",
##                                      "CRP",
##                                      "CXCL10.y",
##                                      "CXCL11.y",
##                                      "EDA",
##                                      "ENPP7",
##                                      "ESM1",
##                                      "F7",
##                                      "FAP",
##                                      "FCER2",
##                                      "FCGR3B",
##                                      "GHR",
##                                      "GNLY",
##                                      "GP1BA",
##                                      "GZMA.y",
##                                      "HGFAC",
##                                      "ICAM5",
##                                      "IDUA",
##                                      "IGFBP1",
##                                      "IGFBP4",
##                                      "IL19",
##                                      "INSR",
##                                      "LGALS3BP", 
##                                      "LGALS4",
##                                      "LTA|LTB",
##                                      "LTF",
##                                      "LY9",
##                                      "LYZ",
##                                      "MIA",
##                                      "MMP1",
##                                      "MMP12",
##                                      "MMP2",
##                                      "MMP9",
##                                      "MPL",
##                                      "MPO",
##                                      "MRC2",
##                                      "MST1",
##                                      "NCAM1",
##                                      "NOTCH1",
##                                      "NTRK3.y",
##                                      "OMD",
##                                      "PAPPA",
##                                      "PIGR",
##                                      "PRSS2",
##                                      "RARRES2",
##                                      "RETN",
##                                      "S100A9",
##                                      "SELE",
##                                      "SELL",
##                                      "SEMA3E",
##                                      "SERPINA3",
##                                      "SERPIND1",
##                                      "SHBG",
##                                      "SLITRK5",
##                                      "SPOCK2",
##                                      "STC1",
##                                      "THBS2",
##                                      "TNFRSF17",
##                                      "TNFRSF1B",
##                                      "TPSB2",
##                                      "VCAM1",
##                                      "WFIKKN2",
##                                      
##                                      # LBC trained proxies
##                                      ##### WHICH STRADL KORA AND STRADL LBC match?
##                                      # GZMA r = 0.71, MMP.1 r = 0.46, CXCL10 r = 0.35, NTRK3 r = 0.26, and CXCL11 r = 0.09
##                                      # x = LBC, y = KORA
##                                      
##                                      "CCL11",
##                                      "CD6",
##                                      "CRTAM",
##                                      "CXCL10.x",
##                                      "CXCL11.x",
##                                      "CXCL9",
##                                      "EN.RAGE",
##                                      "EZR",
##                                      "FcRL2",
##                                      "FGF.21",
##                                      "G.CSF",
##                                      "GDF.8",
##                                      "GZMA.x",
##                                      "HGF",
##                                      "MMP.1",
##                                      "N.CDase",
##                                      "NEP",
##                                      "NMNAT1",
##                                      "NTRK3.x",
##                                      "OSM",
##                                      "SIGLEC1",
##                                      "SKR3",
##                                      "SMPD1",
##                                      "TGF.alpha",
##                                      "VEGFA")
##                  
##                  
##                  # Making a data frame with all combinations of variables
##                  df <- as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
##                    dplyr::rename(DNAm = Var1, metric = Var2)
##                  
##                  # Function to get summary values as tibble from a named list with the info on metric and DNAm
##                  # e.g. try
##                  # info <- list(metric = "ACR", DNAm = "ACY1")
##                  # get_summary_values(info)
##                  # It will return you estimate, std.error and p.value in a tibble
##                  get_summary_values <- function(list) {
##                    metric <- list$metric
##                    DNAm <- list$DNAm
##                    tib <-
##                      broom::tidy(summary(lm(
##                        scale(global_neuroimaging_dataset[[metric]]) ~ 
##                          scale(st_age)
##                        + sex
##                        + site
##                        + batch
##                        + edited
##                        + scale(Standardised_ICV)
##                        + scale(global_neuroimaging_dataset[[DNAm]]),
##                        data = global_neuroimaging_dataset
##                      )))[8, c(2, 3, 5)]
##                    return(tib)
##                  }
##                  
##                  # Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
##                  newcols <- df %>%
##                    split(1:nrow(.)) %>%
##                    purrr::map_dfr(.f = get_summary_values)
##                  
##                  # Now bind on the 3 cols that were created
##                  newdf <- cbind(df, newcols) %>%
##                    
##                    # Create a neurimaging modality column
##                    mutate(modality =
##                             case_when(
##                               str_detect(metric, "global") ~ "global",
##                               TRUE ~ "misc")
##                    ) %>%
##                    
##                    # Create a brain metric column
##                    mutate(brain_metric =
##                             case_when(
##                               
##                               ### Global
##                               metric == "global.cerebral.wm" ~ "global white matter",
##                               metric == "global.total.gm" ~ "global grey matter",
##                               metric == "global.wbv" ~ "total brain volume",
##                               #metric == "global.wbv.vent" ~ "global ventricles volume",
##                               metric == "global_cortical_volume" ~ "global cortical volume",
##                               metric == "global_cortical_thickness" ~ "global cortical thickness",
##                               metric == "global_cortical_surface_area" ~ "global cortical surface area",
##                               
##                               TRUE ~ "misc")
##                    ) %>%
##                    
##                    # Create pFDR column
##                    # group_by(metric, Hemisphere) %>%
##                    group_by(DNAm) %>%
##                    mutate(pFDR = p.adjust(p.value, method = "fdr")) %>%
##                    # Create a column to denote where there are significant hits 
##                    mutate(significance =
##                             case_when(p.value < 0.05 ~ "Yes",
##                                       TRUE ~ "No")) %>%
##                    # Create a column to denote where there are FDR significant hits 
##                    mutate(FDR_significance =
##                             case_when(pFDR < 0.05 ~ "Yes",
##                                       TRUE ~ "No"))
##                  
##                  ###
##
##
## ---------------------------
## RUN SCRIPT FROM HERE 
## ---------------------------

#plot <- read.csv("plot.csv")

plot <- newdf

plot2 <- plot %>% filter(brain_metric == "global grey matter")



### Make colours be related to protein + significance
plot2$concat <- paste(plot2$significance, plot2$DNAm)


plot2 %<>% mutate(concat =
                    case_when(str_detect(concat, "No") ~ "not significant",
                              TRUE ~ concat))

###### REORDER DNAm by effect size
plot3 <-  plot2 %>% filter(brain_metric == "global grey matter" &
                             significance == "Yes")
plot3 %<>% arrange(estimate)


### Only plot significant associations
plot2 %<>% filter(DNAm == "MMP12" |
                    DNAm == "PIGR" |
                    DNAm == "IGFBP4"|
                    DNAm == "THBS2"|
                    DNAm == "VEGFA"|
                    DNAm == "CRP"|
                    DNAm == "RARRES2"|
                    DNAm == "SKR3"|
                    DNAm == "SERPIND1" |
                    DNAm == "PRSS2"|
                    DNAm == "CCL18"|
                    DNAm == "FGF.21"|
                    DNAm == "HGF" |
                    DNAm == "TGF.alpha" |
                    DNAm == "MMP1" |
                    DNAm == "ICAM5"|
                    DNAm == "MMP.1"|
                    DNAm == "CCL17"|
                    DNAm == "MMP9"|
                    DNAm == "CXCL10.y" |
                    DNAm == "G.CSF" | 
                    DNAm == "LGALS4"|
                    DNAm == "BCAM"| 
                    DNAm == "GP1BA"| 
                    DNAm == "GZMA.x"|
                    DNAm == "GNLY"|
                    DNAm == "MMP2"|
                    DNAm == "OMD" |  
                    DNAm == "MRC2"|  
                    DNAm == "NTRK3.x"|
                    DNAm == "ESM1"|
                    DNAm == "SLITRK5"|
                    DNAm == "SELL"|
                    DNAm == "SEMA3E"|
                    DNAm == "NOTCH1"|
                    DNAm == "GDF.8" |
                    DNAm == "CNTN4" |
                    DNAm == "NTRK3.y"|
                    DNAm == "NCAM1") 


# We need to add in ESM1 which WM is significant for, add after NTRK3.x

plot2 %<>% mutate(concat = factor(
  concat,
  levels = c(  "Yes MMP12",
               "Yes PIGR",
               "Yes IGFBP4",
               "Yes THBS2",
               "Yes VEGFA",
               "Yes CRP",
               "Yes RARRES2",  
               "Yes SKR3",
               "Yes SERPIND1",
               "Yes PRSS2",
               "Yes CCL18",
               "Yes FGF.21",
               "Yes HGF",
               "Yes TGF.alpha",
               "Yes MMP1",
               "Yes ICAM5",
               "Yes MMP.1",
               "Yes CCL17",
               "Yes MMP9",
               "Yes CXCL10.y",
               "Yes G.CSF",    
               "Yes LGALS4",
               "Yes BCAM",
               # "Yes ESM1",
               "Yes GP1BA",
               "Yes GZMA.x",
               "Yes GNLY",
               "Yes MMP2",
               "Yes OMD" ,     
               "Yes MRC2",  
               "Yes NTRK3.x",
               
               "Yes SLITRK5",
               "Yes SELL",
               "Yes SEMA3E",
               "Yes NOTCH1",
               "Yes GDF.8" ,   
               "Yes CNTN4" ,
               "Yes NTRK3.y",
               "Yes NCAM1", 
               
               
               "not significant")
)
)


x <- ggplot(
  plot2,
  aes(
    x = reorder(DNAm,-estimate),
    y = estimate,
    colour = concat,
    group = DNAm,
    shape = FDR_significance
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
    width = 0.4,
    colour = "darkgrey",
    alpha = 0.6,
    size = 0.8
  ) +
  
  geom_hline(yintercept = 0, linetype = "dotted") +
  
  coord_flip() +
  theme_classic() +
  xlab("DNAm signatures") +
  ylab("Standardised effect size") +
  facet_wrap( ~ brain_metric, nrow = 1) 



x +
  theme(
    axis.title.x = element_text(
      size = 11,
      face = "bold",
      colour = "black"
    ),
    axis.text.y = element_text(size = 9,
                               colour = "black"),
    axis.text.x = element_text(size = 9,
                               colour = "black"),
    axis.title.y = element_text(
      size = 11,
      face = "bold",
      colour = "black"
    )) +
  scale_shape_manual(values = c(1,
                                16)) +
  
  scale_colour_manual(values = c("#7D1D67", "#841E6A", "#8C1E6D", "#931E6F", "#9A1F72", "#A12074", "#A82276", "#AF2477",
                                 "#B62678", "#BD2879", "#C42B79", "#CA2E7A", "#D13279", "#D73679", "#DD3A78", "#E33E76",
                                 "#E94274", "#EF4772", "#F44C6F", "#F6566E", "#F75F6D", "#F8676C", "#FA6F6C", "#FB766C",
                                 "#FC7D6C", "#FD846D", "#FE8B6F", "#FF9271", "#FF9873", "#FF9E75", "#FFA578", "#FFAB7C",
                                 "#FFB17F", "#FFB783", "#FFBD88", "#FFC28C", "#FFC891", "#FFCE96", 
                                 #"#FFD39B", 
                                 # "#FFD99F",
                                 "darkgrey")) +
  theme(legend.position="none")
