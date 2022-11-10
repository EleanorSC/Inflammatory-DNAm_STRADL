global_neuroimaging_dataset <- merge(STRADL_lifestyle_covariates, 
                                     dataset_n655, 
                                     by = "stradl_ID") 

# Converting multiple varibles into a factor
global_neuroimaging_dataset %<>% mutate_at(c("sex", "site", "edited", "batch"),
                                           as.factor)



# 6000 entries
FULL_DNAm_list <- c(
  "ACY1",
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
  "VEGFA"
)



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
  group_by(DNAm, metric, Hemisphere) %>%
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


# Cortical structures, by hemisphere plot

plot2 <- plot_cortical_volumes

#plot2 %<>% filter(brain_metric == "fusiform")

plot2 %<>% filter(significance == "Yes" |
                    FDR_significance == "Yes")

plot2 %<>% filter(Hemisphere == "Bilateral")
#
#plot2 %<>% filter(Hemisphere == "Left" |
#                    Hemisphere == "Right")

#

### Which DNAm proxy has the most significant hits?
plot2 %<>%
  mutate(number_significant =
           case_when(significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
top_significant_hits <- aggregate(plot2$number_significant,
                                  by = list(brain_metric = plot2$brain_metric),
                                  FUN = sum)

top_significant_hits %<>% arrange(desc(x))


plot2 %<>%  arrange(metric, desc(estimate)) %>%
  group_by(modality)

### Reorder by this
###### REORDER DNAm by no. sig
plot2 %<>% mutate(brain_metric = factor(
  brain_metric,
  levels = c(

                  "parahippocampal", #42
                "inferior temporal", #34
                         "fusiform", #33
                "lateral occipital", #31
                   "supra marginal", #28
            "lateral orbitofrontal", #27
                       "precentral", #26
                "superior parietal", #25
                "superior temporal", #24
                       "entorhinal", #23
                           "insula", #21
                  "middle temporal", #20
                      "postcentral", #20
              "transverse temporal", #20
                          "lingual", #18
           "rostral middle frontal", #18
             "medial orbitofrontal", #17
                   "pars orbitalis", #17
                        "precuneus", #17
       "rostral anterior cingulate", #16
                "inferior parietal", #15
              "posterior cingulate", #14
"banks of superior temporal sulcus", #11
                           "cuneus", #10
                 "pars opercularis", #10
                 "superior frontal", #10
                     "frontal pole", # 9
                "isthmus cingulate", # 8
                    "temporal pole", # 8
            "caudal middle frontal", # 6
                      "paracentral", # 5
        "caudal anterior cingulate", # 3
                "pars triangularis", # 3
                    "pericalcarine"))) # 1

#plot2$brain_metric <- factor(plot2$brain_metric, levels = unique(plot2$brain_metric))

#plot2 %<>%  arrange(metric, estimate) %>%
#  mutate(brain_metric = factor(brain_metric, levels = modality)) %>%


############
############

pal <- colorspace::sequential_hcl(109, palette = "SunsetDark")

##### for effect sizes plot
facetSettings <-
  theme(strip.background = element_rect(
    #fill = "#edf2fb",
    #fill = "#FFFEFC",
    fill = "#F8EEEC",
    colour = "black",
    size = 1
  ))

ggplot(plot2,
       
       aes(
         #  x = reorder(metric, -(estimate)),
         x = reorder(DNAm, -(estimate)),
         y = estimate,
         # group = metric
         group = Hemisphere
       )) +
  
  geom_point(
    aes(
      #col = modality,
      col = reorder(DNAm, (-estimate)),
      alpha = SIG,
      #alpha = reorder(metric, -(estimate)),
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
      vjust = 0.5,
      hjust = 1,
      size = 8,
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
  
#  scale_colour_manual(values = pal
#                      ) +
#  
 viridis::scale_color_viridis(discrete = TRUE,
                              option = "A",
                              direction = -1) +
 
  
# scale_shape_manual(values = c(16,1,5
# )) +
  
  #  scale_colour_manual(values = c(#"#440154FF",
  #    #                             "#404788FF",
  #    "#2D708EFF",
  #    "#238A8DFF",
  #    "#20A387FF"
  #    # "darkgrey" 
  #    )) +
  facet_wrap( ~ brain_metric,
              #scales="free"
              #scales = "free_x",
              nrow = 1) +
  facetSettings +
  scale_alpha_manual(values=c(0.8, 1, 0.2)) +
  scale_shape_manual(values = c(
    # 8,
    # 1,
    16,
    5))


#### try out as a pheWAS


##### SHow less

plot2 %<>% filter(brain_metric == "parahippocampal"|
                    brain_metric == "inferior temporal"|
                  brain_metric == "lateral occipital"|
                    brain_metric == "fusiform"|
                    brain_metric == "superior parietal"|
                    brain_metric == "middle temporal")




##### for effect sizes plot
facetSettings <-
  theme(strip.background = element_rect(
    #fill = "#edf2fb",
    #fill = "#FFFEFC",
    fill = "#F8EEEC",
    colour = "black",
    size = 1
  ))

ggplot(plot2,
       
       aes(
         # x = brain_metric,
         # x = reorder(metric, modality),
         x = reorder(DNAm,
                     (estimate)
                     ), 
         y = -log(p.value) 
         #group = modality
       )
) +
  
  geom_point(aes(
    # col = brain_metric, 
    col = reorder(DNAm,
                  (estimate)
    ), 
    #  col = metric,
    alpha = reorder(DNAm,
                    (estimate)
                    ), 
    #alpha = FDR_significance,
    shape = FDR_significance
    #shape = modality
  ),
  #size = -(estimate)
  size = 2 
  #  alpha = 0.8
  ) +
  
  theme_classic() +
  theme(legend.position = "none") +
  
  theme(
    # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, 
    #                             size = 8)
    axis.text.x = element_blank(),
    strip.text = element_text(
      size = 6,
      face = "bold",
      family="sans",
      colour = "black"),
    axis.text.y = element_text(size = 7),
    #panel.grid.minor = element_line(colour = "grey", linetype = "dashed"),
    axis.ticks = element_blank()
  ) +
  
  labs(
    x = "",
    y = ""
  ) +
  
  ggrepel::geom_text_repel(
    data = . %>% mutate(label = ifelse(p.value < 0.05, 
                                       as.character(DNAm), 
                                       "")
    ),
    aes(label = label),
    segment.size = 0.2, 
    min.segment.length = unit(0.1, "lines"),
    direction = "y",
    size = 1.5,
    #  size = 2.5,
    #nudge_y = 2,
    # nudge_x = 2,
    box.padding = unit(0.9, "lines"),
    max.overlaps = Inf
  ) +
  
  geom_hline(
    yintercept = -log(0.05),
    color = "red",
    linetype = "dashed",
    size = 0.5,
    alpha = 0.5
  ) +
  
  #  geom_hline(
  #    yintercept = -log(0.05),
  #    color = "darkgrey",
  #    size = 1,
  #    alpha = 0.5
  #  ) +
  viridis::scale_color_viridis(discrete = TRUE,
                               option = "A") +
  scale_shape_manual(values = c(
    # 8,
    # 1,
    16,
    5)) +
  # scale_colour_manual(values = c("#7D1D67", 
  #                                "#CA2E7A", 
  #                                "#F75F6D", 
  #                                "darkgrey",
  #                                "#F88379"
  #                               # "#FE8B6F", 
  #                                )) +
  # facet_wrap(~ modality, nrow = 4) +
  # facet_grid(vars(modality), vars(DNAm),
  #            scales="free") +
  facet_wrap(~ brain_metric,
             #scales="free"
             scales="free_x",
             nrow = 2) +
  facetSettings 
#+
#  scale_alpha(range=c(0.8, 1)) 


###################
# REGRESSION PLOT 
###################

#### try out as a pheWAS
plot2 <- plot_cortical_volumes

#plot2 %<>% filter(brain_metric == "fusiform")

plot2 %<>% filter(significance == "Yes" |
                    FDR_significance == "Yes")

plot2 %<>% filter(Hemisphere == "Bilateral")


##### SHow less
#
#plot2 %<>% filter(brain_metric == "parahippocampal"|
#                    brain_metric == "inferior temporal"|
#                    brain_metric == "lateral occipital"|
#                    brain_metric == "fusiform"|
#                    brain_metric == "superior parietal"|
#                    brain_metric == "middle temporal")



###### REORDER DNAm by no. sig
plot2 %<>% mutate(brain_metric = factor(
  brain_metric,
  levels = c(
    
    "parahippocampal", #42
    "inferior temporal", #34
    "fusiform", #33
    "lateral occipital", #31
    "supra marginal", #28
    "lateral orbitofrontal", #27
    "precentral", #26
    "superior parietal", #25
    "superior temporal", #24
    "entorhinal", #23
    "insula", #21
    "middle temporal", #20
    "postcentral", #20
    "transverse temporal", #20
    "lingual", #18
    "rostral middle frontal", #18
    "medial orbitofrontal", #17
    "pars orbitalis", #17
    "precuneus", #17
    "rostral anterior cingulate", #16
    "inferior parietal", #15
    "posterior cingulate", #14
    "banks of superior temporal sulcus", #11
    "cuneus", #10
    "pars opercularis", #10
    "superior frontal", #10
    "frontal pole", # 9
    "isthmus cingulate", # 8
    "temporal pole", # 8
    "caudal middle frontal", # 6
    "paracentral", # 5
    "caudal anterior cingulate", # 3
    "pars triangularis", # 3
    "pericalcarine"))) # 1



##### for effect sizes plot
facetSettings <-
  theme(strip.background = element_rect(
    #fill = "#edf2fb",
    #fill = "#FFFEFC",
    fill = "#F8EEEC",
    colour = "black",
    size = 1
  ))

ggplot(plot2,
       
       aes(x = reorder(DNAm,-(estimate)),
           y = estimate)) +
  
  geom_point(aes(
    col = reorder(DNAm,
                  (estimate)),
    alpha = reorder(DNAm,
                    (estimate)),
    shape = FDR_significance
  ),
  size = 2) +
  
  
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
    axis.text.x = element_text(#angle = 90,
      #vjust = 0.5,
      # hjust=1,
      size = 6),
    #axis.text.x = element_blank(),
    strip.text = element_text(
      size = 6,
      face = "bold",
      family = "sans",
      colour = "black"
    ),
    axis.text.y = element_text(size = 7)
  ) +
  
  labs(x = "",
       y = "") +
  
  geom_hline(
    yintercept = 0,
    color = "lightgrey",
    linetype = "dashed",
    size = 0.3,
    alpha = 0.5
  ) +
  
  viridis::scale_color_viridis(discrete = TRUE,
                               option = "A") +
  scale_shape_manual(values = c(16,
                                5)) +
  facet_wrap( ~ brain_metric,
              #scales="free"
              scales = "free_y",
              nrow = 5) +
  facetSettings
#+
#  scale_alpha(range=c(0.8, 1))

