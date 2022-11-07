## ---------------------------
##
## Script Purpose: Regressions with individual DNAms for  SUBCORTICAL METRICS
##                 Iterates through two lists simultaneously - 
##                 (1) that contains all DNAm proxies &
##                 (2) that contains all neuroimaging metrics of interest
##                 Runs a regression for every DNAm signature & every neuroimaging metric specified 
##
##                
##
##  General notes for neuroimaging data
## ----------------------------

######Load up neuroimaging full dataset
#STRADL_FreeSurfer <- read.csv("STRADL_Measures_FreeSurfer_Main.csv")
## Now link these to STRADL ID
#names(STRADL_FreeSurfer)[names(STRADL_FreeSurfer) == 'id'] <- 'stradl_ID'
#
#
#STRADL_ICV <- read.csv("STRADL_Measures_Standardised_ICV.csv")
## Now link these to STRADL ID
#names(STRADL_ICV)[names(STRADL_ICV) == 'ID'] <- 'stradl_ID'
#
#
#STRADL_MRI <- merge(STRADL_ICV, 
#                    STRADL_FreeSurfer,
#                    by = "stradl_ID")
#
#dataset_n655 <- merge(KORA_LBC_DNAm, 
#                      STRADL_MRI,
#                      by = "stradl_ID")
#
global_neuroimaging_dataset <- merge(STRADL_lifestyle_covariates, 
                                     dataset_n655, 
                                     by = "stradl_ID") 


### Examine global measures, add L and R hemispheres

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


# Converting multiple varibles into a factor
global_neuroimaging_dataset %<>% mutate_at(c("sex", "site", "edited", "batch"),
                                           as.factor)

skimr::skim(global_neuroimaging_dataset)

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
             str_detect(metric, "global") ~ "global",
             TRUE ~ "misc")
  ) %>%
  
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
                     TRUE ~ "No"))

###

#### Work out hippocampal order
## Work out for which values are significant and plot these only
#plot2 %<>% filter(brain_metric == "hippocampus" )
#plot2 %<>% arrange(desc(estimate))
#plot2$DNAm
plot <- newdf


####
plot2 <- plot


### Make colours be related to protein + significance
plot2$concat <- paste(plot2$significance, plot2$DNAm)


plot2 %<>% mutate(concat =
                    case_when(str_detect(concat, "No") ~ "not significant",
                              TRUE ~ concat))
###### REORDER DNAm by estimate hippocampus

plot2 %<>% filter(brain_metric == "hippocampus" & significance == "Yes")

# How many unique colours are needed?
pal <- colorspace::sequential_hcl(25, palette = "SunsetDark")

x <- ggplot(
  plot2,
  aes(
    x = reorder(DNAm,-estimate),
    # x = DNAm,
    y = estimate,
    #colour = reorder(DNAm, estimate),
    colour = brain_metric,
    alpha = reorder(DNAm, -estimate),
    group = DNAm,
    shape = significance
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
  scale_shape_manual(values = c(16,
                                16)) +
  scale_colour_manual(values = c("#429E9D")) +
  #viridis::scale_color_viridis(discrete = TRUE) +
  theme(legend.position="none") +
  scale_y_continuous(limits= c(-0.2, 0.2))


# scale_colour_manual(values = c( 
#                                "darkgrey",
#                                "#7D1D67")) 
#
#theme(legend.position="none")
#plot2$DNAm
plot <- newdf


####
plot2 <- plot


### Make colours be related to protein + significance
plot2$concat <- paste(plot2$significance, plot2$DNAm)


plot2 %<>% mutate(concat =
                    case_when(str_detect(concat, "No") ~ "not significant",
                              TRUE ~ concat))
###### REORDER DNAm by estimate hippocampus

plot2 %<>% filter(brain_metric == "cerebellum (GM)" & significance == "Yes")

plot2 %<>% arrange(desc(estimate))

# How many unique colours are needed?
#pal <- colorspace::sequential_hcl(25, palette = "SunsetDark")

x <- ggplot(
  plot2,
  aes(
    x = reorder(DNAm,-estimate),
    # x = DNAm,
    y = estimate,
    colour = brain_metric,
    alpha = reorder(DNAm, -estimate),
    group = DNAm,
    shape = significance
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
  scale_shape_manual(values = c(16,
                                16)) +
  scale_colour_manual(values = c("#69C8A1")) +
  #viridis::scale_color_viridis(discrete = TRUE) +
  theme(legend.position="none") +
  scale_y_continuous(limits= c(-0.25, 0.2))








#plot2$DNAm
plot <- newdf


####
plot2 <- plot


### Make colours be related to protein + significance
plot2$concat <- paste(plot2$significance, plot2$DNAm)


plot2 %<>% mutate(concat =
                    case_when(str_detect(concat, "No") ~ "not significant",
                              TRUE ~ concat))
###### REORDER DNAm by estimate hippocampus

plot2 %<>% filter(brain_metric == "amygdala" & significance == "Yes")

plot2 %<>% arrange(desc(estimate))

# How many unique colours are needed?
#pal <- colorspace::sequential_hcl(25, palette = "SunsetDark")

x <- ggplot(
  plot2,
  aes(
    x = reorder(DNAm,-estimate),
    # x = DNAm,
    y = estimate,
    #colour = reorder(DNAm, estimate),
    colour = brain_metric,
    alpha = reorder(DNAm, -estimate),
    group = DNAm,
    shape = significance
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
  scale_shape_manual(values = c(16,
                                16)) +
  scale_colour_manual(values = c("#981F71")) +
  #viridis::scale_color_viridis(discrete = TRUE) +
  theme(legend.position="none") +
  scale_y_continuous(limits= c(-0.25, 0.2))








plot <- newdf


####
plot2 <- plot


### Make colours be related to protein + significance
plot2$concat <- paste(plot2$significance, plot2$DNAm)


plot2 %<>% mutate(concat =
                    case_when(str_detect(concat, "No") ~ "not significant",
                              TRUE ~ concat))
###### REORDER DNAm by estimate hippocampus

plot2 %<>% filter(brain_metric == "caudate" & significance == "Yes")

# How many unique colours are needed?


x <- ggplot(
  plot2,
  aes(
    x = reorder(DNAm,-estimate),
    # x = DNAm,
    y = estimate,
    #colour = reorder(DNAm, estimate),
    colour = brain_metric,
    alpha = reorder(DNAm, -estimate),
    group = DNAm,
    shape = significance
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
  scale_shape_manual(values = c(16,
                                16)) +
  scale_colour_manual(values = c("#bdb0d0")) +
  #viridis::scale_color_viridis(discrete = TRUE) +
  theme(legend.position="none") +
  scale_y_continuous(limits= c(-0.2, 0.2))


# scale_colour_manual(values = c( 
#                                "darkgrey",
#                                "#7D1D67")) 
#
#theme(legend.position="none")
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
  scale_shape_manual(values = c(16,
                                16)) +
  scale_colour_manual(values = c("#aed580")) +
  #viridis::scale_color_viridis(discrete = TRUE) +
  theme(legend.position="none") +
  scale_y_continuous(limits= c(-0.2, 0.2))


# scale_colour_manual(values = c( 
#                                "darkgrey",
#                                "#7D1D67")) 
#
#theme(legend.position="none")
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
  scale_shape_manual(values = c(16,
                                16)) +
  scale_colour_manual(values = c("#aed580")) +
  #viridis::scale_color_viridis(discrete = TRUE) +
  theme(legend.position="none") +
  scale_y_continuous(limits= c(-0.2, 0.2))


# scale_colour_manual(values = c( 
#                                "darkgrey",
#                                "#7D1D67")) 
#
#theme(legend.position="none")
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
  scale_shape_manual(values = c(16,
                                16)) +
  scale_colour_manual(values = c("#A45EE9")) +
  #viridis::scale_color_viridis(discrete = TRUE) +
  theme(legend.position="none") +
  scale_y_continuous(limits= c(-0.2, 0.27))


# scale_colour_manual(values = c( 
#                                "darkgrey",
#                                "#7D1D67")) 
#
#theme(legend.position="none")

##### ACCUMBENS
plot <- newdf


####
plot2 <- plot


### Make colours be related to protein + significance
plot2$concat <- paste(plot2$significance, plot2$DNAm)


plot2 %<>% mutate(concat =
                    case_when(str_detect(concat, "No") ~ "not significant",
                              TRUE ~ concat))
###### REORDER DNAm by estimate hippocampus

plot2 %<>% filter(brain_metric == "accumbens" & significance == "Yes")

# How many unique colours are needed?
pal <- colorspace::sequential_hcl(25, palette = "SunsetDark")

x <- ggplot(
  plot2,
  aes(
    x = reorder(DNAm,-estimate),
    # x = DNAm,
    y = estimate,
    #colour = reorder(DNAm, estimate),
    colour = brain_metric,
    alpha = reorder(DNAm, -estimate),
    group = DNAm,
    shape = significance
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
  scale_shape_manual(values = c(16,
                                16)) +
  scale_colour_manual(values = c("#A45EE9")) +
  #viridis::scale_color_viridis(discrete = TRUE) +
  theme(legend.position="none") +
  scale_y_continuous(limits= c(-0.2, 0.27))


# scale_colour_manual(values = c( 
#                                "darkgrey",
#                                "#7D1D67")) 
#
#theme(legend.position="none")
####
#### PUTAMEN

plot2 <- plot


### Make colours be related to protein + significance
plot2$concat <- paste(plot2$significance, plot2$DNAm)


plot2 %<>% mutate(concat =
                    case_when(str_detect(concat, "No") ~ "not significant",
                              TRUE ~ concat))
###### REORDER DNAm by estimate hippocampus

plot2 %<>% filter(brain_metric == "putamen" & significance == "Yes")

# How many unique colours are needed?
pal <- colorspace::sequential_hcl(25, palette = "SunsetDark")

x <- ggplot(
  plot2,
  aes(
    x = reorder(DNAm,-estimate),
    # x = DNAm,
    y = estimate,
    #colour = reorder(DNAm, estimate),
    colour = brain_metric,
    alpha = reorder(DNAm, -estimate),
    group = DNAm,
    shape = significance
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
  scale_shape_manual(values = c(16,
                                16)) +
  scale_colour_manual(values = c("#F88379")) +
  #viridis::scale_color_viridis(discrete = TRUE) +
  theme(legend.position="none") +
  scale_y_continuous(limits= c(-0.2, 0.2))


# scale_colour_manual(values = c( 
#                                "darkgrey",
#                                "#7D1D67")) 
#
#theme(legend.position="none")


#### CEREBELLUM WM
plot <- newdf


####
plot2 <- plot


### Make colours be related to protein + significance
plot2$concat <- paste(plot2$significance, plot2$DNAm)


plot2 %<>% mutate(concat =
                    case_when(str_detect(concat, "No") ~ "not significant",
                              TRUE ~ concat))
###### REORDER DNAm by estimate hippocampus

plot2 %<>% filter(brain_metric == "cerebellum (WM)" & significance == "Yes")

# How many unique colours are needed?
pal <- colorspace::sequential_hcl(25, palette = "SunsetDark")

x <- ggplot(
  plot2,
  aes(
    x = reorder(DNAm,-estimate),
    # x = DNAm,
    y = estimate,
    #colour = reorder(DNAm, estimate),
    colour = brain_metric,
    alpha = reorder(DNAm, -estimate),
    group = DNAm,
    shape = significance
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
  scale_shape_manual(values = c(16,
                                16)) +
  scale_colour_manual(values = c("#69C8A1")) +
  #viridis::scale_color_viridis(discrete = TRUE) +
  theme(legend.position="none") +
  scale_y_continuous(limits= c(-0.2, 0.2))


# scale_colour_manual(values = c( 
#                                "darkgrey",
#                                "#7D1D67")) 
#
#theme(legend.position="none")


# 2.59