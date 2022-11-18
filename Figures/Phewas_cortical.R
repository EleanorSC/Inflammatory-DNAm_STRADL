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
  "SPOCK2")




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

#
#FULL_DNAm_list <- c("CRP", "IGFBP4","PIGR")
#
#FULL_neuroimaging_list <- c("cv.bilat.supramarginal",
#                            "cv.bilat.temporalpole")

# Converting multiple varibles into a factor
Neuroimaging_DNAm %<>% mutate_at(c("sex", "site", "edited", "batch"),
                                 as.factor)

# Making a data frame with all combinations of variables
df <- as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1, 
                metric = Var2)

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


### Define whether markers are good or bad for brain health

plot2 %<>% mutate(direction =
                    as.factor(case_when(estimate < 0 ~ "Bad",
                              TRUE ~ "Good")))


#### Only plot instances where DNAm signature has some significant associations

plot2 %<>% filter(DNAm == "ADAMS13"|
                    DNAm == "ADIPOQ"|
                    DNAm == "BCAM"|
                    DNAm == "BMP1"|
                    DNAm == "C5"|
                    DNAm == "C9"|
                    DNAm == "CCL17"|
                    DNAm == "CCL19"|
                    DNAm == "CCL22"|
                    DNAm == "CCL25"|
                    DNAm == "CNTN4"|
                    DNAm == "CRP"|
                    DNAm == "ENPP7"|
                    DNAm == "FAP"|
                    DNAm == "FCER2"|
                    DNAm == "FGF.21"|
                    DNAm == "G.CSF"|
                    DNAm == "GDF.8"|
                    DNAm == "GNLY"|
                    DNAm == "GZMA_olink"|
                    DNAm == "HGF"|
                    DNAm == "HGFAC"|
                    DNAm == "ICAM5"|
                    DNAm == "IGFBP4"|
                    DNAm == "IL19"|
                    DNAm == "LGALS3BP"|
                    DNAm == "LTF"|
                    DNAm == "MMP.1_olink"|
                    DNAm == "MMP1"|
                    DNAm == "MMP12"|
                    DNAm == "MMP2"|
                    DNAm == "MMP9"|
                    DNAm == "MPL"|
                    DNAm == "MPO"|
                    DNAm == "MRC2"|
                    DNAm == "MST1"|
                    DNAm == "NCAM1"|
                    DNAm == "NEP"|
                    DNAm == "NOTCH1"|
                    DNAm == "NTRK3"|
                    DNAm == "NTRK3_olink"|
                    DNAm == "OMD"|
                    DNAm == "PIGR"|
                    DNAm == "PRSS2"|
                    DNAm == "RARRES2"|
                    DNAm == "RETN"|
                    DNAm == "SELL"|
                    DNAm == "SEMA3E"|
                    DNAm == "SERPIND1"|
                    DNAm == "SKR3"|
                    DNAm == "SLITRK5"|
                    DNAm == "STC1"|
                    DNAm == "THBS2"|
                    DNAm == "VEGFA" 
                    )

# Only plot associations that associate with poor brain health outcomes

#plot2 %<>% filter(direction == "Bad")
#### Only plot instances where DNAm signature has some significant associations

plot2 %<>% filter(
                   #DNAm == "ADAMS13"|
                   #DNAm == "ADIPOQ"|
                    DNAm == "BCAM"|
                   #DNAm == "BMP1"|
                    DNAm == "C5"|
                    DNAm == "C9"|
                    DNAm == "CCL17"|
                    DNAm == "CCL19"|
                    DNAm == "CCL22"|
                    #DNAm == "CCL25"|
                    DNAm == "CNTN4"|
                    DNAm == "CRP"|
                    DNAm == "ENPP7"|
                    #DNAm == "FAP"|
                    #DNAm == "FCER2"|
                    DNAm == "FGF.21"|
                    DNAm == "G.CSF"|
                    #DNAm == "GDF.8"|
                    #DNAm == "GNLY"|
                    #DNAm == "GZMA_olink"|
                    DNAm == "HGF"|
                    DNAm == "HGFAC"|
                    DNAm == "ICAM5"|
                    DNAm == "IGFBP4"|
                    #DNAm == "IL19"|
                    DNAm == "LGALS3BP"|
                    DNAm == "LTF"|
                    DNAm == "MMP.1_olink"|
                    DNAm == "MMP1"|
                    DNAm == "MMP12"|
                    #DNAm == "MMP2"|
                    DNAm == "MMP9"|
                    #DNAm == "MPL"|
                    DNAm == "MPO"|
                    #DNAm == "MRC2"|
                    DNAm == "MST1"|
                    #DNAm == "NCAM1"|
                    DNAm == "NEP"|
                    #DNAm == "NOTCH1"|
                    #DNAm == "NTRK3"|
                    #DNAm == "NTRK3_olink"|
                    #DNAm == "OMD"|
                    DNAm == "PIGR"|
                    DNAm == "PRSS2"|
                    DNAm == "RARRES2"|
                    #DNAm == "RETN"|
                    #DNAm == "SELL"|
                    #DNAm == "SEMA3E"|
                    DNAm == "SERPIND1"|
                    DNAm == "SKR3"|
                    #DNAm == "SLITRK5"|
                    DNAm == "STC1"|
                    DNAm == "THBS2"|
                    DNAm == "VEGFA" 
)

##### PHEWAS PLOT
facetSettings <-
  theme(strip.background = element_rect(
    fill = "#EAE3F2", #purple
    colour = "black",
    size = 1
  ))


viridis::rocket(n = 34)


pals <- c(
 "#03051AFF", "#0D0A21FF", "#170F28FF", "#221331FF" ,"#2E1739FF","#391A41FF",
 "#451C47FF", "#511E4DFF", "#5E1F52FF", "#6A1F56FF" ,"#761F58FF","#841E5AFF",
 "#921C5BFF", "#9E1A5BFF", "#AB185AFF", "#B91657FF" ,"#C51852FF","#D01E4DFF",
 "#D92847FF", "#E23442FF", "#E8413EFF", "#ED513EFF" ,"#F06043FF","#F26E4CFF",
 "#F47C56FF", "#F58A61FF", "#F5976EFF", "#F6A47BFF" ,"#F6B08AFF","#F6BC99FF",
 "#03051AFF", "#0D0A21FF", "#170F28FF", "#221331FF" 
# "#F7C9AAFF", "#F8D4BBFF", "#F9DFCBFF", "#FAEBDDFF"
)


ggplot(plot2,
       aes(
         x = reorder(DNAm, estimate),
         y = -log(p.value),
         shape = FDR_significance,
         colour = brain_metric
         # alpha = brain_metric,
         #colour = DNAm
       )) +
  
  geom_point(aes(#col = DNAm
    alpha = brain_metric),
    #size = -(estimate)
    size = 1.2,
    alpha = 0.8) +
  
  theme_classic() +
  
  theme(
    strip.text = element_text(
      size = 6,
      face = "bold",
      family = "sans",
      colour = "black"
    ),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.grid.minor = element_line(colour = "grey", linetype = "dashed"),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) +
  
  labs(color = "Category",
       #size = "Effect size",
       x = "",
       y = "") +
  
  ggrepel::geom_text_repel(
    data = . %>% mutate(label = ifelse(p.value < 0.05, as.character(DNAm), "")),
    aes(label = label),
    direction = "y",
    size = 2,
    box.padding = unit(0.7, "lines"),
    max.overlaps = Inf
  ) +
  
  
  geom_hline(
    yintercept = -log(0.05),
    color = "darkgrey",
    size = 1,
    alpha = 0.5
  ) +
  
  facet_wrap( ~ brain_metric,
              scales = "free_x") +
  
  scale_shape_manual(values = c(16,
                                8)) +
  
  scale_colour_manual(values = pals)

# viridis::scale_color_viridis(discrete = TRUE,
#                              option = "F")




####### effect sizes plot
# ----------------------------# 
# PLOT 
# ----------------------------#



plot2 <- plot_neuroimaging_methylation

plot2 %<>% filter(significance == "Yes")


# Create a column to denote where there are FDR significant hits, significant, or not
#plot2 %<>% mutate(SIG =
#         case_when(FDR_significance == "Yes" ~ "pFDR",
#                   significance == "Yes" ~ "p",
#                   TRUE ~ "not_significant"))



# ----------------------------# 
# Reorder in terms of most significant hits 
# ----------------------------#
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


#####
###### REORDER DNAm by no. sig
plot2 %<>% mutate(brain_metric = factor(
  brain_metric,
  levels = c(
                   "parahippocampal" , # n = 43
                          "fusiform" , # n = 39
                    "supra marginal" , # n = 31
                 "inferior temporal" , # n = 27
                 "superior parietal" , # n = 27
                        "precentral" , # n = 25
                 "superior temporal" , # n = 25
                 "lateral occipital" , # n = 24
             "lateral orbitofrontal" , # n = 22
                            "insula" , # n = 20
                       "postcentral" , # n = 20
        "rostral anterior cingulate" , # n = 20
                           "lingual" , # n = 19
                        "entorhinal" , # n = 18
                   "middle temporal" , # n = 18
                         "precuneus" , # n = 18
            "rostral middle frontal" , # n = 18
               "transverse temporal" , # n = 18
                 "inferior parietal" , # n = 16
                    "pars orbitalis" , # n = 16
               "posterior cingulate" , # n = 16
 "banks of superior temporal sulcus" , # n = 15
              "medial orbitofrontal" , # n = 14
                      "frontal pole" , # n = 9
             "caudal middle frontal" , # n = 7
                       "paracentral" , # n = 7
                  "pars opercularis" , # n = 7
                  "superior frontal" , # n = 6
                            "cuneus" , # n = 5
                 "pars triangularis" , # n = 4
                     "pericalcarine" , # n = 4
                     "temporal pole" , # n = 4
         "caudal anterior cingulate" , # n = 3
                 "isthmus cingulate"  # n = 3
  ))) 




facetSettings <-
  theme(strip.background = element_rect(
    fill = "#F8EEEC",
    colour = "black",
    size = 1
  ))

ggplot(plot2,
       
       aes(
         x = reorder(DNAm, -(estimate)),
         y = estimate,
         alpha = significance

       )) +
  
  geom_point(
    aes(
      col = reorder(DNAm, (-estimate)),
      alpha = significance,
      #alpha = reorder(metric, -(estimate)),
      shape = FDR_significance 
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
  
  theme(legend.position = "none",
        
    strip.text = element_text(
      size = 6,
      face = "bold",
      family = "sans",
      colour = "black"
    ),
    axis.text.x = element_text(
      #  angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 6,
      #face = "bold",
      family = "sans"
    ),
    axis.text.y = element_text(size = 6)
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
                               option = "A",
                               direction = -1) +
  
  
facet_wrap( ~ brain_metric,
            #scales="free"
            scales = "free_y",
            nrow = 5
           # nrow = 1
            ) +
  facetSettings +
  scale_alpha_manual(values=c(0.8, 1, 0.2)) +
  scale_shape_manual(values = c(
    # 8,
    # 1,
  #  1,
    16,
    5))


# ----------------------------# 
# PLOT 
# ----------------------------#
plot2 <- plot_neuroimaging_methylation

plot2 %<>% filter(FDR_significance == "Yes")



### Which DNAm proxy has the most FDR significant hits?
plot2 %<>%
  mutate(number_significant =
           case_when(FDR_significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
top_significant_hits <- aggregate(plot2$number_significant,
                                  by = list(brain_metric = plot2$brain_metric),
                                  FUN = sum)

top_significant_hits %<>% arrange(desc(x))


#####
###### REORDER DNAm by no. sig
plot2 %<>% mutate(brain_metric = factor(
  brain_metric,
  levels = c(
               "fusiform", #n = 15
         "supra marginal", #n = 13
      "superior parietal", #n = 10
        "parahippocampal", #n =  8
             "precentral", #n =  8
      "inferior temporal", #n =  6
              "precuneus", #n =  6
      "superior temporal", #n =  6
                 "insula", #n =  5
 "rostral middle frontal", #n =  4
             "entorhinal", #n =  3
  "lateral orbitofrontal", #n =  3
                "lingual" #n =  1
  )))




  



facetSettings <-
  theme(strip.background = element_rect(
    fill = "#F8EEEC",
    colour = "black",
    size = 1
  ))

ggplot(plot2,
       
       aes(
         x = reorder(DNAm, -(estimate)),
         y = estimate
       )) +
  
  geom_point(
    aes(
      col = reorder(DNAm, (-estimate)),
      alpha = significance,
      shape = FDR_significance 
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
  
  theme(legend.position = "none",
        
        strip.text = element_text(
          size = 6,
          face = "bold",
          family = "sans",
          colour = "black"
        ),
        axis.text.x = element_text(
          #  angle = 90,
          vjust = 0.5,
          hjust = 1,
          size = 6,
          #face = "bold",
          family = "sans"
        ),
        axis.text.y = element_text(size = 6)
        
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
                               option = "A",
                               direction = -1) +
  
  
  facet_wrap( ~ brain_metric,
              
              #scales = "free",
              nrow = 2) +
  
  facetSettings +
  scale_alpha_manual(values=c(0.8, 1, 0.2)) +
  scale_shape_manual(values = c(5))


# ----------------------------# 
# PLOT - barplot of significant and FDR significant 
# ----------------------------#
plot2 <- plot_neuroimaging_methylation


### Which DNAm proxy has the most FDR significant hits?
plot2 %<>%
  mutate(number_significant =
           case_when(significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
top_significant_hits <- aggregate(plot2$number_significant,
                                  by = list(brain_metric = plot2$brain_metric),
                                  FUN = sum)

top_significant_hits %<>% arrange(desc(x)) %>% rename(p = x)

plot2 %<>%
  mutate(number_significant_FDR =
           case_when(FDR_significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
top_significant_hits_2 <- aggregate(plot2$number_significant_FDR,
                                  by = list(brain_metric = plot2$brain_metric),
                                  FUN = sum) %>% rename(pFDR = x)


top_significant_hits <- merge(top_significant_hits,
                              top_significant_hits_2, 
                              by = "brain_metric")

######



ggplot(top_significant_hits,
       aes(y = reorder(brain_metric, x),
           x = x)) +
  geom_bar(stat="identity")
