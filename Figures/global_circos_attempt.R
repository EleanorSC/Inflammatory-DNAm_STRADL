
# Brain_age, brain_accel
FULL_cognitive_list <- c(
  ### Global measures
  "Brain_age", 
  "brain_accel"
)

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

# Making a data frame with all combinations of variables
df <-
  as.data.frame(expand.grid(FULL_DNAm_list, FULL_cognitive_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1, metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
get_summary_values <- function(list) {
  metric <- list$metric
  DNAm <- list$DNAm
  tib <-
    broom::tidy(summary(
      lm(
        scale(PROTEOMICS_DNAm_DATA_COGNITION[[metric]]) ~
          scale(st_age)
        + sex
        + scale(PROTEOMICS_DNAm_DATA_COGNITION[[DNAm]]),
        data = PROTEOMICS_DNAm_DATA_COGNITION
      )
    ))[4, c(2, 3, 5)]
  return(tib)
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = get_summary_values)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality = "brain_ageing") %>%
  
  # Create a brain metric column
  mutate(
    brain_metric =
      case_when(
        
        metric == "Brain_age" ~ "relative brain age",
        metric == "brain_accel" ~ "brain acceleration",
        
        
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
                     TRUE ~ "No"))%>%
  
  # Create a column to denote whether these are DNAm or proteins
  mutate(omic_type = "DNAm")


plot_brain_age_methylation <- newdf


##########
# COGNTIIVIE
#########
FULL_cognitive_list <- c(
  ### Global measures
  "g", "gf", 
  #"g_STRADL",
  ### individual domains
  "processing_speed", 
  "executive_function",
  "verbal_declarative_memory"
)

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

# Making a data frame with all combinations of variables
df <-
  as.data.frame(expand.grid(FULL_DNAm_list, FULL_cognitive_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1, metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
get_summary_values <- function(list) {
  metric <- list$metric
  DNAm <- list$DNAm
  tib <-
    broom::tidy(summary(
      lm(
        scale(PROTEOMICS_DNAm_DATA_COGNITION[[metric]]) ~
          scale(st_age)
        + sex
        + scale(PROTEOMICS_DNAm_DATA_COGNITION[[DNAm]]),
        data = PROTEOMICS_DNAm_DATA_COGNITION
      )
    ))[4, c(2, 3, 5)]
  return(tib)
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = get_summary_values)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality = "cognition") %>%
  
  # Create a brain metric column
  mutate(
    brain_metric =
      case_when(
        
        metric == "g" ~ "g",
        metric == "gf" ~ "gf",
        # metric == "g_STRADL" ~ "g_STRADL",
        metric == "processing_speed" ~ "processing speed",
        metric == "executive_function" ~ "executive function",
        metric == "verbal_declarative_memory" ~ "verbal declarative memory",
        
        TRUE ~ "misc"
      )
  ) %>%
  
  # Create pFDR column
  # group_by(metric, Hemisphere) %>%
  group_by(DNAm) %>%
  mutate(pFDR = p.adjust(p.value, method = "fdr")) %>%
  # Create a column to denote where there are significant hits
  mutate(significance =
           case_when(p.value < 0.05 ~ "Yes",
                     TRUE ~ "No")) %>%
  # Create a column to denote where there are FDR significant hits
  mutate(FDR_significance =
           case_when(pFDR < 0.05 ~ "Yes",
                     TRUE ~ "No"))%>%
  
  # Create a column to denote whether these are DNAm or proteins
  mutate(omic_type = "DNAm")


plot_cognitive_methylation <- newdf
########

######Load up neuroimaging full dataset
STRADL_FreeSurfer <- read.csv("STRADL_Measures_FreeSurfer_Main.csv")
# Now link these to STRADL ID
names(STRADL_FreeSurfer)[names(STRADL_FreeSurfer) == 'id'] <- 'stradl_ID'


# n =709
Neuroimaging_DNAm <- merge(PROTEOMICS_DNAm_DATA, 
                           STRADL_FreeSurfer,
                           by = "stradl_ID")



Neuroimaging_DNAm %<>% mutate(global_cortical_surface_area = 
                                          hem.lh.csa + hem.rh.csa,
                                        global_cortical_thickness = 
                                          hem.lh.ct + hem.rh.ct,
                                        global_cortical_volume = 
                                          hem.lh.cv + hem.rh.cv) 


skimr::skim(Neuroimaging_DNAm)

FULL_neuroimaging_list <- c(
  ### Global measures
  "global.cerebral.wm",
  "global.total.gm",
  "global.wbv",
  # "global.wbv.vent",
  "global_cortical_volume",
  "global_cortical_thickness",
  "global_cortical_surface_area",
  "gFA",
  "gMD"
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
             str_detect(metric, "global") ~ "global",
             TRUE ~ "misc")
  ) %>%
  
  # Create a brain metric column
  mutate(brain_metric =
           case_when(
             
             ### Global
             metric == "global.cerebral.wm" ~ "global white matter",
             metric == "global.total.gm" ~ "global grey matter",
             metric == "global.wbv" ~ "total brain volume",
             #metric == "global.wbv.vent" ~ "global ventricles volume",
             metric == "global_cortical_volume" ~ "global cortical volume",
             metric == "global_cortical_thickness" ~ "global cortical thickness",
             metric == "global_cortical_surface_area" ~ "global cortical surface area",
             
             metric == "gFA" ~ "gFA",
             metric == "gMD" ~ "gMD",
             
             TRUE ~ "misc")
  ) %>%
  
  # Create pFDR column
  # group_by(metric, Hemisphere) %>%
  group_by(DNAm) %>%
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
plot_neuroimaging_methylation <- newdf











##### REORDER cognitive by global then domains
plot2 <- rbind(plot_brain_age_methylation,
               plot_neuroimaging_methylation,
               plot_cognitive_methylation)
               


plot2 %<>% filter(significance == "Yes")


#plot2 %<>% mutate(brain_metric = factor(
#  brain_metric,
#  levels = c(
#    "relative brain age", 
#    "brain acceleration"
#  )))
#
##### for effect sizes plot
facetSettings <-
  theme(strip.background = element_rect(
    #fill = "#F8EEEC",
    # fill = "#F3FFF9",
    # fill = "#E5E1DF", #light brown
    fill = "#FBE4E3",
    colour = "black",
    size = 1
  ))

ggplot(plot2,
       
       aes(x = reorder(DNAm,(estimate)),
           y = estimate,
           alpha = reorder(DNAm,
                           (estimate)),
           shape = FDR_significance,
           col = brain_metric,
           group = omic_type,
           # alpha = significance
       )
) +
  
  geom_point(position = position_dodge(width = 0.9),
             size = 1.9,
             #colour = "#414487FF",
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
                               option = "H") +
  
  scale_shape_manual(values = c(1,
                                16)
  ) +
  facet_wrap(~ brain_metric,
             #scales="free",
             #scales = "free_y",
             nrow = 1) +
  facetSettings 


###### CHord plot

##### REORDER cognitive by global then domains
plot2 <- rbind(plot_brain_age_methylation,
               plot_neuroimaging_methylation,
               plot_cognitive_methylation)


plot2 %<>% filter(brain_metric == "g" |
                    brain_metric == "gf"|
                    
                    brain_metric == "relative brain age"|

                    brain_metric == "total brain volume"|
                    brain_metric == "global grey matter"|
                    brain_metric == "global white matter"|
                    brain_metric == "global cortical volume"|
                    brain_metric == "global cortical thickness"|
                    brain_metric == "global cortical surface area"|
                    
                    brain_metric == "gMD"|
                    brain_metric == "gFA"
                    )

plot2 %<>% filter(significance == "Yes")


####
plot3 <- plot2 %>% select(DNAm, brain_metric, estimate)

plot3 %<>% 
  dplyr::rename(
    rowname = DNAm,
    key = brain_metric,
    value = estimate)

plot3$rowname <- as.character(plot3$rowname)

#install.packages("circlize")
#library(circlize)
# parameters


circos.clear()
circos.par(
  
 #gap.after = c(rep(5, nrow(plot3)-1), 
 #              15, 
 #              rep(5, ncol(plot3)-1),
 #              15),
  
  circos.par(gap.after = c("g" = 2,
                           "gf" = 1,
                           "relative brain age"= 1,
                           "total brain volume"= 1,
                           "global grey matter"= 1,
                           "global white matter"= 1,
                           "global cortical volume"= 1,
                           "global cortical thickness"= 1,
                           "global cortical surface area"= 1,
                           "gMD"= 1,
                           "gFA"= 2,
                           "MMP12"        = 0.1 ,"PIGR"        = 0.1 ,"SEMA3E"      = 0.1 ,"NCAM1"      = 0.1 , "NTRK3" = 0.1      , "SERPIND1"    =0.1,
                           "NOTCH1"       = 0.1 ,"CRP"         = 0.1 ,"SKR3"        = 0.1 ,"IGFBP4"     = 0.1 , "VEGFA" = 0.1      , "CNTN4"       =0.1,
                           "PRSS2"        = 0.1 ,"FGF.21"      = 0.1 ,"CCL18"       = 0.1 ,"SLITRK5"    = 0.1 , "CCL17" = 0.1      , "RARRES2"     =0.1,
                           "ICAM5"        = 0.1 ,"NTRK3_olink" = 0.1 ,"THBS2"       = 0.1 ,"SELL"       = 0.1 , "OMD"   = 0.1      , "SELE"        =0.1,
                           "MMP.1_olink"  = 0.1 ,"WFIKKN2"     = 0.1 ,"GDF.8"       = 0.1 ,"SIGLEC1"    = 0.1 , "MMP1"  = 0.1      , "TPSB2"       =0.1,
                           "MST1"         = 0.1 ,"HGF"         = 0.1 ,"TGF.alpha"   = 0.1 ,"ACY1"       = 0.1 , "SHBG"  = 0.1      , "MMP9"        =0.1,
                           "B2M"          = 0.1 ,"ENPP7"       = 0.1 ,"ADAMTS13"    = 0.1 ,"AFM"        = 0.1 , "G.CSF" = 0.1      , "CXCL9"       =0.1,
                           "GP1BA"        = 0.1 ,"NEP"         = 0.1 ,"CXCL10_olink"= 0.1 ,"BCAM"       = 0.1 , "CCL25" = 0.1      , "ESM1"        =0.1,
                           "MMP2"         = 0.1 ,"CD209"       = 0.1 ,"LTA|LTB"     = 0.1 ,"MRC2"       = 0.1 , "FAP"   = 0.1      , "GZMA_olink"  =0.1,
                           "IGFBP1"       = 0.1 ,"SPOCK2"      = 0.1 ,"CCL11"       = 0.1 ,"N.CDase"    = 0.1 , "CXCL11"= 0.1      , "CXCL11_olink"=0.1,
                           "LGALS3BP"     = 0.1 ,"TNFRSF17"    = 0.1 ,"EDA"         = 0.1 ,"MIA"        = 0.1 , "CCL22" = 0.1      , "HGFAC"       =0.1,
                           "MPL"          = 0.1 ,"CD5L"        = 0.1 ,"IDUA"        = 0.1 ,"GNLY"       = 0.1 , "FCER2" = 2
                           
                           )),
  
  start.degree = 90, 
  gap.degree = 0.1, 
  track.margin = c(-0.1, 0.1), 
  points.overflow.warning = FALSE
)
#par(mar = rep(0, 4))
#par(mfrow = c(1))
#




#install.packages("colorspace")
#pal <- colorspace::sequential_hcl(71, palette = "SunsetDark")
#
#pal2 <- colorspace::sequential_hcl(11, palette = "Purple-Blu")

pal <- c(
  "#7D1D67","#811D68", "#851E6A" ,"#891E6C" ,"#8D1E6D", "#911E6F", "#951F70", "#991F72", "#9D2073",
  "#A12074","#A52175", "#A92276" ,"#AD2377" ,"#B12477", "#B52578", "#B92779", "#BC2879", "#C02A79",
  "#C42B79","#C72D7A", "#CB2F7A" ,"#CF3179" ,"#D23379", "#D63579", "#D93778", "#DD3978", "#E03C77",
  "#E43E76","#E74075", "#EA4374" ,"#ED4572" ,"#F14871", "#F44B6F", "#F5506E", "#F6556E", "#F65A6D",
  "#F75F6D","#F8646C", "#F8686C" ,"#F96D6C" ,"#FA716C", "#FB756C", "#FB796C", "#FC7D6C", "#FD816D",
  "#FD856D","#FE886E", "#FF8C6F" ,"#FF9070" ,"#FF9471", "#FF9772", "#FF9B74", "#FF9E75", "#FFA277",
  "#FFA579","#FFA97A", "#FFAC7C" ,"#FFAF7F" ,"#FFB381", "#FFB683", "#FFB985", "#FFBD88", "#FFC08A",
  "#FFC38D","#FFC68F", "#FFC992" ,"#FFCD95" ,"#FFD097", "#FFD39A", "#FFD69D", "#FFD99F",
  "#6B0077","#6E3889","#73579B","#7A72AC","#828BBC","#8DA3CA","#9CB8D6","#ADCCE0","#C1DDE8","#D7EAEF","#F1F1F1")


col_fun = function(x) ifelse(x < -0.07 | x > 0.07, pal, #
                             "#ECECEC"
                            # "#00000000"
                             )

transp_fun = function(x) ifelse(x < -0.07 | x > 0.07, 
                                0.25, 
                                0.01)



### arrange circos by size 
### Which DNAm proxy has the most significant hits? 
# Group by sum using dplyr
test <- aggregate(abs(plot3$value), 
          by=list(rowname=plot3$rowname), 
          FUN=sum) %>% 
  arrange(desc(x))

test$rowname


# Base plot
chordDiagram(
  x = plot3, 
  #big.gap = 30,
  order = c(
    
     "MMP12"        ,"PIGR"         ,"SEMA3E"       ,"NCAM1"       , "NTRK3"       , "SERPIND1"    ,
     "NOTCH1"       ,"CRP"          ,"SKR3"         ,"IGFBP4"      , "VEGFA"       , "CNTN4"       ,
     "PRSS2"        ,"FGF.21"       ,"CCL18"        ,"SLITRK5"     , "CCL17"       , "RARRES2"     ,
     "ICAM5"        ,"NTRK3_olink"  ,"THBS2"        ,"SELL"        , "OMD"         , "SELE"        ,
     "MMP.1_olink"  ,"WFIKKN2"      ,"GDF.8"        ,"SIGLEC1"     , "MMP1"        , "TPSB2"       ,
     "MST1"         ,"HGF"          ,"TGF.alpha"    ,"ACY1"        , "SHBG"        , "MMP9"        ,
     "B2M"          ,"ENPP7"        ,"ADAMTS13"     ,"AFM"         , "G.CSF"       , "CXCL9"       ,
     "GP1BA"        ,"NEP"          ,"CXCL10_olink" ,"BCAM"        , "CCL25"       , "ESM1"        ,
     "MMP2"         ,"CD209"        ,"LTA|LTB"      ,"MRC2"        , "FAP"         , "GZMA_olink"  ,
     "IGFBP1"       ,"SPOCK2"       ,"CCL11"        ,"N.CDase"     , "CXCL11"      , "CXCL11_olink",
     "LGALS3BP"     ,"TNFRSF17"     ,"EDA"          ,"MIA"         , "CCL22"       , "HGFAC"       ,
     "MPL"          ,"CD5L"         ,"IDUA"         ,"GNLY"        , "FCER2" ,
            
              "g" ,
              "gf",
              
              "relative brain age",
              
              "total brain volume",
              "global grey matter",
              "global white matter",
              "global cortical volume",
              "global cortical thickness",
              "global cortical surface area",
              
              "gMD",
              "gFA"
             ),
  #col = col_fun,
  grid.col = pal,
  #transparency = transp_fun,
  transparency = 0.25,
  symmetric = TRUE,
  #directional = 1,
  #direction.type = c("arrows", "diffHeight"), 
  #diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.decreasing = TRUE
  #link.largest.ontop = TRUE
)

# Add text and axis
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
#    # Add names to the sector. 
#    circos.text(
#      x = mean(xlim), 
#      y = 4.5, 
#      labels = sector.index, 
#      #adj = c(0, degree(5)), 
#      niceFacing = TRUE,
#      facing = "clockwise",
#      cex = 0.5
#    )
    
    circos.rect(xleft=xlim[1], 
                ybottom=0.3, 
                xright=xlim[2], 
                ytop=0.32, 
                col = "white", 
                border = "white")
    
    
    
    # Attempt to separate labels 
    #circos.axis(h = "top", 
    #            labels.cex = 0.5, 
    #            major.tick.percentage = 0.2, 
    #            sector.index = get.cell.meta.data("sector.index")
    #            
    #            #track.index = 2
    #            )
    #
  }
)
circos.clear()