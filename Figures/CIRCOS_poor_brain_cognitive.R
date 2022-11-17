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


ROTEOMICS_DNAm_DATA <- read.csv("PROTEOMICS_DNAm_DATA.csv")

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
#%>%
#  mutate(model = "Model 1 (n=709)")

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

##########
# COGNTIIVIE
#########
FULL_cognitive_list <- c(
  
  ### Global measures
  "g", 
  "gf", 
  
  ### individual domains
  "processing_speed", 
  "executive_function",
  "verbal_declarative_memory"
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
        metric == "brain_accel" ~ "relative brain age",
        
        metric == "processing_speed" ~ "processing speed",
        metric == "executive_function" ~ "executive function",
        metric == "verbal_declarative_memory" ~ "verbal declarative memory",
        
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


plot_cognitive_methylation <- newdf
##########


plot2 <- rbind(plot_neuroimaging_methylation, plot_cognitive_methylation)

##### For circos, certain associations that associate with POOR brain health neeed to be recoded
plot2 <- plot2 %>% mutate(estimate =
                            case_when(brain_metric == "relative brain age" ~ estimate*-1,
                                      brain_metric == "gMD" ~ estimate*-1,
                                      brain_metric == "WMH" ~ estimate*-1,
                                      TRUE ~ estimate)
                          )

plot2 %<>% filter(significance == "Yes" & estimate < 0)

####
plot3 <- plot2 %>% select(DNAm, brain_metric, estimate)


plot3 %<>% filter(brain_metric == "g" |
                    brain_metric == "gf" |
                    brain_metric == "processing speed" |
                    brain_metric == "relative brain age" |
                    brain_metric == "global subcortical volume" |
                    brain_metric == "global cortical volume" |
                    brain_metric == "WMH" |
                    brain_metric == "gFA" |
                    brain_metric == "gMD" |
                    brain_metric == "global grey matter" |
                    brain_metric == "global white matter" |
                    brain_metric == "total brain volume" )


plot3 %<>% 
  dplyr::rename(
    rowname = DNAm,
    key = brain_metric,
    value = estimate)

plot3$rowname <- as.character(plot3$rowname)

### arrange circos by size 
### Which DNAm proxy has the most significant hits? 
# Group by sum using dplyr
test <- aggregate((plot3$value), 
                  by=list(rowname=plot3$rowname), 
                  FUN=sum) %>% 
  arrange((x))
#arrange(desc(x))

test$rowname

test$rowname <- as.factor(test$rowname)
skimr::skim(test$rowname)
#n=44

#install.packages("circlize")
#library(circlize)
# parameters


circos.clear()


circos.par(gap.after = c("g" = 1.5,
                         "gf" = 1,
                         "processing speed" = 1,
                         "relative brain age"= 1,
                         "WMH" = 1,
                         "total brain volume"= 1,
                         "global grey matter"= 1,
                         "global white matter"= 1,
                         "global cortical volume"= 1,
                         "global subcortical volume"= 1,
                         "gMD"= 1,
                         "gFA"= 1.5,
                         
                          "MMP12"  =0.1  ,   "PIGR"    =0.1 ,   "CRP"         =0.1,"IGFBP4"    =0.1,  "SERPIND1" =0.1 ,  "SKR3"     =0.1, 
                          "VEGFA"  =0.1  ,   "FGF.21"  =0.1 ,   "RARRES2"     =0.1,"CCL18"     =0.1,  "PRSS2"    =0.1 ,  "THBS2"    =0.1, 
                          "CCL17"  =0.1  ,   "ICAM5"   =0.1 ,   "MMP.1_olink" =0.1,"MMP9"      =0.1,  "MMP1"     =0.1 ,  "HGF"      =0.1, 
                          "B2M"    =0.1  ,   "AFM"     =0.1 ,   "SIGLEC1"     =0.1,"TGF.alpha" =0.1,  "SELE"     =0.1 ,  "CXCL9"    =0.1, 
                          "CXCL11" =0.1  ,   "CCL11"   =0.1 ,   "ENPP7"       =0.1,"NEP"       =0.1,  "CCL25"    =0.1 ,  "ACY1"     =0.1, 
                          "MST1"   =0.1  ,   "EN.RAGE" =0.1 ,   "G.CSF"       =0.1,"S100A9"    =0.1,  "RETN"     =0.1 ,  "CLEC11A.1"=0.1, 
                          "FCGR3B" =0.1  ,   "OSM"     =0.1 ,   "TNFRSF17"    =0.1,"IDUA"      =0.1,  "C4A|C4B"  =0.1 ,  "CCL22"    =0.1, 
                          "HGFAC"  =0.1  ,   "CD5L" =0.1
                         
                         
),

start.degree = 90, 
gap.degree = 0.1, 
track.margin = c(-0.1, 0.1), 
points.overflow.warning = FALSE
)
#par(mar = rep(0, 4))
#par(mfrow = c(1))
#




#install.packages("colorspace")


plasma_pal <- c(
  colorspace::sequential_hcl(44, palette = "SunsetDark"),
  #viridis::mako(n = 11)
  colorspace::sequential_hcl(11, palette = "Purple-Blu")
  )


#pal2 <- colorspace::sequential_hcl(11, palette = "Purple-Blu")


col_fun = function(x) ifelse(x < -0.07 | x > 0.07, colorspace::sequential_hcl(41, palette = "SunsetDark"), #
                             "#ECECEC"
                             # "#00000000"
)


# Base plot
chordDiagram(
  x = plot3, 
  #big.gap = 30,
  order = c(
    
    "MMP12"    ,   "PIGR"     ,   "CRP"         ,"IGFBP4"    ,  "SERPIND1"  ,  "SKR3"     , 
    "VEGFA"    ,   "FGF.21"   ,   "RARRES2"     ,"CCL18"     ,  "PRSS2"     ,  "THBS2"    , 
    "CCL17"    ,   "ICAM5"    ,   "MMP.1_olink" ,"MMP9"      ,  "MMP1"      ,  "HGF"      , 
    "B2M"      ,   "AFM"      ,   "SIGLEC1"     ,"TGF.alpha" ,  "SELE"      ,  "CXCL9"    , 
    "CXCL11"   ,   "CCL11"    ,   "ENPP7"       ,"NEP"       ,  "CCL25"     ,  "ACY1"     , 
    "MST1"     ,   "EN.RAGE"  ,   "G.CSF"       ,"S100A9"    ,  "RETN"      ,  "CLEC11A.1", 
    "FCGR3B"   ,   "OSM"      ,   "TNFRSF17"    ,"IDUA"      ,  "C4A|C4B"   ,  "CCL22"    , 
    "HGFAC"    ,   "CD5L",
    
    
    "g",
    "gf",
    "processing speed",
    "relative brain age",
    "total brain volume",
    "global grey matter",
    "global white matter",
    "global cortical volume",
    "global subcortical volume",
    "WMH",
    
    
    "gMD" ,
    "gFA"
    
    
  ),
  # col = col_fun,
  grid.col = plasma_pal,
  transparency = 0.25,
  symmetric = TRUE,
  #directional = 1,
  #direction.type = c("arrows", "diffHeight"), 
  #diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE 
  #link.largest.ontop = TRUE
)

# Add text and axis
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 3.5, 
      labels = sector.index, 
      #adj = c(0, degree(5)), 
      niceFacing = TRUE,
      facing = "clockwise",
      cex = 0.2
    )
    
    circos.rect(xleft=xlim[1], 
                ybottom=0.3, 
                xright=xlim[2], 
                ytop=0.32, 
                col = "white", 
                border = "white")
    
    
  }
)
circos.clear()

#######################
# good brain health plot
#######################

plot2 <- rbind(plot_neuroimaging_methylation, plot_cognitive_methylation)

##### For circos, certain associations that associate with POOR brain health neeed to be recoded
plot2 <- plot2 %>% mutate(estimate =
                            case_when(brain_metric == "relative brain age" ~ estimate*-1,
                                      brain_metric == "gMD" ~ estimate*-1,
                                      brain_metric == "WMH" ~ estimate*-1,
                                      TRUE ~ estimate)
)

plot2 %<>% filter(significance == "Yes" & estimate > 0)

####
plot3 <- plot2 %>% select(DNAm, brain_metric, estimate)


plot3 %<>% filter(brain_metric == "g" |
                    brain_metric == "gf" |
                    brain_metric == "processing speed" |
                    brain_metric == "relative brain age" |
                    brain_metric == "global subcortical volume" |
                    brain_metric == "global cortical volume" |
                    brain_metric == "WMH" |
                    brain_metric == "gFA" |
                    brain_metric == "gMD" |
                    brain_metric == "global grey matter" |
                    brain_metric == "global white matter" |
                    brain_metric == "total brain volume" )


plot3 %<>% 
  dplyr::rename(
    rowname = DNAm,
    key = brain_metric,
    value = estimate)

plot3$rowname <- as.character(plot3$rowname)

### arrange circos by size 
### Which DNAm proxy has the most significant hits? 
# Group by sum using dplyr
test <- aggregate((plot3$value), 
                  by=list(rowname=plot3$rowname), 
                  FUN=sum) %>% 
  #arrange((x))
arrange(desc(x))

test$rowname

test$rowname <- as.factor(test$rowname)
skimr::skim(test$rowname)
#n=33

#install.packages("circlize")
#library(circlize)
# parameters


circos.clear()


circos.par(gap.after = c("g" = 1.5,
                         "gf" = 1,
                         "processing speed" = 1,
                         "relative brain age"= 1,
                         "WMH" = 1,
                         "total brain volume"= 1,
                         "global grey matter"= 1,
                         "global white matter"= 1,
                         "global cortical volume"= 1,
                         "global subcortical volume"= 1,
                         "gMD"= 1,
                         "gFA"= 1.5,
                         
                          "SEMA3E"   =0.1,   "NCAM1"  =0.1  ,   "NOTCH1"  =0.1 ,   "NTRK3"      =0.1, "NTRK3_olink" =0.1,"CNTN4"   =0.1,   
                          "SLITRK5"  =0.1,   "SELL"   =0.1  ,   "TPSB2"   =0.1 ,   "OMD"        =0.1, "GDF.8"       =0.1,"WFIKKN2" =0.1,   
                          "ADAMTS13" =0.1,   "ESM1"   =0.1  ,   "GP1BA"   =0.1 ,   "GZMA_olink" =0.1, "SPOCK2"      =0.1,"SHBG"    =0.1,   
                          "GNLY"     =0.1,   "MRC2"   =0.1  ,   "EZR"     =0.1 ,   "INSR"       =0.1, "VCAM1"       =0.1,"FAP"     =0.1,   
                          "FcRL2"    =0.1,   "CD163"  =0.1  ,   "TNFRSF1B"=0.1 ,   "MMP2"       =0.1, "CD209"       =0.1,"EDA"     =0.1,   
                          "MPL"      =0.1,   "BCAM"   =0.1  ,   "CRTAM" =0.1
                         
                         
),

start.degree = 90, 
gap.degree = 0.1, 
track.margin = c(-0.1, 0.1), 
points.overflow.warning = FALSE
)
#par(mar = rep(0, 4))
#par(mfrow = c(1))
#


#install.packages("colorspace")



plasma_pal <- c(
  colorspace::sequential_hcl(33, palette = "BluYl"),
  #viridis::mako(n = 11)
  colorspace::sequential_hcl(11, palette = "Purple-Blu")
)


#pal2 <- colorspace::sequential_hcl(11, palette = "Purple-Blu")


col_fun = function(x) ifelse(x < -0.07 | x > 0.07, colorspace::sequential_hcl(41, palette = "SunsetDark"), #
                             "#ECECEC"
                             # "#00000000"
)


# Base plot
chordDiagram(
  x = plot3, 
  #big.gap = 30,
  order = c(
    
    "SEMA3E"   ,   "NCAM1"    ,   "NOTCH1"   ,   "NTRK3"      , "NTRK3_olink" ,"CNTN4"   ,   
    "SLITRK5"  ,   "SELL"     ,   "TPSB2"    ,   "OMD"        , "GDF.8"       ,"WFIKKN2" ,   
    "ADAMTS13" ,   "ESM1"     ,   "GP1BA"    ,   "GZMA_olink" , "SPOCK2"      ,"SHBG"    ,   
    "GNLY"     ,   "MRC2"     ,   "EZR"      ,   "INSR"       , "VCAM1"       ,"FAP"     ,   
    "FcRL2"    ,   "CD163"    ,   "TNFRSF1B" ,   "MMP2"       , "CD209"       ,"EDA"     ,   
    "MPL"      ,   "BCAM"     ,   "CRTAM" ,
    
    
    "g",
    "gf",
    "processing speed",
    "relative brain age",
    "total brain volume",
    "global grey matter",
    "global white matter",
    "global cortical volume",
    "global subcortical volume",
    "WMH",
    
    
    "gMD" ,
    "gFA"
    
    
  ),
  # col = col_fun,
  grid.col = plasma_pal,
  transparency = 0.25,
  symmetric = TRUE,
  #directional = 1,
  #direction.type = c("arrows", "diffHeight"), 
  #diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE 
  #link.largest.ontop = TRUE
)

# Add text and axis
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
     # Add names to the sector. 
     circos.text(
       x = mean(xlim), 
       y = 3.5, 
       labels = sector.index, 
       #adj = c(0, degree(5)), 
       niceFacing = TRUE,
       facing = "clockwise",
       cex = 0.2
     )
    
    circos.rect(xleft=xlim[1], 
                ybottom=0.3, 
                xright=xlim[2], 
                ytop=0.32, 
                col = "white", 
                border = "white")
    
    
  }
)
circos.clear()



############


plot2 <- rbind(plot_neuroimaging_methylation, plot_cognitive_methylation)

##### For circos, certain associations that associate with POOR brain health neeed to be recoded
plot2 <- plot2 %>% mutate(estimate =
                            case_when(brain_metric == "relative brain age" ~ estimate*-1,
                                      brain_metric == "gMD" ~ estimate*-1,
                                      brain_metric == "WMH" ~ estimate*-1,
                                      TRUE ~ estimate)
)

plot2 %<>% filter(significance == "Yes" & estimate < 0)

####

plot2 %<>% filter(brain_metric == "g" |
                    brain_metric == "gf" |
                    brain_metric == "processing speed" |
                    brain_metric == "relative brain age" |
                    brain_metric == "global subcortical volume" |
                    brain_metric == "global cortical volume" |
                    brain_metric == "WMH" |
                    brain_metric == "gFA" |
                    brain_metric == "gMD" |
                    brain_metric == "global grey matter" |
                    brain_metric == "global white matter" |
                    brain_metric == "total brain volume" )




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


