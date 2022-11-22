## ---------------------------
##
## Script Purpose: FA / MD  associations with DNAm metrics 
##                
##                (nb.) Neuroimaging models must contain key covariates of (site / edits / batch / estimated ICV) alongside age and sex
## 
##                 This script does the following:
##
##                (1) Creates supplementary table for pFDR significant cortical volume associations;
##                (2) Creates barplot for pFDR significant cortical volume associations facetted by whether they associate with lower or increased volumes
##
## ---------------------------


# neuroimaging data including DTI data
STRADL_WM_tracts_FA <- read.csv("STRADL_Measures_DTI_FA.csv")
names(STRADL_WM_tracts_FA)[names(STRADL_WM_tracts_FA) == 'subjectID'] <- 'stradl_ID'
# drop column 'average FA'
STRADL_WM_tracts_FA <- STRADL_WM_tracts_FA[,-8]

# neuroimaging data including DTI data
STRADL_WM_tracts_MD <- read.csv("STRADL_Measures_DTI_MD.csv")
names(STRADL_WM_tracts_MD)[names(STRADL_WM_tracts_MD) == 'subjectID'] <- 'stradl_ID'
# drop column 'average MD'
STRADL_WM_tracts_MD <- STRADL_WM_tracts_MD[,-8]

skimr::skim(STRADL_WM_tracts_FA)

#colnames(STRADL_WM_tracts_FA) <- paste(colnames(STRADL_WM_tracts_FA),
#                                   "FA",
#                                   sep="_")

# n =709
DTI_DNAm_FA <- merge(PROTEOMICS_DNAm_DATA,
                  STRADL_WM_tracts_FA,
                           by = "stradl_ID")

DTI_DNAm_MD <- merge(PROTEOMICS_DNAm_DATA,
                     STRADL_WM_tracts_MD,
                     by = "stradl_ID")


# Create list to loop through all DNAm signatures
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
  "SPOCK2"
)


#note change FXST -> FX.ST
names(DTI_DNAm_FA)[names(DTI_DNAm_FA) == "FXST"] <- "FX.ST"
names(DTI_DNAm_MD)[names(DTI_DNAm_MD) == "FXST"] <- "FX.ST"
######
FULL_neuroimaging_list <- c("ACR",
                               "ALIC",
                               "BCC",
                               "CC",
                               "CGC",
                               "CGH",
                               "CR",
                               "CST",
                               "EC",
                               "FX",
                               "FX.ST",
                               "GCC",
                               "IC",
                               "IFO",
                               "PCR",
                               "PLIC",
                               "PTR",
                               "RLIC",
                               "SCC",
                               "SCR",
                               "SFO",
                               "SLF",
                               "SS",
                               "UNC")

neuroimaging_regressions <- list()


# Converting multiple varibles into a factor
DTI_DNAm_FA %<>% mutate_at(c("sex", "Site"),
                                 as.factor)

# Making a data frame with all combinations of variables
df <-
  as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1,
                metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
model_output <- function(list) {
  
  metric <- list$metric
  DNAm <- list$DNAm
  
  regression_summary <- summary(
    lm(
      scale(DTI_DNAm_FA[[metric]]) ~
        scale(st_age)
      + sex
      + Site
      + scale(DTI_DNAm_FA[[DNAm]]),
      data = DTI_DNAm_FA
    )
  )
  
  regression_summary_tidy <- broom::tidy(regression_summary)
  regression_summary_tidy_complete <- regression_summary_tidy %>% mutate(r2 = regression_summary$r.squared)
  regression_summary_tidy_DNAm <- regression_summary_tidy_complete[5, c(2, 3, 5, 6)]
  
  return(regression_summary_tidy_DNAm)
  
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- (df) %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = model_output)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality = "FA") %>%
  
  # Create a brain metric column
  mutate(
    brain_metric = metric) %>%
  
  # Create pFDR column
  group_by(brain_metric, modality) %>%
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
  mutate(omic_type = "DNAm") %>%
  mutate(model = "Model 1 (n=683)")

###
plot_FA_methylation <- newdf




#################### MD

# Converting multiple varibles into a factor
DTI_DNAm_MD %<>% mutate_at(c("sex", "Site"),
                           as.factor)

# Making a data frame with all combinations of variables
df <-
  as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1,
                metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
model_output <- function(list) {
  
  metric <- list$metric
  DNAm <- list$DNAm
  
  regression_summary <- summary(
    lm(
      scale(DTI_DNAm_MD[[metric]]) ~
        scale(st_age)
      + sex
      + Site
      + scale(DTI_DNAm_MD[[DNAm]]),
      data = DTI_DNAm_MD
    )
  )
  
  regression_summary_tidy <- broom::tidy(regression_summary)
  regression_summary_tidy_complete <- regression_summary_tidy %>% mutate(r2 = regression_summary$r.squared)
  regression_summary_tidy_DNAm <- regression_summary_tidy_complete[5, c(2, 3, 5, 6)]
  
  return(regression_summary_tidy_DNAm)
  
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- (df) %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = model_output)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality = "MD") %>%
  
  # Create a brain metric column
  mutate(
    brain_metric = metric) %>%
  
  # Create pFDR column
  group_by(brain_metric, modality) %>%
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
  mutate(omic_type = "DNAm") %>%
  mutate(model = "Model 1 (n=683)")

###
plot_MD_methylation <- newdf



# ----------------------------#
# Combine cognitive and brain analyses
# ----------------------------#

plot_DTI_m1 <- rbind(plot_FA_methylation,
                     plot_MD_methylation)

# ----------------------------#
# Examine which DNAm associate across multiple WM tracts
# ----------------------------#

### Which DNAm proxy has the most FDR significant hits?
sigs <- plot_DTI_m1  %>%
  #filter(estimate < 0 & FDR_significance == "Yes") %>%
  filter(significance == "Yes" & modality == "FA" & estimate < 0) %>%
  arrange(estimate) %>%
  group_by(brain_metric)%>%
  group_by(DNAm, brain_metric) %>%
  mutate(number_significant =
           case_when(FDR_significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_sigs <-
  aggregate(sigs$number_significant,
            by = list(DNAm = sigs$DNAm),
            FUN = sum) %>%
  arrange(desc(x)) %>%
  rename(n_pFDR = x)



# Reorder to find which tract has the most sig hits
Top_hits_sigs <-
  aggregate(sigs$number_significant,
            by = list(brain_metric = sigs$brain_metric),
            FUN = sum) %>%
  arrange(desc(x)) %>%
  rename(n_pFDR = x)

x <-Top_hits_sigs$brain_metric
# ----------------------------#
# Examine these individual DNAms
# ----------------------------#

SERPIND_sigs <- sigs %>% filter(DNAm == "SERPIND1" & FDR_significance == "Yes")

# ----------------------------#
# Table of significant DTI regressions for supplementary document
# ----------------------------#
# Reorder to find which tract has the most sig hits

sigs <- plot_DTI_m1  %>%
  #filter(estimate < 0 & FDR_significance == "Yes") %>%
  filter(significance == "Yes" & modality == "FA" & estimate < 0) %>%
  arrange(estimate) %>%
  group_by(brain_metric)%>%
  group_by(DNAm, brain_metric) %>%
  mutate(number_significant_p =
           case_when(significance == "No" ~ 0,
                     TRUE ~ 1))
Top_hits_sigs <-
  aggregate(sigs$number_significant_p,
            by = list(brain_metric = sigs$brain_metric),
            FUN = sum) %>%
  arrange(desc(x)) %>%
  rename(n_p = x)

x <-Top_hits_sigs$brain_metric


# ----------------------------#
# PHEWAS PLOT
# ----------------------------#


plot_DTI_m1_sig  <- plot_DTI_m1  %>%
  filter(FDR_significance == "Yes" & modality == "FA") %>%
  arrange(estimate) %>%
  group_by(brain_metric)


plot2 <- plot_DTI_m1  %>%
  filter(modality == "FA" & estimate < 0) %>%
  mutate(
    SIG = case_when(
      significance == "Yes" & FDR_significance == "Yes" ~ "2",
      significance == "Yes" &
        FDR_significance == "No" ~ "1",
      TRUE ~ "0"
    )
  )


##### REORDER cognitive by global then domains
plot2 %<>% mutate(significance = factor(significance,
                                        levels = c("No", "Yes"))) %>% 
  mutate(DNAm = str_replace(DNAm, "_olink", "x")) %>% 
  
  #### order in terms of FDR significance 
  
mutate(brain_metric = factor(brain_metric,
                            levels = c( "FX"   , "SFO"   ,"SS"   , "CR"   , "ALIC"  ,"PTR"   ,"PCR" ,  "SCR" ,
                                        "EC"   , "GCC"  , "SLF"  , "CGC" ,  "ACR" ,  "BCC"   ,"CC"  ,  "CGH" , 
                                        "CST"  , "FX.ST" ,"IC"   , "IFO"  , "PLIC"  ,"RLIC"  ,"SCC" ,  "UNC"  ) 
                            ))%>% 


#mutate(brain_metric = factor(brain_metric,
#                             levels = c(
#                               "SS"  ,  "FX"   , "CR"    ,"PCR"   ,"SFO"  , "PTR"  , "CC"   , "GCC"   ,"FX.ST",
#                               "SCR" ,  "ACR"  , "BCC" ,  "EC"    ,"SLF" ,  "ALIC" , "CGC"  ,
#                                "IC" ,   "SCC" ,  "RLIC" , "UNC"  , "IFO" ,  "CGH" ,  "CST" ,  "PLIC" ))) %>% 
#  
  # Only plot tracts where there is at least one instance of > 0.05 pFDR
  
  filter(brain_metric == "SS" |
           brain_metric == "FX"|
           brain_metric =="SFO"   |
           brain_metric == "CR"  |
           brain_metric =="ALIC" | 
           brain_metric =="PTR" |  
           brain_metric =="PCR"  | 
           brain_metric =="SCR"  | 
           brain_metric =="EC"   |
           brain_metric =="GCC"  | 
           brain_metric =="SLF"  |
           brain_metric =="CGC" )

ggplot(plot2,
       aes(
         x = reorder(DNAm, estimate),
         y = -log(p.value)
       )) +
  
  geom_point(aes(
    colour = brain_metric,
    alpha = significance,
    shape = SIG),
    size = 1.3
    ) +
    
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
           x = "",
           y = "") +
  
  ggrepel::geom_text_repel(
    data = . %>% mutate(label = ifelse(pFDR < 0.05, as.character(DNAm),
                                        "")),
    aes(label = label),
    size = 1.5,
    colour = "black",
    min.segment.length = Inf,
    max.overlaps = Inf) +
  
  ggrepel::geom_text_repel(
    data = . %>% mutate(label2 = ifelse(p.value < 0.05 & pFDR > 0.05, 
                                       as.character(DNAm),
                                       "")),
    aes(label = label2),
    size = 1.5,
    colour = "darkgrey",
    alpha = 0.8,
    min.segment.length = Inf,
    max.overlaps = Inf) +

      geom_hline(
        yintercept = -log(0.05),
        color = "darkgrey",
        size = 1,
        alpha = 0.5
      ) +
      
      scale_shape_manual(values = c(16,
                                    16,
                                    8)) +
      
      facet_wrap(~ brain_metric
                 #scales = "free_x"
                 ) +
      
      scale_colour_manual(values = c(
        colorspace::qualitative_hcl(24, palette = "Cold") 
        ))

# ----------------------------#
# END PHEWAS PLOT
# ----------------------------#

# ----------------------------#
# EFFECT SIZE RANGE PLOT
# ----------------------------#

plot3 <- plot2 %>% filter(significance == "Yes")

facetSettings <-
  theme(strip.background = element_rect(
    fill = "#EAE3F2",
    colour = "black",
    size = 1
  ))

ggplot(plot3,
       
       aes(
         x = reorder(DNAm, -(estimate)),
         y = estimate,
         alpha = significance
        # col = brain_metric
         
       )) +
  
  geom_point(
    aes(
      col = brain_metric,
      alpha = reorder(DNAm, (-estimate)),
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
  
  scale_colour_manual(values = c(
    colorspace::qualitative_hcl(12, palette = "Cold") 
  )) +
  
  facet_wrap( ~ brain_metric,
              scales = "free_y"

  ) +
  facetSettings +
 # scale_alpha_manual(values=c(0.8, 1, 0.2)) +
  scale_shape_manual(values = c(
    # 8,
    # 1,
    #  1,
    16,
    5))

# ----------------------------#
# END EFFECT SIZE RANGE PLOT
# ----------------------------#


Poor_brain_health <- plot_DTI_m1 %>% filter(significance == "Yes" & estimate < 0 & modality == "FA")

### Which DNAm proxy has the most FDR significant hits?
Poor_brain_health %<>%
  mutate(number_significant =
           case_when(significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Poor_brain_health <-
  aggregate(
    Poor_brain_health$number_significant,
    by = list(brain_metric = Poor_brain_health$brain_metric),
    FUN = sum
  )

Top_hits_Poor_brain_health %<>% arrange(desc(x)) %>% rename(n_p = x)

# Run FDR
Poor_brain_health %<>%
  mutate(number_significant_FDR =
           case_when(FDR_significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Poor_brain_health_2 <-
  aggregate(
    Poor_brain_health$number_significant_FDR,
    by = list(brain_metric = Poor_brain_health$brain_metric),
    FUN = sum
  ) %>% rename(n_pFDR = x)


Poor_brain <- merge(Top_hits_Poor_brain_health,
                    Top_hits_Poor_brain_health_2,
                    by = "brain_metric") %>%
  mutate(direction =
           "Lower white matter integrity (FA)")

Significance_plot <- Poor_brain

# ----------------------------#
# Code for Barplot of count of significant hits
# ----------------------------#

####
ggplot(Significance_plot,
  
  aes(
    y = reorder(brain_metric, n_p),
    x = n_p
    
  )) +

  
  geom_col(aes(y = reorder(brain_metric, n_p),
               x = n_p,
               fill = reorder(brain_metric, n_p))) +
  
  geom_col(aes(
    y = reorder(brain_metric, n_p),
    x = n_pFDR,
    colour = c("grey"),
    alpha = 0.95
  )) +
  
  
  theme_classic() +
  
  theme(
    legend.position = "none",
    
    strip.text = element_text(
      size = 6,
      face = "bold",
      family = "sans",
      colour = "black"
    ),
    axis.text.x = element_text(
      vjust = 0.5,
      hjust = 1,
      size = 6,
      family = "sans"
    ),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(
      size = 7,
      face = "bold",
      family = "sans"
    ),
    axis.title.x = element_text(
      size = 7,
      face = "bold",
      family = "sans"
    )
  ) +
  
  labs(x = "number of significant assocations with DNAm signatures",
       y = "global brain metrics") +
  
  facet_wrap( ~ direction) +
  scale_fill_manual(values = c(
    # viridis::magma(n = 19)
    #colorspace::qualitative_hcl(24, palette = "Cold") 
    colorspace::sequential_hcl(24, palette = "Purple-Blu") 
    #colorspace::sequential_hcl(24, palette = "SunsetDark")
  )) +
  scale_colour_manual(values = c("#808080"))


####
## ----------------------------# 
# MODEL 3 - LIFESTYLE
## ----------------------------#

DTI_DNAm_FA_lifestyle <- merge(Lifestyle_DNAm, 
                               STRADL_WM_tracts_FA,
                                     by = "stradl_ID")

names(DTI_DNAm_FA_lifestyle)[names(DTI_DNAm_FA_lifestyle) == "FXST"] <- "FX.ST"

# Making a data frame with all combinations of variables
df <-
  as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1,
                metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
model_output <- function(list) {
  
  metric <- list$metric
  DNAm <- list$DNAm
  
  regression_summary <- summary(
    lm(
      scale(DTI_DNAm_FA_lifestyle[[metric]]) ~
        scale(st_age)
      + sex
      + Site
      
      + hypertension
      + CurrentSmoker
      + CurrentDrinker
      + scale(bmi)
      
      + scale(DTI_DNAm_FA_lifestyle[[DNAm]]),
      data = DTI_DNAm_FA_lifestyle
    )
  )
  
  regression_summary_tidy <- broom::tidy(regression_summary)
  regression_summary_tidy_complete <- regression_summary_tidy %>% mutate(r2 = regression_summary$r.squared)
  regression_summary_tidy_DNAm <- regression_summary_tidy_complete[9, c(2, 3, 5, 6)]
  
  return(regression_summary_tidy_DNAm)
  
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- (df) %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = model_output)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality = "FA") %>%
  
  # Create a brain metric column
  mutate(
    brain_metric = metric) %>%
  
  # Create pFDR column
  group_by(brain_metric, modality) %>%
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
  mutate(omic_type = "DNAm") %>%
  mutate(model = "Model 2 (n=683)")

###
plot_FA_methylation_lifestyle <- newdf


Poor_brain_health <- plot_FA_methylation_lifestyle %>% filter(significance == "Yes" & estimate < 0 & modality == "FA")

### Which DNAm proxy has the most FDR significant hits?
Poor_brain_health %<>%
  mutate(number_significant =
           case_when(significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Poor_brain_health <-
  aggregate(
    Poor_brain_health$number_significant,
    by = list(brain_metric = Poor_brain_health$brain_metric),
    FUN = sum
  )

Top_hits_Poor_brain_health %<>% arrange(desc(x)) %>% rename(n_p = x)

# Run FDR
Poor_brain_health %<>%
  mutate(number_significant_FDR =
           case_when(FDR_significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Poor_brain_health_2 <-
  aggregate(
    Poor_brain_health$number_significant_FDR,
    by = list(brain_metric = Poor_brain_health$brain_metric),
    FUN = sum
  ) %>% rename(n_pFDR = x)


Poor_brain <- merge(Top_hits_Poor_brain_health,
                    Top_hits_Poor_brain_health_2,
                    by = "brain_metric") %>%
  mutate(direction =
           "Lower white matter integrity (FA)")

Significance_plot <- Poor_brain

# ----------------------------#
# Code for Barplot of count of significant hits
# ----------------------------#
Significance_plot %<>%
  
  mutate(brain_metric = factor(brain_metric,
                               levels = c(     
                                          
                                           
                                           "CGH",
                                           "CST","PLIC","IFO","UNC",   "RLIC",  "CGC","IC", "SCC", "ALIC", "ACR", "BCC", "EC","SLF", 
                                           "FX.ST", "SCR","CC",  "GCC",  "PTR","PCR","SFO","CR","FX","SS" 
                                           ))) 

####
ggplot(Significance_plot,
       
       aes(
         y = brain_metric,
         x = n_p
         
       )) +
  
  
  geom_col(aes(
    y = brain_metric,
               x = n_p,
    fill = brain_metric
    )) +
  
  geom_col(aes(
    y = brain_metric,
    x = n_pFDR,
    colour = c("grey"),
    alpha = 0.95
  )) +
  
  
  theme_classic() +
  
  theme(
    legend.position = "none",
    
    strip.text = element_text(
      size = 6,
      face = "bold",
      family = "sans",
      colour = "black"
    ),
    axis.text.x = element_text(
      vjust = 0.5,
      hjust = 1,
      size = 6,
      family = "sans"
    ),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(
      size = 7,
      face = "bold",
      family = "sans"
    ),
    axis.title.x = element_text(
      size = 7,
      face = "bold",
      family = "sans"
    )
  ) +
  
  labs(x = "number of significant assocations with DNAm signatures",
       y = "global brain metrics") +
  
  facet_wrap( ~ direction) +
  scale_fill_manual(values = c(
    # viridis::magma(n = 19)
    #colorspace::qualitative_hcl(24, palette = "Cold") 
    colorspace::sequential_hcl(24, palette = "Purple-Blu") 
    #colorspace::sequential_hcl(24, palette = "SunsetDark")
  )) +
  scale_colour_manual(values = c("#808080")) +
  xlim(0, 22)


