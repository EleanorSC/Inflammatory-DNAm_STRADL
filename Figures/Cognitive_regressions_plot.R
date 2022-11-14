# examine data

STRADL_cogn_full <- readxl::read_excel("STRADL_Laura_variables.xlsx",sheet=4)

STRADL_cogn <- STRADL_cogn_full %>%
  
  dplyr::select(ID,
                mema, #Logical Memory I: Story A immediate recall
                memdela, #Logical Memory I: Story A delayed recall
                digsym, #WAIS III - Digit Symbol-coding total score
                vftot, #Verbal Fluency total score
                mhv, #Mill Hill Vocabulary total score
                mrtotc #Matrix Reasoning total correct
                # mrtottime #Matrix Reasoning total time taken (seconds)
  ) %>%
  
  rename(stradl_ID = ID, 
         mema_STRADL = mema, 
         memdela_STRADL = memdela, 
         Digit_Symbol_coding = digsym, 
         Verbal_Fluency = vftot, 
         Mill_Hill_vocab = mhv,
         Matrix_Reasoning = mrtotc
         #  matrix_reasoning = mrtottime
  ) %>%
  mutate(Logical_Memory = mema_STRADL+memdela_STRADL)

#Examine data
skimr::skim(STRADL_cogn)

# prepare data
cognition_STRADL <- STRADL_cogn %>%
  dplyr::select(stradl_ID, 
                Digit_Symbol_coding, # Processing speed
                Verbal_Fluency, # Executive function
                Mill_Hill_vocab,
                Matrix_Reasoning,
                Logical_Memory ## n.b also verbal declarative memory
  )

#Examine data
skimr::skim(cognition_STRADL)

#skimr::skim(PROTEOMICS_DNAm_DATA)


PROTEOMICS_DNAm_DATA_COGNITION <- merge(PROTEOMICS_DNAm_DATA, 
                        cognition_STRADL,
                        by = "stradl_ID")

#PROTEOMICS_DNAm_DATA_COGNITION <- merge(PROTEOMICS_DNAm_DATA_COGNITION, 
#                        pca_component,
#                        by = "stradl_ID")

## ----------------------------#

skimr::skim(PROTEOMICS_DNAm_DATA_COGNITION)

#glimpse(PROTEOMICS_DNAm_DATA_COGNITION$pca_factor)

PROTEOMICS_DNAm_DATA_COGNITION %<>% rename( 
 # g_STRADL = pca_factor, 
  processing_speed = Digit_Symbol_coding,
  executive_function = Verbal_Fluency, 
  verbal_declarative_memory = Logical_Memory
) 


# Converting multiple varibles into a factor
PROTEOMICS_DNAm_DATA_COGNITION %<>% mutate_at(c("sex", "Site"),
                              as.factor)


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

### ----------------------------# 
## EXAMINING PROTEIN ASSOCIATIONS WITH COG FUNCTIONS
### ----------------------------#
#select_proteins <- protein_data_only %>% 
#  select(
#    stradl_ID,
#    ACY1    ,
#    ADAMTS13  ,
#    ADIPOQ,
#    AFM       ,
#    B2M      ,
#    BCAM     ,
#    BMP1   ,
#    `C4A|C4B` ,
#    C5      ,
#    C9        ,
#    CCL17     ,
#    CCL18     ,
#    CCL21    ,
#    CCL22    ,
#    CCL25  ,
#    CD163   ,
#    CD209   ,
#    CD48      ,
#    CD5L      ,
#    CHIT1     ,
#    CLEC11A  ,
#    CLEC11A.1,
#    CNTN4  ,
#    CRP     ,
#    CXCL10  ,
#    CXCL11    ,
#    EDA       ,
#    ENPP7     ,
#    ESM1     ,
#    F7       ,
#    FAP    ,
#    FCER2   ,
#    FCGR3B  ,
#    GHR       ,
#    GNLY      ,
#    GP1BA     ,
#    GZMA     ,
#    HGFAC    ,
#    ICAM5  ,
#    IDUA    ,
#    IGFBP1  ,
#    IGFBP4    ,
#    IL19      ,
#    INSR      ,
#    LGALS3BP ,
#    LGALS4   ,
#    `LTA|LTB`,
#    LTF     ,
#    LY9     ,
#    LYZ       ,
#    MIA       ,
#    MMP1     ,
#    MMP12    ,
#    MMP2   ,
#    MMP9    ,
#    MPL     ,
#    MPO       ,
#    MRC2      ,
#    MST1      ,
#    NCAM1    ,
#    NOTCH1   ,
#    NTRK3  ,
#    OMD     ,
#    PAPPA   ,
#    PIGR      ,
#    PRSS2     ,
#    RARRES2   ,
#    RETN     ,
#    S100A9   ,
#    SELE   ,
#    SELL    ,
#    SEMA3E  ,
#    SERPINA3  ,
#    SERPIND1  ,
#    SHBG      ,
#    SLITRK5  ,
#    SPOCK2   ,
#    STC1   ,
#    THBS2   ,
#    TNFRSF17  ,
#    TNFRSF1B  ,
#    TPSB2     ,
#    VCAM1    ,
#    WFIKKN2
#  )
#
#
#select_proteins <- merge(STRADL_lifestyle_covariates,
#                         select_proteins,
#                         by = "stradl_ID")
#
#PROTEINS_COGNITION <- merge(select_proteins , 
#                            cognition_STRADL,
#                            by = "stradl_ID")
#
#PROTEINS_COGNITION <- merge(PROTEINS_COGNITION, 
#                            pca_component,
#                            by = "stradl_ID")
#
### ----------------------------#
#
##skimr::skim(PROTEOMICS_DNAm_DATA_COGNITION)
#
##glimpse(PROTEOMICS_DNAm_DATA_COGNITION$pca_factor)
#
#PROTEINS_COGNITION %<>% rename( 
#  g_STRADL = pca_factor, 
#  processing_speed = Digit_Symbol_coding,
#  executive_function = Verbal_Fluency, 
#  verbal_declarative_memory = Logical_memory
#) 
#
#
## Converting multiple varibles into a factor
#PROTEINS_COGNITION %<>% mutate_at(c("sex", "Site"),
#                                  as.factor)
#
#
#
## Making a data frame with all combinations of variables
#df <-
#  as.data.frame(expand.grid(FULL_DNAm_list, FULL_cognitive_list, stringsAsFactors = F)) %>%
#  dplyr::rename(DNAm = Var1, metric = Var2)
#
## Function to get summary values as tibble from a named list with the info on metric and DNAm
#get_summary_values <- function(list) {
#  metric <- list$metric
#  DNAm <- list$DNAm
#  tib <-
#    broom::tidy(summary(
#      lm(
#        scale(PROTEINS_COGNITION[[metric]]) ~
#          scale(st_age)
#        + sex
#        + scale(PROTEINS_COGNITION[[DNAm]]),
#        data = PROTEINS_COGNITION
#      )
#    ))[4, c(2, 3, 5)]
#  return(tib)
#}
#
## Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
#newcols <- df %>%
#  split(1:nrow(.)) %>%
#  purrr::map_dfr(.f = get_summary_values)
#
## Now bind on the 3 cols that were created
#newdf <- cbind(df, newcols) %>%
#  
#  # Create a neurimaging modality column
#  mutate(modality = "cognition") %>%
#  
#  # Create a brain metric column
#  mutate(
#    brain_metric =
#      case_when(
#        
#        metric == "g" ~ "g",
#        metric == "gf" ~ "gf",
#        metric == "g_STRADL" ~ "g_STRADL",
#        metric == "processing_speed" ~ "processing speed",
#        metric == "executive_function" ~ "executive function",
#        metric == "verbal_declarative_memory" ~ "verbal declarative memory",
#        
#        TRUE ~ "misc"
#      )
#  ) %>%
#  
#  # Create pFDR column
#  # group_by(metric, Hemisphere) %>%
#  group_by(DNAm) %>%
#  mutate(pFDR = p.adjust(p.value, method = "fdr")) %>%
#  # Create a column to denote where there are significant hits
#  mutate(significance =
#           case_when(p.value < 0.05 ~ "Yes",
#                     TRUE ~ "No")) %>%
#  # Create a column to denote where there are FDR significant hits
#  mutate(FDR_significance =
#           case_when(pFDR < 0.05 ~ "Yes",
#                     TRUE ~ "No"))%>%
#  
#  # Create a column to denote whether these are DNAm or proteins
#  mutate(omic_type = "protein")
#
#
#plot_cognitive_proteins <- newdf
#
### ----------------------------# 
## COMBINE ANALYSES TABLES
### ----------------------------#
#
#plot <- rbind(plot_cognitive_methylation, plot_cognitive_proteins)

## ----------------------------# 
# PLOT ANALYSES
## ----------------------------#

plot <- plot_cognitive_methylation 

plot2 <- plot


##### REORDER cognitive by global then domains
plot2 %<>% mutate(brain_metric = factor(
  brain_metric,
  levels = c(#"g_STRADL",
             "gf", "g",
             "processing speed", 
             "executive function",
             "verbal declarative memory")))

##### for effect sizes plot
facetSettings <-
  theme(strip.background = element_rect(
    #fill = "#F8EEEC",
    fill = "#F3FFF9",
    colour = "black",
    size = 1
  ))

ggplot(plot2,
       
       aes(x = reorder(DNAm,-(estimate)),
           y = estimate,
           alpha = reorder(DNAm,
                         (-estimate)),
           shape = significance,
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
                               option = "G") +
  scale_shape_manual(values = c(1,
                                16)) +
  facet_wrap( ~ brain_metric,
              #scales="free",
              #scales = "free_y",
              nrow = 1) +
  facetSettings 
#+
#  scale_alpha_continuous(range = c(0.3, 1)) +
  
# 11 x 25

## ----------------------------# 
# PLOT ANALYSES - ONLY SIGNIFICANT ASSOCIATIONS
## ----------------------------#
plot <- plot_cognitive_methylation

### Which DNAm proxy has the most significant hits? I should be able to do this...
plot %<>%
  mutate(
    number_significant =
      case_when(significance == "No" ~ 0,
                TRUE ~ 1))

# Group by sum using dplyr
aggregate(plot$number_significant, 
          by=list(brain_metric=plot$brain_metric), 
          FUN=sum)


plot2 <- plot %>% filter(significance == "Yes")

##### REORDER cognitive by global then domains
plot2 %<>% mutate(brain_metric = factor(
  brain_metric,
  levels = c(#"g_STRADL",
    "gf", "g",
    "processing speed", 
    "executive function",
    "verbal declarative memory")))

##### for effect sizes plot
facetSettings <-
  theme(strip.background = element_rect(
    #fill = "#F8EEEC",
    #fill = "#F3FFF9",
    fill = "#f3f1ff",
    colour = "black",
    size = 1
  ))

ggplot(plot2,
       
       aes(x = reorder(DNAm,-(estimate)),
           y = estimate,
           alpha = reorder(DNAm,
                           (-estimate)),
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
                               option = "G") +
  scale_shape_manual(values = c(16, 
                                5
                                #16
                                )) +
  facet_wrap( ~ brain_metric,
              scales="free_y",
              #scales = "free_y",
              nrow = 1) +
  facetSettings 


### ----------------------------# 
## PLOT individually then cocatenate 
### ----------------------------#


### ----------------------------# 
## PLOT ANALYSES - ONLY SIGNIFICANT ASSOCIATIONS
### ----------------------------#
#plot <- rbind(plot_cognitive_methylation, plot_cognitive_proteins)
#
#plot2 <- plot %>% filter(significance == "Yes")
#
###### REORDER cognitive by global then domains
#plot2 %<>% mutate(brain_metric = factor(
#  brain_metric,
#  levels = c("g_STRADL", "g", "gf",
#             "processing speed", 
#             "executive function",
#             "verbal declarative memory")))
#
#
#ggplot(plot2,
#       
#       aes(x = reorder(DNAm,-(estimate)),
#           y = estimate,
#           col = reorder(DNAm,
#                         (estimate)),
#           shape = omic_type,
#           group = omic_type,
#           alpha = FDR_significance
#       )
#) +
#  
#  geom_point(position = position_dodge(width = 0.9),
#             size = 2.2,
#             #colour = "#414487FF",
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
#    axis.text.x = element_text(#angle = 90,
#      #vjust = 0.5,
#      # hjust=1,
#      size = 6),
#    #axis.text.x = element_blank(),
#    strip.text = element_text(
#      size = 6,
#      face = "bold",
#      family = "sans",
#      colour = "black"
#    ),
#    axis.text.y = element_text(size = 7)
#  ) +
#  
#  labs(x = "",
#       y = "") +
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
#                               option = "A") +
#  scale_shape_manual(values = c(16,
#                                1)) +
#  facet_wrap( omic_type ~ brain_metric,
#              #scales="free"
#              scales = "free_y",
#              nrow = 1) +
#  facetSettings