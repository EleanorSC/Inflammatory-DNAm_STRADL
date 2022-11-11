FULL_DNAm_list <- c(
  "ACY1"    ,
  "ADAMTS13"  ,
  "ADIPOQ",
  "AFM"       ,
  "B2M"      ,
  "BCAM"     ,
  "BMP1"   ,
  "C4A|C4B" ,
  "C5"      ,
  "C9"        ,
  "CCL17"     ,
  "CCL18"     ,
  "CCL21"    ,
  "CCL22"    ,
  "CCL25"  ,
  "CD163"   ,
  "CD209"   ,
  "CD48"      ,
  "CD5L"      ,
  "CHIT1"     ,
  "CLEC11A"  ,
  "CLEC11A.1",
  "CNTN4"  ,
  "CRP"     ,
  "CXCL10"  ,
  "CXCL11"    ,
  "EDA"       ,
  "ENPP7"     ,
  "ESM1"     ,
  "F7"       ,
  "FAP"    ,
  "FCER2"   ,
  "FCGR3B"  ,
  "GHR"       ,
  "GNLY"      ,
  "GP1BA"     ,
  "GZMA"     ,
  "HGFAC"    ,
  "ICAM5"  ,
  "IDUA"    ,
  "IGFBP1"  ,
  "IGFBP4"    ,
  "IL19"      ,
  "INSR"      ,
  "LGALS3BP" ,
  "LGALS4"   ,
  "LTA|LTB",
  "LTF"     ,
  "LY9"     ,
  "LYZ"       ,
  "MIA"       ,
  "MMP1"     ,
  "MMP12"    ,
  "MMP2"   ,
  "MMP9"    ,
  "MPL"     ,
  "MPO"       ,
  "MRC2"      ,
  "MST1"      ,
  "NCAM1"    ,
  "NOTCH1"   ,
  "NTRK3"  ,
  "OMD"     ,
  "PAPPA"   ,
  "PIGR"      ,
  "PRSS2"     ,
  "RARRES2"   ,
  "RETN"     ,
  "S100A9"   ,
  "SELE"   ,
  "SELL"    ,
  "SEMA3E"  ,
  "SERPINA3"  ,
  "SERPIND1"  ,
  "SHBG"      ,
  "SLITRK5"  ,
  "SPOCK2"   ,
  "STC1"   ,
  "THBS2"   ,
  "TNFRSF17"  ,
  "TNFRSF1B"  ,
  "TPSB2"     ,
  "VCAM1"    ,
  "WFIKKN2"
)

FULL_neuroimaging_list <- c(
  ### Global measures
  "global.cerebral.wm",
  "global.total.gm",
  "global.wbv",
  # "global.wbv.vent",
  "global_cortical_volume",
  "global_cortical_thickness",
  "global_cortical_surface_area"
)

## Load up DNAm data
STRADL_KORA <- read.csv("STRADL_scores_trained_KORA.csv", check.names=FALSE)

newnames <- c()
for(colname in names(STRADL_KORA)){
  if(colname %in% annotation_data$SeqId){
    egname <- annotation_data[annotation_data$SeqId == colname,]$`Entrez Gene Name`
    newnames <- append(newnames, egname)
  } else {
    newnames <- append(newnames, colname)
  }
}

print(newnames)
names(STRADL_KORA) <- newnames

# Now link these to STRADL ID
names(STRADL_KORA)[names(STRADL_KORA) == 'ID'] <- 'meth_ID'

# Load ID linkages
ID_link_attempt <- read.csv("STRADL_DNAm_target_REM_17April2020.csv")
Meth_link <- ID_link_attempt %>% select("stradl_ID", "meth_ID")

STRADL_KORA <- merge(Meth_link, 
                     STRADL_KORA,
                     by = "meth_ID")

# rename the ID columns 
names(STRADL_KORA)[names(STRADL_KORA) == 'meth_ID_DNAm'] <- 'meth_ID'
names(STRADL_KORA)[names(STRADL_KORA) == 'stradl_ID_DNAm'] <- 'stradl_ID'


dataset_KORA <- merge(STRADL_KORA, 
                      STRADL_lifestyle_covariates,
                      by = "stradl_ID")

global_neuroimaging_dataset <- merge(dataset_KORA,
                                     STRADL_MRI,
                                     by = "stradl_ID")

# Converting multiple varibles into a factor
global_neuroimaging_dataset %<>% mutate_at(c("sex", "site", "edited", "batch"),
                                           as.factor)

### Examine global measures, add L and R hemispheres

global_neuroimaging_dataset %<>% mutate(global_cortical_surface_area =
                                          hem.lh.csa + hem.rh.csa)

global_neuroimaging_dataset %<>% mutate(global_cortical_thickness =
                                          hem.lh.ct + hem.rh.ct)

global_neuroimaging_dataset %<>% mutate(global_cortical_volume =
                                          hem.lh.cv + hem.rh.cv)

# Making a data frame with all combinations of variables
df <-
  as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1, metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
get_summary_values <- function(list) {
  metric <- list$metric
  DNAm <- list$DNAm
  tib <-
    broom::tidy(summary(
      lm(
        scale(global_neuroimaging_dataset[[metric]]) ~
          scale(st_age)
        + sex
        + site
        + batch
        + edited
        + scale(Standardised_ICV)
        + scale(global_neuroimaging_dataset[[DNAm]]),
        data = global_neuroimaging_dataset
      )
    ))[8, c(2, 3, 5)]
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
           case_when(str_detect(metric, "global") ~ "global",
                     TRUE ~ "misc")) %>%
  
  # Create a brain metric column
  mutate(
    brain_metric =
      case_when(
        ### Global
        metric == "global.cerebral.wm" ~ "global white matter",
        metric == "global.total.gm" ~ "global grey matter",
        metric == "global.wbv" ~ "total brain volume",
        #metric == "global.wbv.vent" ~ "global ventricles volume",
        metric == "global_cortical_volume" ~ "global cortical volume",
        metric == "global_cortical_thickness" ~ "global cortical thickness",
        metric == "global_cortical_surface_area" ~ "global cortical surface area",
        
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


plot_global_methylation <- newdf

##### Proteins

protein_data_full <- read.csv("prot_file_150621.csv", check.names=FALSE)
annotation_data <- read.csv("Annotation_data.csv", check.names=FALSE)

## ---------------------------

# clean protein data
newnames <- c()
for(colname in names(protein_data_full)){
  if(colname %in% annotation_data$SeqId){
    egname <- annotation_data[annotation_data$SeqId == colname,]$`Entrez Gene Name`
    newnames <- append(newnames, egname)
  } else {
    newnames <- append(newnames, colname)
  }
}

print(newnames)

names(protein_data_full) <- newnames

# Load main data
names(protein_data_full)[names(protein_data_full) == 'SampleId'] <- 'stradl_ID'
protein_data_full <- protein_data_full[ -c(1:6,8:32) ]
protein_data_only <- protein_data_full[c(1:4236)]


select_proteins <- protein_data_only %>% 
  select(
    stradl_ID,
    ACY1    ,
    ADAMTS13  ,
    ADIPOQ,
    AFM       ,
    B2M      ,
    BCAM     ,
    BMP1   ,
    `C4A|C4B` ,
    C5      ,
    C9        ,
    CCL17     ,
    CCL18     ,
    CCL21    ,
    CCL22    ,
    CCL25  ,
    CD163   ,
    CD209   ,
    CD48      ,
    CD5L      ,
    CHIT1     ,
    CLEC11A  ,
    CLEC11A.1,
    CNTN4  ,
    CRP     ,
    CXCL10  ,
    CXCL11    ,
    EDA       ,
    ENPP7     ,
    ESM1     ,
    F7       ,
    FAP    ,
    FCER2   ,
    FCGR3B  ,
    GHR       ,
    GNLY      ,
    GP1BA     ,
    GZMA     ,
    HGFAC    ,
    ICAM5  ,
    IDUA    ,
    IGFBP1  ,
    IGFBP4    ,
    IL19      ,
    INSR      ,
    LGALS3BP ,
    LGALS4   ,
    `LTA|LTB`,
    LTF     ,
    LY9     ,
    LYZ       ,
    MIA       ,
    MMP1     ,
    MMP12    ,
    MMP2   ,
    MMP9    ,
    MPL     ,
    MPO       ,
    MRC2      ,
    MST1      ,
    NCAM1    ,
    NOTCH1   ,
    NTRK3  ,
    OMD     ,
    PAPPA   ,
    PIGR      ,
    PRSS2     ,
    RARRES2   ,
    RETN     ,
    S100A9   ,
    SELE   ,
    SELL    ,
    SEMA3E  ,
    SERPINA3  ,
    SERPIND1  ,
    SHBG      ,
    SLITRK5  ,
    SPOCK2   ,
    STC1   ,
    THBS2   ,
    TNFRSF17  ,
    TNFRSF1B  ,
    TPSB2     ,
    VCAM1    ,
    WFIKKN2
  )

select_proteins <- merge(STRADL_lifestyle_covariates,
                         select_proteins,
                         by = "stradl_ID")

global_neuroimaging_dataset <- merge(select_proteins,
                                     STRADL_MRI,
                                     by = "stradl_ID")

# Converting multiple varibles into a factor
global_neuroimaging_dataset %<>% mutate_at(c("sex", "site", "edited", "batch"),
                                           as.factor)



### Examine global measures, add L and R hemispheres

global_neuroimaging_dataset %<>% mutate(global_cortical_surface_area =
                                          hem.lh.csa + hem.rh.csa)

global_neuroimaging_dataset %<>% mutate(global_cortical_thickness =
                                          hem.lh.ct + hem.rh.ct)

global_neuroimaging_dataset %<>% mutate(global_cortical_volume =
                                          hem.lh.cv + hem.rh.cv)

# Making a data frame with all combinations of variables
df <-
  as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1, metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
get_summary_values <- function(list) {
  metric <- list$metric
  DNAm <- list$DNAm
  tib <-
    broom::tidy(summary(
      lm(
        scale(global_neuroimaging_dataset[[metric]]) ~
          scale(st_age)
        + sex
        + site
        + batch
        + edited
        + scale(Standardised_ICV)
        + scale(global_neuroimaging_dataset[[DNAm]]),
        data = global_neuroimaging_dataset
      )
    ))[8, c(2, 3, 5)]
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
           case_when(str_detect(metric, "global") ~ "global",
                     TRUE ~ "misc")) %>%
  
  # Create a brain metric column
  mutate(
    brain_metric =
      case_when(
        ### Global
        metric == "global.cerebral.wm" ~ "global white matter",
        metric == "global.total.gm" ~ "global grey matter",
        metric == "global.wbv" ~ "total brain volume",
        #metric == "global.wbv.vent" ~ "global ventricles volume",
        metric == "global_cortical_volume" ~ "global cortical volume",
        metric == "global_cortical_thickness" ~ "global cortical thickness",
        metric == "global_cortical_surface_area" ~ "global cortical surface area",
        
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
  mutate(omic_type = "protein")


plot_global_proteins <- newdf

plot <- rbind(plot_global_methylation, plot_global_proteins)


###### PLOT

#plot2 <- plot 

plot2 <- plot %>% filter(brain_metric == "global grey matter")

plot2$omic_type <- as.factor(plot2$omic_type)

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
           y = estimate,
           col = reorder(DNAm,
                         (estimate)),
           shape = omic_type,
           group = omic_type
       )
) +
  
  geom_point(position = position_dodge(width = 0.9),
             size = 2.2,
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
                                1)) +
  facet_wrap( omic_type ~ brain_metric,
              #scales="free"
              scales = "free_y",
              nrow = 1) +
  facetSettings
#+
#  scale_alpha(range=c(0.8, 1))





