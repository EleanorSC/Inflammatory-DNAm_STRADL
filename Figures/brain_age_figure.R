## Brain age associations

## ---------------------------
## START plot multiple histograms 
## ----------------------------



df <- Neuroimaging_DNAm %>% select(stradl_ID,
                                   Brain_age,
                                   st_age,
                                   brain_accel) %>%
  rename(brain_age = Brain_age,
         chronological_age = st_age,
         accelerated_brain_age = brain_accel)

#Convert into long format

data_long <- reshape2::melt(df,
                            # ID variables - all the variables to keep but not split apart on
                            id.vars=c("stradl_ID"),
                            # The source columns
                            measure.vars=c("brain_age", 
                                           "chronological_age", 
                                           "accelerated_brain_age"),
                            # Name of the destination column that will identify the original
                            # column that the measurement came from
                            variable.name="age_metric",
                            value.name="measurement")


# Calculate the mean of each
## Mean and SD
msd <- na.omit(data_long) %>% 
  group_by(age_metric) %>% 
  summarise(mean= mean(measurement), 
            sd=sd(measurement),
            n = n()
  )

## ----------------------------# 
# PLOT 1
## ----------------------------#

df <- Neuroimaging_DNAm %>% select(stradl_ID,
                                   Brain_age,
                                   st_age,
                                   brain_accel) %>%
  rename(brain_age = Brain_age,
         chronological_age = st_age,
         accelerated_brain_age = brain_accel)

#Convert into long format

data_long2 <- reshape2::melt(df,
                            # ID variables - all the variables to keep but not split apart on
                            id.vars=c("stradl_ID"),
                            # The source columns
                            measure.vars=c("brain_age", 
                                           "chronological_age"),
                            # Name of the destination column that will identify the original
                            # column that the measurement came from
                            variable.name="age_metric",
                            value.name="measurement")


# Calculate the mean of each
## Mean and SD
msd <- na.omit(data_long2) %>% 
  group_by(age_metric) %>% 
  summarise(mean= mean(measurement), 
            sd=sd(measurement),
            n = n()
  )
  
ggplot(data_long2, 
       aes(measurement, fill = age_metric)
) + 
  
  geom_histogram(alpha = 0.5, 
                 aes(y = ..density..), 
                 #  binwidth = 30,
                 position = 'identity') +
  
  geom_vline(data = msd, 
             aes(xintercept = mean, 
                 color = age_metric), 
             size = 1,
             alpha = 0.6
  ) +
  
  geom_vline(data = msd, 
             aes(xintercept = mean - sd, 
                 color = age_metric), 
             size = 0.5,
             alpha = 0.6,
             linetype="dotdash") +
  
  geom_vline(data = msd, 
             aes(xintercept = mean + sd, 
                 color = age_metric), 
             size = 0.5,
             alpha = 0.6,
             linetype="dotdash") +
  
  facet_wrap(~age_metric,
            # scales = "free_y",
             nrow = 1) +
  
  labs (x = "",
        y = "") +
  
  theme_classic() +
  theme(
    axis.title.x =
      element_text(
        size = 11,
        face = "bold",
        colour = "black",
        family = "sans"
      ),
    
    strip.text = element_text(
      size = 10,
      face = "bold",
      family = "sans",
      colour = "black"
    ),
    
    
    axis.text.y = element_text(size = 9,
                               colour = "black"
                               
    ),
    
    axis.text.x = element_text(size = 9,
                               colour = "black"
                               
    ),
    
    axis.title.y = element_text(
      size = 11,
      face = "bold",
      colour = "black",
      family = "sans"
    )
  ) +
  geom_density(alpha = 0.2) +
  
  scale_colour_manual(values = c("#DCD0FF",
                                 "#3E356BFF")) +
  
  scale_fill_manual(values = c("#DCD0FF",
                               "#3E356BFF")) 
  
  
  
#  scale_colour_manual(values = c("#3E356BFF",
#                                 "#357BA2FF",
#                                 "#49C1ADFF")) +
#  
#  scale_fill_manual(values = c("#3E356BFF",
#                                 "#357BA2FF",
#                                 "#49C1ADFF"))
#

#viridis::scale_fill_viridis(discrete = TRUE,
#                            option = "G",
#                            direction = -1) +
#  
# viridis::scale_colour_viridis(discrete = TRUE,
#                                option = "G",
#                                direction = -1) +

#values = c(viridis::mako(n = 5))
## ---------------------------
## END 
## ----------------------------

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


FULL_neuroimaging_list <- c("Brain_age", "brain_accel")

# Converting multiple varibles into a factor
Neuroimaging_DNAm %<>% mutate_at(c("sex", "site", "edited", "batch"),
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
      scale(Neuroimaging_DNAm[[metric]]) ~
        scale(st_age)
      + sex
      #+ site
      #+ batch
      #+ edited
      #+ scale(est.icv.BAD)
      + scale(Neuroimaging_DNAm[[DNAm]]),
      data = Neuroimaging_DNAm
    )
  )
  
  regression_summary_tidy <- broom::tidy(regression_summary)
  regression_summary_tidy_complete <- regression_summary_tidy %>% mutate(r2 = regression_summary$r.squared)
  regression_summary_tidy_DNAm <- regression_summary_tidy_complete[4, c(2, 3, 5, 6)]
  
  return(regression_summary_tidy_DNAm)
  
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = model_output)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality = 
           case_when(
             metric == "Brain_age" ~ "brain ageing",
             metric == "brain_accel" ~ "brain ageing",
             TRUE ~ "misc")) %>%
  
  # Create a brain metric column
  mutate(
    brain_metric =
      case_when(
        metric == "Brain_age" ~ "brain age",
        metric == "brain_accel" ~ "brain age acceleration",
        TRUE ~ "misc"
      )
  ) %>%
  
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
  mutate(model = "Model 1 (n=709)")

###
plot_brain_ageing_methylation <- newdf

#### proteomics

# Create list to loop through all DNAm signatures
FULL_proteins_list <- c(
  "ACY1_proteomics"    ,
  "ADAMTS13_proteomics"  ,
  "ADIPOQ_proteomics",
  "AFM_proteomics"       ,
  "B2M_proteomics"      ,
  "BCAM_proteomics"     ,
  "BMP1_proteomics"   ,
  "C4A.C4B_proteomics" ,
  "C5_proteomics"      ,
  "C9_proteomics"       ,
  "CCL17_proteomics"     ,
  "CCL18_proteomics"   ,
  "CCL21_proteomics"    ,
  "CCL22_proteomics"  ,
  "CCL25_proteomics"  ,
  "CD163_proteomics" ,
  "CD209_proteomics"   ,
  "CD48_proteomics"      ,
  "CD5L_proteomics"      ,
  "CHIT1_proteomics"     ,
  "CLEC11A_proteomics"  ,
  "CLEC11A.1_proteomics",
  "CNTN4_proteomics"  ,
  "CRP_proteomics"     ,
  "CXCL10_proteomics"  ,
  "CXCL11_proteomics"    ,
  "EDA_proteomics"       ,
  "ENPP7_proteomics"     ,
  "ESM1_proteomics"     ,
  "F7_proteomics"       ,
  "FAP_proteomics"    ,
  "FCER2_proteomics"   ,
  "FCGR3B_proteomics"  ,
  "GHR_proteomics"       ,
  "GZMA_proteomics" ,
  "GNLY_proteomics"      ,
  "GP1BA_proteomics"     ,
  "HGFAC_proteomics"    ,
  "ICAM5_proteomics"  ,
  "IDUA_proteomics"    ,
  "IGFBP1_proteomics"  ,
  "IGFBP4_proteomics"    ,
  "IL19_proteomics"      ,
  "INSR_proteomics"      ,
  "LGALS3BP_proteomics" ,
  "LGALS4_proteomics"   ,
  "LTA.LTB_proteomics",
  "LTF_proteomics"   ,
  "LY9_proteomics"     ,
  "LYZ_proteomics"     ,
  "MIA_proteomics"       ,
  "MMP1_proteomics"     ,
  "MMP12_proteomics"    ,
  "MMP2_proteomics"   ,
  "MMP9_proteomics"    ,
  "MPL_proteomics"     ,
  "MPO_proteomics"       ,
  "MRC2_proteomics"      ,
  "MST1_proteomics"      ,
  "NCAM1_proteomics"    ,
  "NOTCH1_proteomics"   ,
  "NTRK3_proteomics"  ,
  "OMD_proteomics"     ,
  "PAPPA_proteomics"   ,
  "PIGR_proteomics"      ,
  "PRSS2_proteomics"     ,
  "RARRES2_proteomics"   ,
  "RETN_proteomics"     ,
  "S100A9_proteomics"   ,
  "SELE_proteomics"   ,
  "SELL_proteomics"   ,
  "SEMA3E_proteomics"  ,
  "SERPINA3_proteomics"  ,
  "SERPIND1_proteomics" ,
  "SHBG_proteomics"      ,
  "SLITRK5_proteomics"  ,
  "SPOCK2_proteomics"   ,
  "STC1_proteomics"   ,
  "THBS2_proteomics"   ,
  "TNFRSF17_proteomics"  ,
  "TNFRSF1B_proteomics" ,
  "TPSB2_proteomics"     ,
  "VCAM1_proteomics"    ,
  "WFIKKN2_proteomics",
  
  # LBC trained
  
  "CRTAM_proteomics"    ,
  "EZR_proteomics"      ,
  "NMNAT1_proteomics"   ,
  "SMPD1_proteomics"    ,
  "CCL11_proteomics"  ,
  "CXCL9_proteomics"    ,
  "HGF_proteomics"      ,
  "OSM_proteomics"      ,
  "VEGFA_proteomics" )

#  FcRL2    , #NO PROTEIN EQUIVALENT
#  G.CSF    , #NO PROTEIN EQUIVALENT
#  GDF.8    , #NO PROTEIN EQUIVALENT
#  N.CDase  , #NO PROTEIN EQUIVALENT
#  NEP      , #NO PROTEIN EQUIVALENT
#  SIGLEC1  , #NO PROTEIN EQUIVALENT
#  SKR3     , #NO PROTEIN EQUIVALENT
#  CD6      , #NO PROTEIN EQUIVALENT
#  EN.RAGE  , #NO PROTEIN EQUIVALENT
#  FGF.21   , #NO PROTEIN EQUIVALENT
#  TGF.alpha, #NO PROTEIN EQUIVALENT

# Making a data frame with all combinations of variables
df <-
  as.data.frame(expand.grid(FULL_proteins_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1,
                metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
model_output <- function(list) {
  
  metric <- list$metric
  DNAm <- list$DNAm
  
  regression_summary <- summary(
    lm(
      scale(Neuroimaging_DNAm[[metric]]) ~
        scale(st_age)
      + sex
      #+ site
      #+ batch
      #+ edited
      #+ scale(est.icv.BAD)
      + scale(Neuroimaging_DNAm[[DNAm]]),
      data = Neuroimaging_DNAm
    )
  )
  
  regression_summary_tidy <- broom::tidy(regression_summary)
  regression_summary_tidy_complete <- regression_summary_tidy %>% mutate(r2 = regression_summary$r.squared)
  regression_summary_tidy_DNAm <- regression_summary_tidy_complete[4, c(2, 3, 5, 6)]
  
  return(regression_summary_tidy_DNAm)
  
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = model_output)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality = 
           case_when(
             metric == "Brain_age" ~ "brain ageing",
             metric == "brain_accel" ~ "brain ageing",
             TRUE ~ "misc")) %>%
  
  # Create a brain metric column
  mutate(
    brain_metric =
      case_when(
        metric == "Brain_age" ~ "brain age",
        metric == "brain_accel" ~ "brain age acceleration",
        TRUE ~ "misc"
      )
  ) %>%
  
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
  mutate(omic_type = "proteomics") %>%
  mutate(model = "Model 1 (n=709)")

###
plot_brain_ageing_proteins <- newdf

# ----------------------------#
# LIFESTYLE
# ----------------------------#

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
      scale(Neuroimaging_DNAm_lifestyle[[metric]]) ~
        scale(st_age)
      + sex

      + hypertension
      + CurrentSmoker
      + CurrentDrinker
      + scale(bmi)
      
      + scale(Neuroimaging_DNAm_lifestyle[[DNAm]]),
      data = Neuroimaging_DNAm_lifestyle
    )
  )
  
  regression_summary_tidy <- broom::tidy(regression_summary)
  regression_summary_tidy_complete <- regression_summary_tidy %>% mutate(r2 = regression_summary$r.squared)
  regression_summary_tidy_DNAm <- regression_summary_tidy_complete[8, c(2, 3, 5, 6)]
  
  return(regression_summary_tidy_DNAm)
  
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = model_output)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality = 
           case_when(
             metric == "Brain_age" ~ "brain ageing",
             metric == "brain_accel" ~ "brain ageing",
             TRUE ~ "misc")) %>%
  
  # Create a brain metric column
  mutate(
    brain_metric =
      case_when(
        metric == "Brain_age" ~ "brain age",
        metric == "brain_accel" ~ "brain age acceleration",
        TRUE ~ "misc"
      )
  ) %>%
  
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
  mutate(model = "Model 3 lifestyle (n=709)")

###
plot_brain_ageing_methylation_lifestyle <- newdf

df <-
  as.data.frame(expand.grid(FULL_proteins_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1,
                metric = Var2)

# Function to get summary values as tibble from a named list with the info on metric and DNAm
model_output <- function(list) {
  
  metric <- list$metric
  DNAm <- list$DNAm
  
  regression_summary <- summary(
    lm(
      scale(Neuroimaging_DNAm_lifestyle[[metric]]) ~
        scale(st_age)
      + sex

      + hypertension
      + CurrentSmoker
      + CurrentDrinker
      + scale(bmi)
      
      + scale(Neuroimaging_DNAm_lifestyle[[DNAm]]),
      data = Neuroimaging_DNAm_lifestyle
    )
  )
  
  regression_summary_tidy <- broom::tidy(regression_summary)
  regression_summary_tidy_complete <- regression_summary_tidy %>% mutate(r2 = regression_summary$r.squared)
  regression_summary_tidy_DNAm <- regression_summary_tidy_complete[8, c(2, 3, 5, 6)]
  
  return(regression_summary_tidy_DNAm)
  
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = model_output)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  
  # Create a neurimaging modality column
  mutate(modality = 
           case_when(
             metric == "Brain_age" ~ "brain ageing",
             metric == "brain_accel" ~ "brain ageing",
             TRUE ~ "misc")) %>%
  
  # Create a brain metric column
  mutate(
    brain_metric =
      case_when(
        metric == "Brain_age" ~ "brain age",
        metric == "brain_accel" ~ "brain age acceleration",
        TRUE ~ "misc"
      )
  ) %>%
  
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
  mutate(omic_type = "proteomics") %>%
  mutate(model = "Model 3 lifestyle (n=709)")

###
plot_brain_ageing_proteins_lifestyle <- newdf

# ----------------------------#
# Supplementary table plot
# ----------------------------#

test1 <- plot_brain_ageing_methylation %>% filter(FDR_significance=="Yes" & brain_metric == "brain age")

test <- plot_brain_ageing_methylation_lifestyle %>% filter(FDR_significance=="Yes" & brain_metric == "brain age")

## ----------------------------#
# Combining both sets of models to see percentage attenuation
## ----------------------------#

plot_brain_ageing_methylation_lifestyle_table <- plot_brain_ageing_methylation_lifestyle %>%
  select(estimate, std.error, r2, p.value, pFDR) %>%
  rename(
    estimate_lifestyle = estimate,
    r2_lifestyle = r2,
    p.value_lifestyle = p.value,
    pFDR_lifestyle = pFDR,
    std.error_lifestyle = std.error
  )

plot_brain_ageing_methylation_lifestyle_table <- subset(plot_brain_ageing_methylation_lifestyle_table, select = -c(brain_metric, modality))

table_new <- cbind(plot_brain_ageing_methylation, plot_brain_ageing_methylation_lifestyle_table)

table_new %<>% mutate(percentage_increase_decrease =
                        100 * (estimate - estimate_lifestyle)) %>%
  
  filter(brain_metric == "brain age") %>%
  
  group_by(DNAm) %>%
  
  arrange(desc(estimate),
          by_group = TRUE) %>%
  
  mutate(
    CI_lower =
      estimate - (1.96 * std.error),
    CI_upper =
      estimate + (1.96 * std.error),
    
    CI_lower_lifestyle =
      estimate_lifestyle - (1.96 * std.error_lifestyle),
    CI_upper_lifestyle =
      estimate_lifestyle + (1.96 * std.error_lifestyle),
    
    better_model =
      case_when(r2_lifestyle > r2 ~ "Yes",
                TRUE ~ "No")) %>%
  
  filter(FDR_significance == "Yes") %>%
  
  select(
    brain_metric,
    DNAm,
    # Model 1
    estimate,
    CI_lower,
    CI_upper,
    p.value,
    pFDR,
    # Model 2 
    estimate_lifestyle,
    CI_lower_lifestyle,
    CI_upper_lifestyle,
    p.value_lifestyle,
    pFDR_lifestyle,
    # R2
    r2,
    r2_lifestyle,
    better_model,
    percentage_increase_decrease
  ) 


# ----------------------------#
# write to .csv
# ----------------------------#
write.csv(table_new, "brain_age_methylation_models.csv")


mean(abs(table_new$percentage_increase_decrease))

test <- table_new %>% filter(pFDR_lifestyle < 0.05)
# ----------------------------#
# EFFECT SIZE RANGE PLOT
# ----------------------------#






plot4 <- rbind(plot_brain_ageing_methylation,
               plot_brain_ageing_proteins) %>% 
  mutate(DNAm =
           stringr::str_remove(DNAm, "_proteomics")) %>% 
  mutate(DNAm = stringr::str_replace(DNAm, "_olink", ".x")) %>% 
  #filter(significance == "Yes") %>%
  filter(brain_metric == "brain age")

# ----------------------------#
# POOR BRAIN HEALTH
# We only want cases of DNAm significant, aka n = 26
# ----------------------------#


plot2 <- rbind(plot_brain_ageing_methylation,
               plot_brain_ageing_proteins) %>%
  filter(brain_metric == "brain age" & estimate > 0 & significance == "Yes" & omic_type == "DNAm")%>%
  arrange(estimate)

sig_DNAms <- c(plot2$DNAm)

plot4 <- rbind(plot_brain_ageing_methylation,
               plot_brain_ageing_proteins) %>% 
  mutate(DNAm =
           stringr::str_remove(DNAm, "_proteomics")) %>% 
  mutate(DNAm = stringr::str_replace(DNAm, "_olink", ".x")) %>% 
  filter(brain_metric == "brain age" & estimate > 0)


##### REORDER by effect size

plot4 <- plot4 %>% filter(
                          DNAm == "GHR"          | DNAm == "MMP.1_olink" | DNAm == "CCL17" |DNAm == "MST1"    |DNAm==    "MMP1"     | DNAm ==    "MMP9"  | DNAm ==       "TPSB2"  |    
                          DNAm == "CXCL10_olink" | DNAm == "SIGLEC1"     | DNAm == "ICAM5" |DNAm == "SERPIND1"|DNAm==    "LGALS3BP" | DNAm ==    "VEGFA" | DNAm ==       "CRP"    |     
                          DNAm == "SKR3"         | DNAm == "CCL18"       | DNAm == "THBS2" |DNAm == "NEP"     |DNAm==    "SELE"     | DNAm ==    "PIGR"  | DNAm ==       "MMP12"  |      
                          DNAm == "STC1"         | DNAm == "FGF.21"      | DNAm == "ACY1"  |DNAm == "PRSS2"   |DNAm==    "IGFBP4"
                          )

plot4 %<>% mutate(DNAm = factor(DNAm,
                                        levels = sig_DNAms))


facetSettings <-
  theme(strip.background = element_rect(
    fill = "#EAE3F2",
    colour = "black",
    size = 1
  ))

ggplot(plot4,
       
       aes(
         x = reorder(DNAm, (estimate)),
         y = estimate,
         colour = omic_type,
         alpha = significance, 
         group = omic_type
         
       )) +
  
  geom_point(
    aes(
      #col = brain_metric,
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
          vjust = 0.5,
          hjust = 1,
          size = 6,
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
  
#  facet_wrap( ~ omic_type 
#             # scales = "free_y"
#              
#  ) +
  
  facetSettings +

  scale_shape_manual(values = c(
    16,
    5))

# ----------------------------#
# END EFFECT SIZE RANGE PLOT
# ----------------------------#

# ----------------------------#
# BARPLOT
# ----------------------------#

####### poor brain ageing


plot2 <- rbind(plot_brain_ageing_methylation,
               plot_brain_ageing_proteins) %>%
  filter(brain_metric == "brain age" & estimate > 0 & FDR_significance == "Yes" & omic_type == "DNAm")%>%
  arrange(estimate)

sig_DNAms <- c(plot2$DNAm)


plot4 <- rbind(plot_brain_ageing_methylation,
               plot_brain_ageing_proteins) %>% 
  
  mutate(DNAm =
           stringr::str_remove(DNAm, "_proteomics")) %>% 
  mutate(DNAm = stringr::str_replace(DNAm, "_olink", ".x")) %>% 
  filter(brain_metric == "brain age" & estimate > 0) %>% 
  filter(
    DNAm ==       "TPSB2"  |    
      DNAm == "CXCL10.x" | 
      DNAm == "SIGLEC1"     | 
      DNAm == "ICAM5" |
      DNAm == "SERPIND1"|
      DNAm==    "LGALS3BP" | 
      DNAm ==    "VEGFA" | 
      DNAm ==       "CRP"    |     
      DNAm == "SKR3"         | 
      DNAm == "CCL18"       | 
      DNAm == "THBS2" |
      DNAm == "NEP"     |
      DNAm==    "SELE"     | 
      DNAm ==    "PIGR"  | 
      DNAm ==       "MMP12"  |      
      DNAm == "STC1"         | 
      DNAm == "FGF.21"      | 
      DNAm == "ACY1"  |
      DNAm == "PRSS2"   |
      DNAm==    "IGFBP4"
  )

plot4 %<>% mutate(my_alpha = case_when(omic_type == "proteomics" ~
                                         0.3, 
                                       TRUE ~ 1))



order_levels <- c("TPSB2" ,      "CXCL10.x" ,
                  "SIGLEC1"    ,  "ICAM5"    ,    
                  "SERPIND1"   ,  "LGALS3BP"  ,   "VEGFA"      , 
                  "CRP"        , 
                  "SKR3"  ,   
                  "CCL18"    ,  
                  "THBS2"  ,    
                  "NEP"  ,     
                  "SELE"   ,  
                  "PIGR"  ,  
                  "MMP12"  ,    
                  "STC1",
                  "FGF.21"  ,    
                  "ACY1" ,       "PRSS2"   ,     "IGFBP4")

plot4 %<>% mutate(DNAm = factor(DNAm,
                                levels = rev(order_levels)
)
)

plot5 <- plot4 %>% arrange(estimate)

ggplot(plot4,
       
       aes(
         x = DNAm,
         y = estimate,
         group = omic_type
         
       )) +
  
  geom_col(
    aes(x = DNAm,
        y = estimate,
        fill = DNAm,
        alpha = plot4$my_alpha, 
    ),
    colour = "black",
    position=position_dodge2(preserve = "single")
  ) +
  
  
  geom_errorbar(data =plot4, 
                aes(x = DNAm, 
                    ymin = estimate + std.error*1.96,
                    ymax = estimate
                ),
                position = position_dodge2(.9, 
                                           padding = .6, 
                                           preserve = "single"),
                size = 0.6
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
          angle = 45,
          hjust = 1,
          size = 6,
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
  
  facetSettings +
  scale_fill_manual(values = c(
    colorspace::sequential_hcl(26, 
                               palette = "Purple-Blu")
  ))

####### good brain ageing


plot2 <- rbind(plot_brain_ageing_methylation,
               plot_brain_ageing_proteins) %>%
  filter(brain_metric == "brain age" & estimate < 0 & significance == "Yes" & omic_type == "DNAm")%>%
  arrange(estimate)

sig_DNAms <- c(plot2$DNAm)


plot4 <- rbind(plot_brain_ageing_methylation,
               plot_brain_ageing_proteins) %>% 
  
  mutate(DNAm =
           stringr::str_remove(DNAm, "_proteomics")) %>% 
  mutate(DNAm = stringr::str_replace(DNAm, "_olink", ".x")) %>% 
  filter(brain_metric == "brain age" & estimate < 0) %>% 
filter( DNAm =="CNTN4"  | DNAm ==     "NCAM1"| DNAm == "SHBG"   | DNAm ==     "IGFBP1"| DNAm ==      "WFIKKN2"|DNAm ==     "NOTCH1"|DNAm ==      "NTRK3"       |DNAm =="OMD"    | DNAm==    "SLITRK5" |   
        DNAm =="SEMA3E"| DNAm ==      "MIA"   | DNAm == "SELL"    | DNAm ==    "GP1BA"  | DNAm ==     "CD209"     |DNAm ==  "FCER2"     |DNAm ==  "NTRK3_olink" |DNAm =="ADAMTS13" |DNAm==   "C4A.C4B" )

plot4 %<>% mutate(my_alpha = case_when(omic_type == "proteomics" ~
                                         0.3, 
                                       TRUE ~ 1)
)
plot4 %<>% mutate(DNAm = factor(DNAm,
                                levels = sig_DNAms))

plot5 <- plot4 %>% arrange(estimate)

ggplot(plot4,
       
       aes(
         x = DNAm,
         y = -estimate,
         #colour = significance, 
         group = omic_type
         
       )) +
  
  geom_col(
    aes(x = DNAm,
        y = -estimate,
        fill = DNAm,
        alpha = plot4$my_alpha, 
        
        
    ),colour = "black",
    position=position_dodge2(preserve = "single")
  ) +
  
  
  geom_errorbar(data =plot4, 
    aes(x = DNAm, 
      ymax = -estimate + std.error*1.96,
      ymin = -estimate
    ),
    position = position_dodge2(.9, 
                               padding = .6, 
                               preserve = "single"),
    size = 0.6
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
          angle = 45,
          hjust = 1,
          size = 6,
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
  
  facetSettings +
  scale_fill_manual(values = c(
    colorspace::sequential_hcl(17, 
                               palette = "Dark Mint") 
  )) +
  scale_colour_manual(values = c("grey", "black"))