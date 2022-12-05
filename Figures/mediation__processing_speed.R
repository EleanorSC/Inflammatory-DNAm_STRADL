## ----------------------------# 
# Test: where are significant associations with g?
## ----------------------------#
test <- plot_global_methylation %>% filter(brain_metric == "processing speed" & pFDR <0.05)

test$DNAm
## ----------------------------# 
# The issue here is I don't think it is generating new values for each cognitive-metric-DNAm combo
## ----------------------------#

FULL_DNAm_list <- c("NTRK3_olink" ,"SIGLEC1"    , "SKR3"      ,  "CXCL9"  ,     "FGF.21",      "NTRK3" ,      "RARRES2" ,   
                    "PIGR"        ,"THBS2"      , "TPSB2"     ,  "B2M"    ,     "CRP" ,        "MMP12"  ,     "NCAM1"   ,   
                     "AFM"        , "GP1BA"     ,  "NOTCH1"   ,   "SEMA3E"  )


FULL_neuroimaging_list <- c(
  "global.cerebral.wm",
  "global.total.gm",
  "global.wbv",
  
  "global_cortical_surface_area",
  "global_cortical_thickness",
  "global_cortical_volume",
  "global_subcortical_volume",
  
  "gFA",
  "gMD", 
  
  "Fazekas_Score_Total",
  "Brain_age"
)



df <-
  as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1,
                metric = Var2)


## ----------------------------# 
#  generate indirect  effect from the mediation model [IDE]
## ----------------------------#

# Function to get summary values as tibble from a named list with the info on metric and DNAm
SEM_model_output <- function(list) {
  
  metric <- list$metric
  DNAm <- list$DNAm
  
  SEM_summary_model_test <- '

            # direct effect (a path)
            COGNITION ~ cprime*DNAM + b*BRAIN + AGE + SEX

            # mediator (b path)
            BRAIN  ~ a*DNAM + AGE + SEX + SITE + EDITS + BATCH + ICV

            # indirect effect (a*b)
            indirect_effect := a*b

            # total effect
            c := cprime + (a*b)
  '
  
  TEMP_DATA = data.frame(
    COGNITION = scale(SEM_DATA$processing_speed),
    DNAM = scale(SEM_DATA[[DNAm]]),
    BRAIN = scale(SEM_DATA[[metric]]),
    SEX = SEM_DATA$sex,
    AGE = scale(SEM_DATA$st_age),
    ICV = scale(SEM_DATA$est.icv.BAD),
    BATCH = SEM_DATA$batch,
    SITE = SEM_DATA$site,
    EDITS = SEM_DATA$edited
  )
  
  Single_SEM_raw <- sem(SEM_summary_model_test,
                        data = TEMP_DATA,
                        missing = "fiml")
  
  Single_SEM <- parameterEstimates(Single_SEM_raw,
                                   standardized = TRUE)
  
  
  Single_SEM_IDE <-
    Single_SEM %>% filter(label == "indirect_effect")
  
  
  return(Single_SEM_IDE)
  
  
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = SEM_model_output)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  select(DNAm, metric, label, est, ci.lower, ci.upper, pvalue) %>%
  group_by(metric) %>%
  #group_by(metric, DNAm) %>%
  mutate(
    pFDR = p.adjust(pvalue, method = "fdr"),
    significance = case_when(pvalue < 0.05 ~ "Yes",
                             TRUE ~ "No"),
    FDR_significance = case_when(pFDR < 0.05 ~ "Yes",
                                 TRUE ~ "No")
  ) 

baseline_model_IDE <- newdf

## ----------------------------# 
#  EXAMINE significants 
## ----------------------------#
examine_quick <- baseline_model_IDE %>% filter(pFDR < 0.05)
## ----------------------------# 
#  generate total effect from the mediation model [c]
## ----------------------------#
# Function to get summary values as tibble from a named list with the info on metric and DNAm
SEM_model_output <- function(list) {
  
  metric <- list$metric
  DNAm <- list$DNAm
  
  SEM_summary_model_test <- '

            # direct effect (a path)
            COGNITION ~ cprime*DNAM + b*BRAIN + AGE + SEX

            # mediator (b path)
            BRAIN  ~ a*DNAM + AGE + SEX + SITE + EDITS + BATCH + ICV

            # indirect effect (a*b)
            indirect_effect := a*b

            # total effect
            c := cprime + (a*b)
  '
  
  TEMP_DATA = data.frame(
    COGNITION = scale(SEM_DATA$processing_speed),
    DNAM = scale(SEM_DATA[[DNAm]]),
    BRAIN = scale(SEM_DATA[[metric]]),
    SEX = SEM_DATA$sex,
    AGE = scale(SEM_DATA$st_age),
    ICV = scale(SEM_DATA$est.icv.BAD),
    BATCH = SEM_DATA$batch,
    SITE = SEM_DATA$site,
    EDITS = SEM_DATA$edited
  )
  
  Single_SEM_raw <- sem(SEM_summary_model_test,
                        data = TEMP_DATA,
                        missing = "fiml")
  
  Single_SEM <- parameterEstimates(Single_SEM_raw,
                                   standardized = TRUE)
  
  
  Single_SEM_c <-
    Single_SEM %>% filter(label == "c")
  
  
  return(Single_SEM_c)
  
  
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = SEM_model_output)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  select(DNAm, metric, label, est, ci.lower, ci.upper, pvalue) %>%
  group_by(metric) %>%
 # group_by(metric, DNAm) %>%
  mutate(
    pFDR = p.adjust(pvalue, method = "fdr"),
    significance = case_when(pvalue < 0.05 ~ "Yes",
                             TRUE ~ "No"),
    FDR_significance = case_when(pFDR < 0.05 ~ "Yes",
                                 TRUE ~ "No")
  ) 

baseline_model_c <- newdf

## ----------------------------# 
#  generate cprime from the mediation model [cprime]
## ----------------------------#
# Function to get summary values as tibble from a named list with the info on metric and DNAm
SEM_model_output <- function(list) {
  
  metric <- list$metric
  DNAm <- list$DNAm
  
  SEM_summary_model_test <- '

            # direct effect (a path)
            COGNITION ~ cprime*DNAM + b*BRAIN + AGE + SEX

            # mediator (b path)
            BRAIN  ~ a*DNAM + AGE + SEX + SITE + EDITS + BATCH + ICV

            # indirect effect (a*b)
            indirect_effect := a*b

            # total effect
            c := cprime + (a*b)
  '
  
  TEMP_DATA = data.frame(
    COGNITION = scale(SEM_DATA$processing_speed),
    DNAM = scale(SEM_DATA[[DNAm]]),
    BRAIN = scale(SEM_DATA[[metric]]),
    SEX = SEM_DATA$sex,
    AGE = scale(SEM_DATA$st_age),
    ICV = scale(SEM_DATA$est.icv.BAD),
    BATCH = SEM_DATA$batch,
    SITE = SEM_DATA$site,
    EDITS = SEM_DATA$edited
  )
  
  Single_SEM_raw <- sem(SEM_summary_model_test,
                        data = TEMP_DATA,
                        missing = "fiml")
  
  Single_SEM <- parameterEstimates(Single_SEM_raw,
                                   standardized = TRUE)
  
  
  Single_SEM_cprime <-
    Single_SEM %>% filter(label == "cprime")
  
  
  return(Single_SEM_cprime)
  
  
}

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = SEM_model_output)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  select(DNAm, metric, label, est, ci.lower, ci.upper, pvalue) %>%
  group_by(metric) %>%
  #group_by(metric, DNAm) %>%
  mutate(
    pFDR = p.adjust(pvalue, method = "fdr"),
    significance = case_when(pvalue < 0.05 ~ "Yes",
                             TRUE ~ "No"),
    FDR_significance = case_when(pFDR < 0.05 ~ "Yes",
                                 TRUE ~ "No")
  ) 

baseline_model_cprime <- newdf

## ----------------------------# 
#  Pull together for supplementary table
## ----------------------------#

baseline_model_c_table <- baseline_model_c %>%
  select(DNAm, metric, est, ci.lower, ci.upper, pvalue, pFDR) %>% 
  rename(est_c = est,
         ci.lower_c = ci.lower, 
         ci.upper_c = ci.upper, 
         pvalue_c = pvalue, 
         pFDR_c = pFDR)

baseline_model_IDE_table <- baseline_model_IDE %>%
  select(est, ci.lower, ci.upper, pvalue, pFDR) %>% 
  subset(select = -c(metric)) %>%
  rename(est_IDE = est,
         ci.lower_IDE = ci.lower, 
         ci.upper_IDE = ci.upper, 
         pvalue_IDE = pvalue, 
         pFDR_IDE = pFDR) 

baseline_model_cprime_table <- baseline_model_cprime %>%
  select(est, ci.lower, ci.upper, pvalue, pFDR) %>% 
  subset(select = -c(metric)) %>%
  rename(est_cprime = est,
         ci.lower_cprime = ci.lower, 
         ci.upper_cprime = ci.upper, 
         pvalue_cprime = pvalue, 
         pFDR_cprime = pFDR)

supplementary_table_mediation_baseline_model <- cbind(baseline_model_c_table, 
                                                      baseline_model_IDE_table, 
                                                      baseline_model_cprime_table) 


supplementary_table_mediation_baseline_model_pFDR <- supplementary_table_mediation_baseline_model %>%
  
  mutate(est_c2 = abs(est_c),
         est_cprime2 = abs(est_cprime),
         difference =
           case_when(
             (est_c2 > est_cprime2) ~ (est_c2 - est_cprime2),
             (est_c2 < est_cprime2) ~ (est_cprime2 - est_c2),
             TRUE ~ 0),
         percentage_increase_decrease = ((difference/est_c2)*100),
         increase_or_decrease = 
           case_when(
             (est_c2 > est_cprime2) ~ "Attenuation",
             (est_c2 < est_cprime2) ~ "Increase",
             TRUE ~ "misc" )
  )%>%
  
  #mutate(attenuation = 
  #         ((est_c - est_cprime)/est_c)*100) %>%
  #
  group_by(DNAm, metric) %>%
  arrange(est_IDE, by_group = TRUE) %>%
  filter(pFDR_IDE < 0.05) #n.b when adding lifestyle also remove this



supplementary_table_mediation_baseline_model %<>%
  
mutate(FDR_sig =
         case_when(pFDR_IDE <0.05 ~ "Yes", 
                   TRUE ~ "No"),
       pvalue_sig =
         case_when(pvalue_IDE < 0.05 ~ "Yes", 
                   TRUE ~ "No"),
       # Create a brain metric column

         metric =
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
             metric == "Brain_age" ~ "relative brain age",
             
             TRUE ~ "misc"
           )
       )


supplementary_table_mediation_baseline_model_CRP <- supplementary_table_mediation_baseline_model %>% filter(DNAm == "CRP")


## ----------------------------# 
#  Examine PIGR for figure
## ----------------------------#

supplementary_table_mediation_baseline_model_PIGR <- supplementary_table_mediation_baseline_model %>% 
  filter(DNAm == "PIGR" & pvalue_IDE <0.05)

## ----------------------------# 
#  PLOT
## ----------------------------#


ggplot(supplementary_table_mediation_baseline_model,
       
       aes(#x = reorder(brain_metric,-(est_IDE)),
         x = reorder(metric,-(est_IDE)),
         y = est_IDE,
         shape = pvalue_sig,
         alpha =pvalue_sig,
         col = reorder(metric,-(est_IDE)),
         #col = reorder(brain_metric,-(est_IDE)),
         group = DNAm
       )
) +
  
  geom_point(position = position_dodge(width = 0.9),
             size = 1.6,
             stroke = 0.9) +
  
  geom_errorbar(
    aes(
      ymin = ci.lower_IDE,
      ymax = ci.upper_IDE
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
       x = "") +
  
  geom_hline(
    yintercept = 0,
    color = "lightgrey",
    linetype = "dashed",
    size = 0.3,
    alpha = 0.5
  ) +
 # scale_colour_manual(values = c(
 #   colorspace::sequential_hcl(20, 
 #                              palette = "SunsetDark")))  +
 # 
  scale_colour_manual(values = c(
    colorspace::sequential_hcl(15, 
                               palette = "Dark Mint")))  +
  scale_shape_manual(values = c(16, 5)) +
  facetSettings +
  facet_wrap(~DNAm) +
  scale_alpha_discrete(range = c(0.5,1))
