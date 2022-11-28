
SEM_DATA <- Neuroimaging_DNAm


## ----------------------------# 
# First, generate total effect from the mediation model [c]
## ----------------------------#

# Function to get summary values as tibble from a named list with the info on metric and DNAm
SEM_model_output <- function(list) {
  #cog_metric <- list$cog_metric
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
  mutate(outcome = "processing speed") %>%
  select(DNAm, metric, outcome, label, est, ci.lower, ci.upper, pvalue) %>%
  # Create a brain metric column
  mutate(
    brain_metric =
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
        TRUE ~ "misc"
      )
  ) %>%
  group_by(metric) %>%
  mutate(
    pFDR = p.adjust(pvalue, method = "fdr"),
    significance = case_when(pvalue < 0.05 ~ "Yes",
                             TRUE ~ "No"),
    FDR_significance = case_when(pFDR < 0.05 ~ "Yes",
                                 TRUE ~ "No")
  ) 
 
cortical_mediators_c <- newdf

## ----------------------------# 
# /END , generate total effect from the mediation model [c]
## ----------------------------#

## ----------------------------# 
# (2) , generate total indirect effect from the mediation model [IDE]
## ----------------------------#

# Function to get summary values as tibble from a named list with the info on metric and DNAm
SEM_model_output <- function(list) {
  #cog_metric <- list$cog_metric
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
  mutate(outcome = "processing speed") %>%
  select(DNAm, metric, outcome, label, est, ci.lower, ci.upper, pvalue) %>%
  # Create a brain metric column
  mutate(
    brain_metric =
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
        TRUE ~ "misc"
      )
  ) %>%
  group_by(metric) %>%
  mutate(
    pFDR = p.adjust(pvalue, method = "fdr"),
    significance = case_when(pvalue < 0.05 ~ "Yes",
                             TRUE ~ "No"),
    FDR_significance = case_when(pFDR < 0.05 ~ "Yes",
                                 TRUE ~ "No")
  ) 

cortical_mediators_IDE <- newdf

## ----------------------------# 
# \END (2) generate total indirect effect from the mediation model [IDE]
## ----------------------------#

## ----------------------------# 
# \Start (3) generate cprime effect from the mediation model [c']
## ----------------------------#
# Function to get summary values as tibble from a named list with the info on metric and DNAm
SEM_model_output <- function(list) {
  #cog_metric <- list$cog_metric
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
  mutate(outcome = "processing speed") %>%
  select(DNAm, metric, outcome, label, est, ci.lower, ci.upper, pvalue) %>%
  # Create a brain metric column
  mutate(
    brain_metric =
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
        TRUE ~ "misc"
      )
  ) %>%
  group_by(metric) %>%
  mutate(
    pFDR = p.adjust(pvalue, method = "fdr"),
    significance = case_when(pvalue < 0.05 ~ "Yes",
                             TRUE ~ "No"),
    FDR_significance = case_when(pFDR < 0.05 ~ "Yes",
                                 TRUE ~ "No")
  ) 

cortical_mediators_cprime <- newdf

## ----------------------------# 
# \END (3) generate cprime effect from the mediation model [c']
## ----------------------------#


cortical_mediators_c <- cortical_mediators_c %>%
  select(DNAm, brain_metric, est, ci.lower, ci.upper, pvalue, pFDR)%>% 
  subset(select = -c(metric)) %>% 
  rename(est_c = est,
         ci.lower_c = ci.lower, 
         ci.upper_c = ci.upper, 
         pvalue_c = pvalue, 
         pFDR_c = pFDR)

cortical_mediators_IDE <- cortical_mediators_IDE %>%
  select(est, ci.lower, ci.upper, pvalue, pFDR) %>% 
  subset(select = -c(metric)) %>%
  rename(est_IDE = est,
         ci.lower_IDE = ci.lower, 
         ci.upper_IDE = ci.upper, 
         pvalue_IDE = pvalue, 
         pFDR_IDE = pFDR) 

cortical_mediators_cprime <- cortical_mediators_cprime %>%
  select(est, ci.lower, ci.upper, pvalue, pFDR) %>% 
  subset(select = -c(metric)) %>%
  rename(est_cprime = est,
         ci.lower_cprime = ci.lower, 
         ci.upper_cprime = ci.upper, 
         pvalue_cprime = pvalue, 
         pFDR_cprime = pFDR)

supplementary_table_mediation <- cbind(cortical_mediators_c, cortical_mediators_IDE, cortical_mediators_cprime) %>%
  mutate(attenuation = 
           ((est_c - est_cprime)/est_c)*100) %>%
  group_by(DNAm, brain_metric) %>%
  arrange(est_IDE, by_group = TRUE) %>%
  filter(pFDR_IDE < 0.05 & brain_metric == "superior temporal")


write.csv(supplementary_table_mediation, "supplementary_table_mediation.csv")
######

## ----------------------------# 
# MODEL 3 - LIFESTYLE
## ----------------------------#

phenotypes_STRADL <- read.csv("phenotypes.csv")
ID_link_attempt <- read.csv("STRADL_DNAm_target_REM_17April2020.csv")
GS_link <- ID_link_attempt %>% select("stradl_ID", "GS_id")
names(phenotypes_STRADL)[names(phenotypes_STRADL) == 'id'] <- 'GS_id'
phenotypes_STRADL <- merge(GS_link, phenotypes_STRADL, by = "GS_id")


sensitivity_analysis_lifestyle <- phenotypes_STRADL %>% select(stradl_ID, 
                                                               units, 
                                                               drink_status, 
                                                               hypertension_category, 
                                                               bmi, 
                                                               whr, 
                                                               body_fat)
#STRADL_main_data
### Clean up this lifestyle data
sensitivity_analysis_lifestyle %<>% 
  mutate(hypertension = case_when(
    hypertension_category == 1 ~ 1,
    hypertension_category == 2 ~ 1,
    hypertension_category == 3 ~ 1,
    TRUE ~ 0),
    CurrentDrinker = case_when(
      drink_status == 1 ~ 1,
      TRUE ~ 0))

# n = 778
Lifestyle_DNAm <- merge(PROTEOMICS_DNAm_DATA, 
                        sensitivity_analysis_lifestyle,
                        by = "stradl_ID")

#####
Neuroimaging_DNAm_lifestyle <- merge(Lifestyle_DNAm, 
                                     STRADL_FreeSurfer,
                                     by = "stradl_ID")

Neuroimaging_DNAm_lifestyle <- merge(Neuroimaging_DNAm_lifestyle,
                                     Cognitive_data_select,
                                     by = "stradl_ID")

# Converting multiple varibles into a factor
Neuroimaging_DNAm_lifestyle %<>% mutate_at(c("sex", "site", "edited", "batch",
                                             "hypertension", "CurrentSmoker", "CurrentDrinker"),
                                           as.factor)


SEM_DATA_LIFESTYLE <- Neuroimaging_DNAm_lifestyle

# Function to get summary values as tibble from a named list with the info on metric and DNAm
SEM_lifestyle_model_output <- function(list) {
  #cog_metric <- list$cog_metric
  metric <- list$metric
  DNAm <- list$DNAm
  
  SEM_summary_model_test <- '

            # direct effect (a path)
            COGNITION ~ cprime*DNAM + b*BRAIN + AGE + SEX

            # mediator (b path)
            BRAIN  ~ a*DNAM + AGE + SEX + SITE + EDITS + BATCH + ICV + HYPERTENSION + SMOKER + DRINKER + BMI

            # indirect effect (a*b)
            indirect_effect := a*b

            # total effect
            c := cprime + (a*b)
  '
  
  TEMP_DATA = data.frame(
    COGNITION = scale(SEM_DATA_LIFESTYLE$processing_speed),
    DNAM = scale(SEM_DATA_LIFESTYLE[[DNAm]]),
    BRAIN = scale(SEM_DATA_LIFESTYLE[[metric]]),
    SEX = SEM_DATA_LIFESTYLE$sex,
    AGE = scale(SEM_DATA_LIFESTYLE$st_age),
    ICV = scale(SEM_DATA_LIFESTYLE$est.icv.BAD),
    BATCH = SEM_DATA_LIFESTYLE$batch,
    SITE = SEM_DATA_LIFESTYLE$site,
    EDITS = SEM_DATA_LIFESTYLE$edited,
    
    HYPERTENSION = SEM_DATA_LIFESTYLE$hypertension,
    SMOKER = SEM_DATA_LIFESTYLE$CurrentSmoker,
    DRINKER = SEM_DATA_LIFESTYLE$CurrentDrinker,
    BMI = scale(SEM_DATA_LIFESTYLE$bmi)
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
  purrr::map_dfr(.f = SEM_lifestyle_model_output)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  mutate(outcome = "processing speed") %>%
  select(DNAm, metric, outcome, label, est, ci.lower, ci.upper, pvalue) %>%
  # Create a brain metric column
  mutate(
    brain_metric =
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
        
        
        TRUE ~ "misc"
      )
  ) %>%
  group_by(metric) %>%
  mutate(
    pFDR = p.adjust(pvalue, method = "fdr"),
    significance = case_when(pvalue < 0.05 ~ "Yes",
                             TRUE ~ "No"),
    FDR_significance = case_when(pFDR < 0.05 ~ "Yes",
                                 TRUE ~ "No")
  )

######

cortical_mediators_IDE_lifestyle <- newdf

#### c
# Function to get summary values as tibble from a named list with the info on metric and DNAm
SEM_lifestyle_model_output <- function(list) {
  #cog_metric <- list$cog_metric
  metric <- list$metric
  DNAm <- list$DNAm
  
  SEM_summary_model_test <- '

            # direct effect (a path)
            COGNITION ~ cprime*DNAM + b*BRAIN + AGE + SEX

            # mediator (b path)
            BRAIN  ~ a*DNAM + AGE + SEX + SITE + EDITS + BATCH + ICV + HYPERTENSION + SMOKER + DRINKER + BMI

            # indirect effect (a*b)
            indirect_effect := a*b

            # total effect
            c := cprime + (a*b)
  '
  
  TEMP_DATA = data.frame(
    COGNITION = scale(SEM_DATA_LIFESTYLE$processing_speed),
    DNAM = scale(SEM_DATA_LIFESTYLE[[DNAm]]),
    BRAIN = scale(SEM_DATA_LIFESTYLE[[metric]]),
    SEX = SEM_DATA_LIFESTYLE$sex,
    AGE = scale(SEM_DATA_LIFESTYLE$st_age),
    ICV = scale(SEM_DATA_LIFESTYLE$est.icv.BAD),
    BATCH = SEM_DATA_LIFESTYLE$batch,
    SITE = SEM_DATA_LIFESTYLE$site,
    EDITS = SEM_DATA_LIFESTYLE$edited,
    
    HYPERTENSION = SEM_DATA_LIFESTYLE$hypertension,
    SMOKER = SEM_DATA_LIFESTYLE$CurrentSmoker,
    DRINKER = SEM_DATA_LIFESTYLE$CurrentDrinker,
    BMI = scale(SEM_DATA_LIFESTYLE$bmi)
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
  purrr::map_dfr(.f = SEM_lifestyle_model_output)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  mutate(outcome = "processing speed") %>%
  select(DNAm, metric, outcome, label, est, ci.lower, ci.upper, pvalue) %>%
  # Create a brain metric column
  mutate(
    brain_metric =
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
        
        
        TRUE ~ "misc"
      )
  ) %>%
  group_by(metric) %>%
  mutate(
    pFDR = p.adjust(pvalue, method = "fdr"),
    significance = case_when(pvalue < 0.05 ~ "Yes",
                             TRUE ~ "No"),
    FDR_significance = case_when(pFDR < 0.05 ~ "Yes",
                                 TRUE ~ "No")
  )

######

cortical_mediators_c_lifestyle <- newdf

#####
#### c
# Function to get summary values as tibble from a named list with the info on metric and DNAm
SEM_lifestyle_model_output <- function(list) {
  #cog_metric <- list$cog_metric
  metric <- list$metric
  DNAm <- list$DNAm
  
  SEM_summary_model_test <- '

            # direct effect (a path)
            COGNITION ~ cprime*DNAM + b*BRAIN + AGE + SEX

            # mediator (b path)
            BRAIN  ~ a*DNAM + AGE + SEX + SITE + EDITS + BATCH + ICV + HYPERTENSION + SMOKER + DRINKER + BMI

            # indirect effect (a*b)
            indirect_effect := a*b

            # total effect
            c := cprime + (a*b)
  '
  
  TEMP_DATA = data.frame(
    COGNITION = scale(SEM_DATA_LIFESTYLE$processing_speed),
    DNAM = scale(SEM_DATA_LIFESTYLE[[DNAm]]),
    BRAIN = scale(SEM_DATA_LIFESTYLE[[metric]]),
    SEX = SEM_DATA_LIFESTYLE$sex,
    AGE = scale(SEM_DATA_LIFESTYLE$st_age),
    ICV = scale(SEM_DATA_LIFESTYLE$est.icv.BAD),
    BATCH = SEM_DATA_LIFESTYLE$batch,
    SITE = SEM_DATA_LIFESTYLE$site,
    EDITS = SEM_DATA_LIFESTYLE$edited,
    
    HYPERTENSION = SEM_DATA_LIFESTYLE$hypertension,
    SMOKER = SEM_DATA_LIFESTYLE$CurrentSmoker,
    DRINKER = SEM_DATA_LIFESTYLE$CurrentDrinker,
    BMI = scale(SEM_DATA_LIFESTYLE$bmi)
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
  purrr::map_dfr(.f = SEM_lifestyle_model_output)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>%
  mutate(outcome = "processing speed") %>%
  select(DNAm, metric, outcome, label, est, ci.lower, ci.upper, pvalue) %>%
  # Create a brain metric column
  mutate(
    brain_metric =
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
        
        
        TRUE ~ "misc"
      )
  ) %>%
  group_by(metric) %>%
  mutate(
    pFDR = p.adjust(pvalue, method = "fdr"),
    significance = case_when(pvalue < 0.05 ~ "Yes",
                             TRUE ~ "No"),
    FDR_significance = case_when(pFDR < 0.05 ~ "Yes",
                                 TRUE ~ "No")
  )

######

cortical_mediators_cprime_lifestyle <- newdf

## ----------------------------# 
# \END (3) generate cprime effect from the mediation model [c']
## ----------------------------#


cortical_mediators_c_lifestyle <- cortical_mediators_c_lifestyle %>%
  select(DNAm, brain_metric, est, ci.lower, ci.upper, pvalue, pFDR)%>% 
  subset(select = -c(metric)) %>% 
  rename(est_c = est,
         ci.lower_c = ci.lower, 
         ci.upper_c = ci.upper, 
         pvalue_c = pvalue, 
         pFDR_c = pFDR)

cortical_mediators_IDE_lifestyle <- cortical_mediators_IDE_lifestyle %>%
  select(est, ci.lower, ci.upper, pvalue, pFDR) %>% 
  subset(select = -c(metric)) %>%
  rename(est_IDE = est,
         ci.lower_IDE = ci.lower, 
         ci.upper_IDE = ci.upper, 
         pvalue_IDE = pvalue, 
         pFDR_IDE = pFDR) 

cortical_mediators_cprime_lifestyle <- cortical_mediators_cprime_lifestyle %>%
  select(est, ci.lower, ci.upper, pvalue, pFDR) %>% 
  subset(select = -c(metric)) %>%
  rename(est_cprime = est,
         ci.lower_cprime = ci.lower, 
         ci.upper_cprime = ci.upper, 
         pvalue_cprime = pvalue, 
         pFDR_cprime = pFDR)

supplementary_table_mediation_lifestyle <- cbind(cortical_mediators_c_lifestyle, 
                                       cortical_mediators_IDE_lifestyle, 
                                       cortical_mediators_cprime_lifestyle) %>%
  mutate(attenuation = 
           ((est_c - est_cprime)/est_c)*100) %>%
  group_by(DNAm, brain_metric) %>%
  arrange(est_IDE, by_group = TRUE) %>%
  filter(pFDR_IDE < 0.05 
         #& brain_metric == "superior temporal"
         )


#write.csv(supplementary_table_mediation, "supplementary_table_mediation.csv")

###############
# PLOT
###############
## ----------------------------# 
# (2) , generate total indirect effect from the mediation model [IDE]
## ----------------------------#

# Function to get summary values as tibble from a named list with the info on metric and DNAm
SEM_model_output <- function(list) {
  #cog_metric <- list$cog_metric
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
  mutate(outcome = "processing speed") %>%
  select(DNAm, metric, outcome, label, est, ci.lower, ci.upper, pvalue) %>%
  # Create a brain metric column
  mutate(
    brain_metric =
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
        TRUE ~ "misc"
      )
  ) %>%
  group_by(metric) %>%
  mutate(
    pFDR = p.adjust(pvalue, method = "fdr"),
    significance = case_when(pvalue < 0.05 ~ "Yes",
                             TRUE ~ "No"),
    FDR_significance = case_when(pFDR < 0.05 ~ "Yes",
                                 TRUE ~ "No")
  ) 

IDE_plot <- newdf


x <- ggplot(
  IDE_plot,
  aes(
    y = est,
    x = reorder(DNAm,-est),
    color = reorder(brain_metric,-est),
    fill = reorder(brain_metric,-est),
    alpha = FDR_significance,
    group = brain_metric
  )
) +
  
  geom_col(aes(
    y = est,
    x = reorder(DNAm,-est),
    group = brain_metric
  )) +
  
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(axis.text.x = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(size = 12)) +
  
  coord_flip() +
  xlab("") +
  ylab("") 
  #scale_y_continuous(limits = c(-0.05, 0.05))

#Visualise

x +
  theme_classic() +
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
    )
  )