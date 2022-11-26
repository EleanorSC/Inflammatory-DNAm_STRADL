### Mediation models
skimr::skim(Neuroimaging_DNAm)


#install.packages("lavaan")
#library(lavaan)

# ----------------------------#
# Examine individual DNAms
# ----------------------------#

# n = 702
test2 <- Neuroimaging_DNAm %>% drop_na(processing_speed)  

model1_TBV <-
  '
            # direct effect (a path)
            COGNITION ~ cprime*RARRES2 + b*TB + AGE + SEX
            # mediator (b path)
            TB  ~ a*RARRES2 + AGE + SEX + SITE + EDITS + BATCH + ICV
            #(c prime path)
            #COGNITION ~ b*TB
            # indirect effect (a*b)
            TBIDE := a*b
            # total effect
            c := cprime + (a*b)
            '
tmp = data.frame(
  COGNITION = scale(test2$processing_speed),
  RARRES2 = scale(test2$RARRES2),
  TB = scale(test2$global.wbv),
  SEX = test2$sex,
  AGE = scale(test2$st_age),
  ICV = scale(test2$est.icv.BAD),
  BATCH = test2$batch,
  SITE = test2$site,
  EDITS = test2$edited
)

Single_SEM <- sem(model1_TBV, data = tmp, missing = "fiml")

summary(
  Single_SEM,
  rsq = T,
  standardized = TRUE,
  fit.measures = TRUE
)

##############outputting MEDIATION results into a  table#################
table <- parameterEstimates(Single_SEM, standardized = TRUE)
#table <-performance::model_performance(Single_SEM, metrics = "all", verbose = TRUE)
table2 <- table %>% filter(label == "TBIDE" | label == "c" | label == "cprime")
table2 %<>% select(label, est, se, pvalue, ci.lower, ci.upper)

TB <- table2 %>% filter(label == "TBIDE") %>%
  mutate(DNAm = "RARRES2",
         outcome = "processing speed",
         effect = "indirect effect") %>%
  rename(mediator = label)

############## GM #################

model1_GM <-
  '
            # direct effect (a path)
            COGNITION ~ cprime*RARRES2 + b*GM + AGE + SEX
            # mediator (b path)
            GM  ~ a*RARRES2 + AGE + SEX + SITE + EDITS + BATCH + ICV
            #(c prime path)
            #COGNITION ~ b*GM
            # indirect effect (a*b)
            GMIDE := a*b
            # total effect
            c := cprime + (a*b)
            '
tmp = data.frame(
  COGNITION = scale(test2$processing_speed),
  RARRES2 = scale(test2$RARRES2),
  GM = scale(test2$global.total.gm),
  SEX = test2$sex,
  AGE = scale(test2$st_age),
  ICV = scale(test2$est.icv.BAD),
  BATCH = test2$batch,
  SITE = test2$site,
  EDITS = test2$edited
)

Single_SEM <- sem(model1_GM, data = tmp, missing = "fiml")

summary(
  Single_SEM,
  rsq = T,
  standardized = TRUE,
  fit.measures = TRUE
)

##############outputting MEDIATION results into a  table#################
table <- parameterEstimates(Single_SEM, standardized = TRUE)
#table <-performance::model_performance(Single_SEM, metrics = "all", verbose = TRUE)
table2 <- table %>% filter(label == "GMIDE" | label == "c" | label == "cprime")
table2 %<>% select(label, est, se, pvalue, ci.lower, ci.upper)

GM <- table2 %>% filter(label == "GMIDE") %>%
  mutate(DNAm = "RARRES2",
         outcome = "processing speed",
         effect = "indirect effect") %>%
  rename(mediator = label)


#Single_SEM_neuroimaging_ab <- rbind(TB, GM, WM, gFA, gMD)

# ----------------------------#
# AUTOMATE SEM with a function
# ----------------------------#
Neuroimaging_DNAm <- merge(PROTEOMICS_DNAm_DATA,
                           STRADL_FreeSurfer,
                           by = "stradl_ID")


# add in some cognitive metrics

Cognitive_data_select <- STRADL_main_data %>% 
  select(stradl_ID, digit_symbol, verbal_total, vocabulary, LM, mr_correct) %>% 
  rename(processing_speed = digit_symbol,
         executive_function = verbal_total,
         verbal_declarative_memory = LM,
         matrix_reasoning = mr_correct)


Neuroimaging_DNAm <- merge(Neuroimaging_DNAm,
                           Cognitive_data_select,
                           by = "stradl_ID")

Neuroimaging_DNAm %<>% mutate(scv.cerebellum_GM =  vol.cerebellum.lh.gm + vol.cerebellum.rh.gm,
                              scv.cerebellum_WM =  vol.cerebellum.lh.wm + vol.cerebellum.rh.wm)


Neuroimaging_DNAm %<>% mutate(global_cortical_surface_area = 
                                hem.lh.csa + hem.rh.csa,
                              global_cortical_thickness = 
                                hem.lh.ct + hem.rh.ct,
                              global_cortical_volume = 
                                hem.lh.cv + hem.rh.cv,
                              global_subcortical_volume =
                                scv.bilat.accumbens + scv.bilat.amygdala + scv.bilat.caudate +
                                scv.bilat.hippocampus + scv.bilat.pallidum + scv.bilat.putamen +
                                scv.bilat.thalamus + scv.bilat.ventraldc + scv.cerebellum_GM + scv.cerebellum_WM + vol.brainstem
)

# Invert the polarity of these measures
Neuroimaging_DNAm %<>% mutate(gFA =
                                gFA*-1,
                              gMD = 
                                gMD*-1)

# Converting multiple varibles into a factor
Neuroimaging_DNAm %<>% mutate_at(c("sex", "site", "edited", "batch"),
                                 as.factor)

test2 <- Neuroimaging_DNAm %>% drop_na(processing_speed,
                                       global_cortical_volume)  


FULL_DNAm_list <- c("RARRES2", 
                    "CRP", 
                    "IGFBP4", 
                    "PIGR")



FULL_neuroimaging_list <- c(
  
  ### Global measures
  "global.cerebral.wm",
  "global.total.gm",
  "global.wbv"
  
# "global_cortical_surface_area",
# "global_cortical_thickness",
# "global_cortical_volume"
# "global_subcortical_volume",
# 
# "gFA",
# "gMD"
# 
# "Fazekas_Score_Total"
  
)

df <-
  as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1,
                metric = Var2)


# Function to get summary values as tibble from a named list with the info on metric and DNAm
model_output <- function(list) {
  
  #cog_metric <- list$cog_metric
  metric <- list$metric
  DNAm <- list$DNAm
  
  SEM_summary_model_test <- '
  
            # direct effect (a path)
            COGNITION ~ cprime*DNAM + b*BRAIN + AGE + SEX
            
            # mediator (b path)
            BRAIN  ~ a*DNAM + AGE + SEX + SITE + EDITS + BATCH + ICV
            
            # indirect effect (a*b)
            BRAINIDE := a*b
            
            # total effect
            c := cprime + (a*b)
  '
  
  TEMP_DATA = data.frame(
    COGNITION = scale(test2$processing_speed),
    DNAM = scale(test2[[DNAm]]),
    BRAIN = scale(test2[[metric]]),
    SEX = test2$sex,
    AGE = scale(test2$st_age),
    ICV = scale(test2$est.icv.BAD),
    BATCH = test2$batch,
    SITE = test2$site,
    EDITS = test2$edited
    )
  
  Single_SEM <- parameterEstimates(sem(SEM_summary_model_test, 
                                       data = TEMP_DATA, 
                                       missing = "fiml"),
                        standardized = TRUE)
    
    
  SEM_summary_tidy <- broom::tidy(Single_SEM)
    
    
  return(SEM_summary_tidy)
}
  

# Making 3 new columns to add to data frame, 1 is estimate, 1 is std. error, 1 is p.value
newcols <- df %>%
  split(1:nrow(.)) %>%
  purrr::map_dfr(.f = model_output)

# Now bind on the 3 cols that were created
newdf <- cbind(df, newcols) %>% 
  mutate(outcome = "processing speed") %>%
  select(DNAm, metric, outcome, label, est, ci.lower, ci.upper, pvalue) %>%
  filter(label == "cprime" | label == "b"|label == "a" | label == "c" | label == "GMIDE")

