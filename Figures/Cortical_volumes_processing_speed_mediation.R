## ---------------------------
##
## Script Purpose: Mediation of [DNAm-processing speed association] via regional cortical volumes
##                
## ---------------------------

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

# Converting multiple varibles into a factor
Neuroimaging_DNAm %<>% mutate_at(c("sex", "site", "edited", "batch"),
                                 as.factor)

SEM_DATA <- Neuroimaging_DNAm

###### select only instances where [global cortical volume] was a significant mediator

FULL_DNAm_list <- c("MMP12", "PIGR", "SKR3", "RARRES2", "THBS2", "NOTCH1", "NTRK3", "SEMA3E", "NCAM1")


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



df <-
  as.data.frame(expand.grid(FULL_DNAm_list, FULL_neuroimaging_list, stringsAsFactors = F)) %>%
  dplyr::rename(DNAm = Var1,
                metric = Var2)


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

######

plot <- newdf


x <- ggplot(
  plot,
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
  ylab("") +
  scale_y_continuous(limits = c(-0.05, 0.05))

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
