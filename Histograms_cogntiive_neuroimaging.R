## ---------------------------
##
## Script Purpose: Histogram plots of cognitive and brain metrics
##
##                
##
##                
##
##  
## ----------------------------
### Histogram plots

STRADL_cogn_full <- read_excel("STRADL_Laura_variables.xlsx",sheet=4)
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


# remove missing values
#cognition_STRADL  <- na.omit(cognition_STRADL)

dataset_KORA <- merge(STRADL_KORA, 
                      STRADL_lifestyle_covariates,
                      by = "stradl_ID")


cognition_STRADL <- merge(dataset_KORA,
                          cognition_STRADL,
                          by = "stradl_ID")


#plot multiple histograms

df <- cognition_STRADL %>% select(stradl_ID,
                                  Digit_Symbol_coding, # Processing speed
                                  Verbal_Fluency, # Executive function
                                  Mill_Hill_vocab,
                                  Matrix_Reasoning,
                                  Logical_Memory,
                                  g,
                                  gf)

#Convert into long format

data_long <- reshape2::melt(df,
                  # ID variables - all the variables to keep but not split apart on
                  id.vars=c("stradl_ID"),
                  # The source columns
                  measure.vars=c("Digit_Symbol_coding", 
                                 "Verbal_Fluency", 
                                 "Mill_Hill_vocab",
                                 "Matrix_Reasoning",
                                 "Logical_Memory",
                                 "g",
                                 "gf"),
                  # Name of the destination column that will identify the original
                  # column that the measurement came from
                  variable.name="cognitive_test",
                  value.name="measurement")


# Calculate the mean of each
## Mean and SD
msd <- na.omit(data_long) %>% 
  group_by(cognitive_test) %>% 
  summarise(mean= mean(measurement), 
            sd=sd(measurement),
            n = n()
            )

## ----------------------------# 
# PLOT 1
## ----------------------------#

ggplot(data_long, 
       aes(measurement, fill = cognitive_test)
       ) + 
  
  geom_histogram(alpha = 0.5, 
                 aes(y = ..density..), 
               #  binwidth = 30,
                 position = 'identity') +
  
  geom_vline(data = msd, 
             aes(xintercept = mean, 
                 color = cognitive_test), 
             size = 1,
             alpha = 0.6
             ) +
  
  geom_vline(data = msd, 
             aes(xintercept = mean - sd, 
                 color = cognitive_test), 
             size = 0.5,
             alpha = 0.6,
             linetype="dotdash") +
  
  geom_vline(data = msd, 
             aes(xintercept = mean + sd, 
                 color = cognitive_test), 
             size = 0.5,
             alpha = 0.6,
             linetype="dotdash") +
  
  facet_wrap(~cognitive_test,
             scales = "free",
             nrow = 2) +
  
  labs (x = "test score",
       y = "frequency") +
  
  viridis::scale_fill_viridis(discrete = TRUE,
                              option = "G",
                              direction = -1) +
  
  viridis::scale_colour_viridis(discrete = TRUE,
                              option = "G",
                              direction = -1) +
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
  geom_density(alpha = 0.2)

## ----------------------------# 
# Neuroimaging metrics histograms  
## ----------------------------#

## Next load up methylation and proteomics data
#protein_data_full <- read.csv("prot_file_150621.csv", check.names=FALSE)
#annotation_data <- read.csv("Annotation_data.csv", check.names=FALSE)
## Methylation:
#STRADL_KORA <- read.csv("STRADL_scores_trained_KORA.csv", check.names=FALSE)
#EpiScores <- read.csv("STRADL_scores_trained_LBC.csv")
### ---------------------------
#
## clean protein data
#newnames <- c()
#for(colname in names(protein_data_full)){
#  if(colname %in% annotation_data$SeqId){
#    egname <- annotation_data[annotation_data$SeqId == colname,]$`Entrez Gene Name`
#    newnames <- append(newnames, egname)
#  } else {
#    newnames <- append(newnames, colname)
#  }
#}
#
#print(newnames)
#
#names(protein_data_full) <- newnames
#
## Load main data
#names(protein_data_full)[names(protein_data_full) == 'SampleId'] <- 'stradl_ID'
#protein_data_full <- protein_data_full[ -c(1:6,8:32) ]
#protein_data_only <- protein_data_full[c(1:4236)]
#STRADL_main_data <- protein_data_full[c(1, 4237:4286)]


### Link what we need from STRADL main data into phenotypes + disease_covariates
STRADL_covariates <- STRADL_main_data %>% select(stradl_ID, sex, st_age,
                                                 APOE, apoe, Brain_age, brain_accel, CurrentSmoker,
                                                 Site, Global_GM_Volume, Cerebrum_WM_Volume,
                                                 g, gf, ethnicity,
                                                 Fazekas_Score_Total,
                                                 gFA,
                                                 gMD,
                                                 WBV_No_Ventricles,
                                                 Estimated_ICV,
                                                 CurrentSmoker,
                                                 BMI)


######Load up neuroimaging full dataset
STRADL_FreeSurfer <- read.csv("STRADL_Measures_FreeSurfer_Main.csv")
names(STRADL_FreeSurfer)[names(STRADL_FreeSurfer) == 'id'] <- 'stradl_ID'
STRADL_ICV <- read.csv("STRADL_Measures_Standardised_ICV.csv")
names(STRADL_ICV)[names(STRADL_ICV) == 'ID'] <- 'stradl_ID'


STRADL_MRI <- merge(STRADL_ICV, 
                    STRADL_FreeSurfer,
                    by = "stradl_ID")

# model H1, standard models
# DNAm + neuroimaging without ICV
# n = 709
dataset_n709 <- merge(KORA_LBC_DNAm, 
                      STRADL_FreeSurfer,
                      by = "stradl_ID")

# model H2,  models controlling for ICV
# DNAm + neuroimaging WITH ICV
# n = 655
dataset_n655 <- merge(KORA_LBC_DNAm, 
                      STRADL_MRI,
                      by = "stradl_ID")

#n709_noICV <- merge(STRADL_lifestyle_covariates, dataset_n709, by = "stradl_ID") 

n709_noICV <- merge(STRADL_covariates, dataset_n709, by = "stradl_ID") 

ICV_corrected <- merge(STRADL_lifestyle_covariates, dataset_n655, by = "stradl_ID") 



### Examine global measures, add L and R hemispheres

n709_noICV %<>% mutate(global_cortical_surface_area = hem.lh.csa + hem.rh.csa,
                       global_cortical_thickness = hem.lh.ct + hem.rh.ct,
                       global_cortical_volume = hem.lh.cv + hem.rh.cv) 



#plot multiple histograms

df <- n709_noICV %>% select(stradl_ID,
                            # Global_GM_Volume, 
                            # Cerebrum_WM_Volume,
                            Fazekas_Score_Total,
                            gFA,
                            gMD,
                            Brain_age,
                            brain_accel,
                            # WBV_No_Ventricles,
                            # Estimated_ICV,
                            
                            global.total.gm,
                            global.cerebral.wm,
                            global.wbv,
                            global_cortical_volume,
                            global_cortical_thickness,
                            global_cortical_surface_area,
                            est.icv.BAD
)

skimr::skim(df)

#Convert into long format

data_long <- reshape2::melt(df,
                            # ID variables - all the variables to keep but not split apart on
                            id.vars=c("stradl_ID"),
                            
                            # The source columns
                            measure.vars=c(
                              # "Global_GM_Volume", 
                              # "Cerebrum_WM_Volume",
                              "Fazekas_Score_Total",
                              "gFA",
                              "gMD",
                              "Brain_age",
                              "brain_accel",
                              # "WBV_No_Ventricles",
                              # "Estimated_ICV",
                              
                              "global.total.gm",
                              "global.cerebral.wm",
                              "global.wbv",
                              "global_cortical_volume",
                              "global_cortical_thickness",
                              "global_cortical_surface_area",
                              "est.icv.BAD"),
                            
                            # Name of the destination column that will identify the original
                            # column that the measurement came from
                            variable.name="neuroimaging_metric",
                            value.name="measurement")


# Calculate the mean of each
## Mean and SD
msd <- na.omit(data_long) %>% 
  group_by(neuroimaging_metric) %>% 
  summarise(mean= mean(measurement), 
            sd=sd(measurement),
            n = n()
  )

## ----------------------------# 
# PLOT 3 - HISTOGRAM TRADITIONAL 
## ----------------------------#

ggplot(data_long,
       aes(measurement, 
           fill = neuroimaging_metric)
) +
  
  geom_histogram(alpha = 0.5,
                 aes(y = ..density..),
                 #  binwidth = 30,
                 position = 'identity') +
  
  geom_vline(
    data = msd,
    aes(xintercept = mean,
        color = neuroimaging_metric),
    size = 1,
    alpha = 0.6
  ) +
  
  geom_vline(
    data = msd,
    aes(xintercept = mean - sd,
        color = neuroimaging_metric),
    size = 0.5,
    alpha = 0.6,
    linetype = "dotdash"
  ) +
  
  geom_vline(
    data = msd,
    aes(xintercept = mean + sd,
        color = neuroimaging_metric),
    size = 0.5,
    alpha = 0.6,
    linetype = "dotdash"
  ) +
  
  facet_wrap( ~ neuroimaging_metric,
              scales = "free",
              nrow = 3) +
  
  labs (x = "units",
        y = "frequency") +
  
  viridis::scale_fill_viridis(discrete = TRUE,
                              option = "D") +
  
  viridis::scale_colour_viridis(discrete = TRUE,
                                option = "D") +
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
      size = 8,
      face = "bold",
      family = "sans",
      colour = "black"
    ),
    
    axis.text.y = element_text(size = 9,
                               colour = "black"),
    
    axis.text.x = element_text(size = 6,
                               colour = "black"),
    
    axis.title.y = element_text(
      size = 11,
      face = "bold",
      colour = "black",
      family = "sans"
    )
  ) +
  geom_density(alpha = 0.2) 


##### The same but for Standardised ICV
# n = 655
dataset_n655 <- merge(KORA_LBC_DNAm, 
                      STRADL_MRI,
                      by = "stradl_ID")

ICV_corrected <- merge(STRADL_covariates, 
                       dataset_n655, by = "stradl_ID") 



### Examine global measures, add L and R hemispheres

ICV_corrected  %<>% mutate(global_cortical_surface_area = hem.lh.csa + hem.rh.csa,
                           global_cortical_thickness = hem.lh.ct + hem.rh.ct,
                           global_cortical_volume = hem.lh.cv + hem.rh.cv) 



#plot multiple histograms

df <- ICV_corrected  %>% select(stradl_ID,
                                # Global_GM_Volume, 
                                # Cerebrum_WM_Volume,
                                Fazekas_Score_Total,
                                gFA,
                                gMD,
                                Brain_age,
                                brain_accel,
                                # WBV_No_Ventricles,
                                # Estimated_ICV,
                                
                                global.total.gm,
                                global.cerebral.wm,
                                global.wbv,
                                global_cortical_volume,
                                global_cortical_thickness,
                                global_cortical_surface_area,
                                #  est.icv.BAD,
                                Standardised_ICV
)

skimr::skim(df)

#Convert into long format

data_long <- reshape2::melt(df,
                            # ID variables - all the variables to keep but not split apart on
                            id.vars=c("stradl_ID"),
                            
                            # The source columns
                            measure.vars=c(
                              # "Global_GM_Volume", 
                              # "Cerebrum_WM_Volume",
                              "Fazekas_Score_Total",
                              "gFA",
                              "gMD",
                              "Brain_age",
                              "brain_accel",
                              # "WBV_No_Ventricles",
                              # "Estimated_ICV",
                              
                              "global.total.gm",
                              "global.cerebral.wm",
                              "global.wbv",
                              "global_cortical_volume",
                              "global_cortical_thickness",
                              "global_cortical_surface_area",
                              # "est.icv.BAD",
                              "Standardised_ICV"),
                            
                            # Name of the destination column that will identify the original
                            # column that the measurement came from
                            variable.name="neuroimaging_metric",
                            value.name="measurement")


# Calculate the mean of each
## Mean and SD
msd <- na.omit(data_long) %>% 
  group_by(neuroimaging_metric) %>% 
  summarise(mean= mean(measurement), 
            sd=sd(measurement),
            n = n()
  )

## ----------------------------# 
# PLOT 3 - HISTOGRAM TRADITIONAL 
## ----------------------------#

ggplot(data_long,
       aes(measurement, 
           fill = neuroimaging_metric)
) +
  
  geom_histogram(alpha = 0.5,
                 aes(y = ..density..),
                 #  binwidth = 30,
                 position = 'identity') +
  
  geom_vline(
    data = msd,
    aes(xintercept = mean,
        color = neuroimaging_metric),
    size = 1,
    alpha = 0.6
  ) +
  
  geom_vline(
    data = msd,
    aes(xintercept = mean - sd,
        color = neuroimaging_metric),
    size = 0.5,
    alpha = 0.6,
    linetype = "dotdash"
  ) +
  
  geom_vline(
    data = msd,
    aes(xintercept = mean + sd,
        color = neuroimaging_metric),
    size = 0.5,
    alpha = 0.6,
    linetype = "dotdash"
  ) +
  
  facet_wrap( ~ neuroimaging_metric,
              scales = "free",
              nrow = 3) +
  
  labs (x = "units",
        y = "frequency") +
  
  viridis::scale_fill_viridis(discrete = TRUE,
                              option = "D") +
  
  viridis::scale_colour_viridis(discrete = TRUE,
                                option = "D") +
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
      size = 8,
      face = "bold",
      family = "sans",
      colour = "black"
    ),
    
    axis.text.y = element_text(size = 9,
                               colour = "black"),
    
    axis.text.x = element_text(size = 6,
                               colour = "black"),
    
    axis.title.y = element_text(
      size = 11,
      face = "bold",
      colour = "black",
      family = "sans"
    )
  ) +
  geom_density(alpha = 0.2) 



