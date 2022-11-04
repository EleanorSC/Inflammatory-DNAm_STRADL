## ---------------------------
##
## Script Purpose: Data loading and data cleaning
##   
##  General notes for neuroimaging data
----------------------------
  ##  
  ##          1. Estimated ICVs were standardised separately per site (Aberdeen / Dundee) in
  ##         file STRADL_Measures_Standardised_ICV.csv;
  ##
  ##         2. Quality control-related covariates in the STRADL_Measures_FreeSurfer_Main.csv file
  ##          (site / edited / batch) are recommended to be included in all analyses with the
  ##          FreeSurfer measures;
  ##
  ##         3. Global white matter integrity measures in file STRADL_Measures_DTI_Global.csv
  ##         were derived from all tracts in STRADL_Measures_DTI_Tracts_Readme.txt except
  ##         for corpus callosum, corona radiata and internal capsule.
##
##
## Notes: you will need the following data: 
##        (a) a .csv file containing: proteomics data [protein_data_full <- read.csv("prot_file_150621.csv")]
##        (b) a .csv file containing: LBC_trained DNAm data [EpiScores <- read.csv("STRADL_scores_trained_LBC.csv")]
##        (a) a .csv file containing: KORA_trained DNAm data [STRADL_KORA <- read.csv("STRADL_scores_trained_KORA.csv")]
##        (a) a .csv file containing: STRADL main neuroimaging data [STRADL_FreeSurfer <- read.csv("STRADL_Measures_FreeSurfer_Main.csv")]
##        (a) a .csv file containing: STRADL ICV data [STRADL_ICV <- read.csv("STRADL_Measures_Standardised_ICV.csv")
##        (a) a .csv file containing: STRADL dMRI data [STRADL_WM_tracts_FA <- read.csv("STRADL_Measures_DTI_FA.csv")]
##        (a) a .csv file containing: STRADL global dMRI data [gFA_gMD <- read.csv("STRADL_Measures_DTI_Global.csv")]
##        (a) a .csv file containing: STRADL main phenotype data (demographics etc)
##
## 
## Notes: the aim is to have these dataframes to work on:
##        (a) Standard_covariates   (full DNAm, neuroimaging and covariate data)
##        (a) Disease_covariates    (full DNAm, neuroimaging and covariate data + disease status) 
##
## ---------------------------


# First load up demographics and phenotypes
phenotypes_STRADL <- read.csv("phenotypes.csv")
disease_covariates <- read.csv("Depression_Inflammation_DNAm_final.csv")
ID_link_attempt <- read.csv("STRADL_DNAm_target_REM_17April2020.csv")

# Next load up methylation and proteomics data
protein_data_full <- read.csv("prot_file_150621.csv", check.names=FALSE)
annotation_data <- read.csv("Annotation_data.csv", check.names=FALSE)
# Methylation:
STRADL_KORA <- read.csv("STRADL_scores_trained_KORA.csv", check.names=FALSE)
EpiScores <- read.csv("STRADL_scores_trained_LBC.csv")
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
STRADL_main_data <- protein_data_full[c(1, 4237:4286)]

## disease_covariates <- read.csv("Depression_Inflammation_DNAm_final.csv")
## Add in some phenotypes for morbidity status
names(disease_covariates)[names(disease_covariates) == 'ID'] <- 'GS_id'

# Load ID linkages
GS_link <- ID_link_attempt %>% select("stradl_ID", "GS_id")

disease_covariates <- merge(GS_link, 
                            disease_covariates,
                            by = "GS_id")

disease_covariates %<>% select(stradl_ID, 
                               diabetes_Y, osteo_arthritis_Y, rheum_arthritis_Y, depression_Y,
                               asthma_Y, COPD_Y, inflam_disease,
                               CD8T,CD4T,NK, Bcell, Mono, Gran)

### Clean up this disease data
disease_covariates %<>% 
  mutate(inflam_disease = case_when(
    inflam_disease == 1 ~ 1,
    inflam_disease == 2 ~ 1,
    inflam_disease == 3 ~ 1,
    TRUE ~ 0))

disease_covariates$inflam_disease <- as.factor(disease_covariates$inflam_disease)
skimr::skim(disease_covariates$inflam_disease)

### Link what we need from STRADL main data into phenotypes + disease_covariates
STRADL_covariates <- STRADL_main_data %>% select(stradl_ID, sex, st_age,
                                                 APOE, apoe, Brain_age, brain_accel, CurrentSmoker,
                                                 Site, Global_GM_Volume, Cerebrum_WM_Volume,
                                                 g, gf, ethnicity)


names(phenotypes_STRADL)[names(phenotypes_STRADL) == 'id'] <- 'GS_id'
phenotypes_STRADL <- merge(GS_link, phenotypes_STRADL, by = "GS_id")

STRADL_alcohol <- phenotypes_STRADL %>% select(stradl_ID, units, drink_status, 
                                               hypertension_category, bmi, whr, body_fat)

STRADL_lifestyle_covariates <- merge(STRADL_alcohol, STRADL_covariates, by = "stradl_ID")

# CLEAN UP ALCOHOL DATA: use a case_when clause for current and not
# alc <- STRADL_alcohol %>% filter(drink_status == 1)

STRADL_lifestyle_covariates %<>% 
  mutate(CurrentDrinker = case_when(
    drink_status == 1 ~ 1,
    TRUE ~ 0))

STRADL_lifestyle_covariates$CurrentDrinker <- as.factor(STRADL_lifestyle_covariates$CurrentDrinker)
skimr::skim(STRADL_lifestyle_covariates$CurrentDrinker)


### DNA methylation data load
######Load up KORA DNAm proxies
#STRADL_KORA <- read.csv("STRADL_scores_trained_KORA.csv", check.names=FALSE)

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

# add DNAm to the end
#colnames(STRADL_KORA) <- paste(colnames(STRADL_KORA),"DNAm",sep="_")

# rename the ID columns 
names(STRADL_KORA)[names(STRADL_KORA) == 'meth_ID_DNAm'] <- 'meth_ID'
names(STRADL_KORA)[names(STRADL_KORA) == 'stradl_ID_DNAm'] <- 'stradl_ID'


##########
# Load DNAm EpiScores data trained in STRADL
#EpiScores <- read.csv("STRADL_scores_trained_LBC.csv")
names(EpiScores)[names(EpiScores) == 'ID'] <- 'meth_ID'

STRADL_LBC <- merge(ID_link_attempt, 
                    EpiScores,
                    by = "meth_ID")

# Remove columns we don't need
STRADL_LBC <- STRADL_LBC[c(2,7:31)]

# add DNAm to the end
#colnames(STRADL_LBC) <- paste(colnames(STRADL_LBC),"DNAm",sep="_")

# rename the ID columns 
names(STRADL_LBC)[names(STRADL_LBC) == 'stradl_ID_DNAm'] <- 'stradl_ID'

########
# Merge KORA trained and STRADL_trained
KORA_LBC_DNAm <- merge(STRADL_LBC, STRADL_KORA, by = "stradl_ID")


######Load up neuroimaging full dataset
STRADL_FreeSurfer <- read.csv("STRADL_Measures_FreeSurfer_Main.csv")
# Now link these to STRADL ID
names(STRADL_FreeSurfer)[names(STRADL_FreeSurfer) == 'id'] <- 'stradl_ID'


STRADL_ICV <- read.csv("STRADL_Measures_Standardised_ICV.csv")
# Now link these to STRADL ID
names(STRADL_ICV)[names(STRADL_ICV) == 'ID'] <- 'stradl_ID'


STRADL_MRI <- merge(STRADL_ICV, 
                    STRADL_FreeSurfer,
                    by = "stradl_ID")



# neuroimaging data including DTI data
STRADL_WM_tracts_FA <- read.csv("STRADL_Measures_DTI_FA.csv")
names(STRADL_WM_tracts_FA)[names(STRADL_WM_tracts_FA) == 'subjectID'] <- 'stradl_ID'
# drop column 'average FA'
STRADL_WM_tracts_FA <- STRADL_WM_tracts_FA[,-8]

gFA_gMD <- read.csv("STRADL_Measures_DTI_Global.csv")
names(gFA_gMD)[names(gFA_gMD) == 'ID'] <- 'stradl_ID'

STRADL_DTI <- merge(gFA_gMD, 
                    STRADL_WM_tracts_FA,
                    by = "stradl_ID")


#####

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

###########
STRADL_age_sex_only_covariates <- STRADL_main_data %>% select(stradl_ID, sex, st_age)

# n = 595
STRADL_lifestyle_disease_covariates <- merge(STRADL_lifestyle_covariates, disease_covariates, by = "stradl_ID")

disease_sensitivity <- merge(STRADL_lifestyle_disease_covariates, dataset_n709, by = "stradl_ID")

test <- merge(STRADL_age_sex_only_covariates, dataset_n709, by = "stradl_ID")

n709_noICV <- merge(STRADL_lifestyle_covariates, dataset_n709, by = "stradl_ID") 

ICV_corrected <- merge(STRADL_lifestyle_covariates, dataset_n655, by = "stradl_ID") 


### n = 778
COGNTIIVE_DATA <- merge(KORA_LBC_DNAm, STRADL_lifestyle_covariates, by = "stradl_ID")

STRADL_DTI <- merge(gFA_gMD, 
                    STRADL_WM_tracts_FA,
                    by = "stradl_ID")

test3 <- merge(ICV_corrected, STRADL_DTI , by = "stradl_ID") 

DTI_2 <-merge(STRADL_lifestyle_covariates,STRADL_DTI, by = "stradl_ID")
