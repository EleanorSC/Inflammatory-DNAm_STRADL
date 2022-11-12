## ---------------------------
##
## Script Purpose: Derive a general factor of cognitive ability from the 5 following variables:
##
##                (1a)  mema, #Logical Memory I: Story A immediate recall
##                (1b)  memdela, #Logical Memory I: Story A delayed recall
##                (1a + 1b)  Logical_memory, #Logical Memory I
##                (2) digsym, #WAIS III - Digit Symbol-coding total score
##                (3) vftot, #Verbal Fluency total score
##                (4) mhv, #Mill Hill Vocabulary total score
##                (5) mrtotc #Matrix Reasoning total correct
##                 
##                 
##                 
##
##                In this script, we will derive g from dataset where DNAm data is available (n=778)
##
##  
## ----------------------------

#cog <- read.csv("cognitive.csv")
#skimr::skim(cog)

# cognition ####
library(readxl)
library(readr)
library(tidyr)
library(dplyr)


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
         Mill_Hall_vocab = mhv,
         matrix_reasoning = mrtotc
       #  matrix_reasoning = mrtottime
         ) %>%
  mutate(Logical_memory = mema_STRADL+memdela_STRADL)

#Examine data
skimr::skim(STRADL_cogn)

# prepare data
cognition_STRADL <- STRADL_cogn %>%
  dplyr::select(stradl_ID, 
                Digit_Symbol_coding, # Processing speed
                Verbal_Fluency, # Executive function
                Mill_Hall_vocab,
                matrix_reasoning,
                Logical_memory ## n.b also verbal declarative memory
                )

#Examine data
skimr::skim(cognition_STRADL)


# remove missing values
cognition_STRADL  <- na.omit(cognition_STRADL)

dataset_KORA <- merge(STRADL_KORA, 
                      STRADL_lifestyle_covariates,
                      by = "stradl_ID")

cognition_STRADL <- merge(dataset_KORA,
                          cognition_STRADL,
                                     by = "stradl_ID")

cognition_STRADL_DNAm <- cognition_STRADL %>% 
  select(
  stradl_ID,
  Digit_Symbol_coding,
  Verbal_Fluency,
  Mill_Hall_vocab,
  matrix_reasoning,
  Logical_memory
)

skimr::skim(cognition_STRADL_DNAm)

### scree plot attempt
#install.packages("factoextra")
#library(factoextra)

cognition_STRADL_PCA <- cognition_STRADL_DNAm %>%
  column_to_rownames('stradl_ID')


# Run PCA
res.pca <- prcomp(cognition_STRADL_PCA, scale = TRUE)

# Print scree plot
fviz_eig(res.pca)
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
#make eigen value plot for each metric
fviz_eig(res.pca, 
         geom=c("line"),
         main = "PCA scree plot cognition")

#extract the eigen values and % explained variance
eig_val <- get_eigenvalue(res.pca)
eig_val <- tibble::rownames_to_column(eig_val, "Component")

########################################
#extract the rotation/loadings
loadings <- as.data.frame(res.pca$rotation)
loadings <- tibble::rownames_to_column(loadings, "Cognitive test")


########################################
#extract the individual PCA results
pca_ind <- get_pca_ind(res.pca)

#take out the coordinates (the first dimension is the general factor of this metric for each subject)
ind_coord <- as.data.frame(pca_ind$coord)

rownames(ind_coord) <- cognition_STRADL_DNAm$stradl_ID
ind_coord <- tibble::rownames_to_column(ind_coord, "stradl_ID")


#######################################
#extract the variable PCA results
pca_var <- get_pca_var(res.pca)

#coordinates for each tract
var_coord <- as.data.frame(pca_var$coord)
var_coord <- tibble::rownames_to_column(var_coord, "Cognitive test")

#### check out loadings inversion

pca_component <- ind_coord %>% select(stradl_ID, Dim.1)

if (sum(var_coord$Dim.1) < 0){
  pca_component$pca_factor <- -1*(pca_component$Dim.1)
} else {
  pca_component$pca_factor <- (pca_component$Dim.1)
}

pca_component <- pca_component %>% select(-Dim.1)

# derive g (n=771)
pca_component %<>%
rename(g_STRADL = pca_component)



#######################################
#extract the variable PCA results
pca_var <- get_pca_var(pca)

#coordinates for each tract
var_coord <- as.data.frame(pca_var$coord)
var_coord <- tibble::rownames_to_column(var_coord, "Tract")

####################################
#print the PCA tract correlations with the 1st component
print.data.frame(var_coord[, 1:2])


# Extract the results for variables
var <- get_pca_var(res.pca)
var

# Coordinates of variables
head(var$coord)

# Contribution of variables
head(var$contrib)

# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)

#####
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

## Extract eigenvalues/variances
get_eig(res.pca)

### Match up to main dataset

#pca_component
#cognition_STRADL

dataset_KORA <- merge(STRADL_KORA, 
                      STRADL_lifestyle_covariates,
                      by = "stradl_ID")

KORA_COGNITION <- merge(dataset_KORA, 
                        cognition_STRADL,
                      by = "stradl_ID")

KORA_COGNITION <- merge(KORA_COGNITION, 
                        pca_component,
                        by = "stradl_ID")

