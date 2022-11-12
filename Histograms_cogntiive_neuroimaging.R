## ---------------------------
##
## Script Purpose: Histogram plots of cognitive and brain metrics
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
#cognition_STRADL  <- na.omit(cognition_STRADL)

dataset_KORA <- merge(STRADL_KORA, 
                      STRADL_lifestyle_covariates,
                      by = "stradl_ID")


cognition_STRADL <- merge(dataset_KORA,
                          cognition_STRADL,
                          by = "stradl_ID")

hist(cognition_STRADL$matrix_reasoning)

#plot multiple histograms

df <- cognition_STRADL %>% select(stradl_ID,
                                  Digit_Symbol_coding, # Processing speed
                                  Verbal_Fluency, # Executive function
                                  Mill_Hall_vocab,
                                  matrix_reasoning,
                                  Logical_memory)

#Convert into long format

data_long <- reshape2::melt(df,
                  # ID variables - all the variables to keep but not split apart on
                  id.vars=c("stradl_ID"),
                  # The source columns
                  measure.vars=c("Digit_Symbol_coding", 
                                 "Verbal_Fluency", 
                                 "Mill_Hall_vocab",
                                 "matrix_reasoning",
                                 "Logical_memory"),
                  # Name of the destination column that will identify the original
                  # column that the measurement came from
                  variable.name="cognitive_test",
                  value.name="measurement")


# Calculate the mean of each
stats <- na.omit(data_long) %>% 
  group_by(cognitive_test) %>% 
  summarise(mean = mean(measurement),
            n = n()
            )


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
#
#ggplot(data_long, 
#       aes(measurement, fill = cognitive_test)
#       ) + geom_density(alpha = 0.2) +
#  
#  facet_wrap(~cognitive_test,
#             scales = "free_x",
#             nrow = 1) +
#  
#  viridis::scale_fill_viridis(discrete = TRUE,
#                              option = "D") +
#  
#  theme_classic() 

## ----------------------------# 
# PLOT 2 
## ----------------------------#

#ggplot(data_long, 
#       aes(measurement, 
#           fill = cognitive_test)
#       ) + 
#  
#  geom_bar(pos="dodge",
#           alpha = 0.8) +
#  
#  geom_vline(data = stats, 
#             aes(xintercept = mean, 
#                 color = cognitive_test), 
#             size = 2) +
#  
#  facet_wrap(~cognitive_test,
#             scales = "free_x",
#             nrow = 1) +
#  viridis::scale_fill_viridis(discrete = TRUE,
#                               option = "D") +
#  theme_classic() 

## ----------------------------# 
# PLOT 3 - HISTOGRAM TRADITIONAL 
## ----------------------------#

ggplot(data_long, 
       aes(measurement, fill = cognitive_test)
       ) + 
  
  geom_histogram(alpha = 0.5, 
                 aes(y = ..density..), 
               #  binwidth = 30,
                 position = 'identity') +
  
  geom_vline(data = stats, 
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
             nrow = 1) +
  
  theme(axis.text.y = element_text(size = 7)) +
  
  labs (x = "test score",
       y = "frequency") +
  
  viridis::scale_fill_viridis(discrete = TRUE,
                              option = "D") +
  
  viridis::scale_colour_viridis(discrete = TRUE,
                              option = "D") +
  theme_classic() +
  geom_density(alpha = 0.2)
