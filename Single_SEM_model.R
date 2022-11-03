## ---------------------------
##
## Script Purpose: Single SEM with individual proteins, first CRP
##
##
#load packages you may need
#library(lavaan)

#gFA_gMD <- read.csv("STRADL_Measures_DTI_Global.csv")
#names(gFA_gMD)[names(gFA_gMD) == 'ID'] <- 'stradl_ID'
#test2 <- merge(gFA_gMD, test, by = "stradl_ID")
#write.csv(test2, "n640_full_neuroimaging_data_inc_gFA_gMD.csv")

test2 <- read.csv("n640_full_neuroimaging_data_inc_gFA_gMD.csv")

test2 <- test2 %>% drop_na(g)  

# Converting multiple varibles into a factor
test2 %<>% mutate_at(c("sex", "site", "edited", "batch"),
                     as.factor)

skimr::skim(test2$sex)
#n = 640 for mediation 

model1_TBV <-
  '
            # direct effect (a path)
            COGNITION ~ cprime*CRP + b*TB + AGE + SEX
            # mediator (b path)
            TB  ~ a*CRP + AGE + SEX + SITE + EDITS + BATCH + ICV
            #(c prime path)
            #COGNITION ~ b*TB
            # indirect effect (a*b)
            TBIDE := a*b
            # total effect
            c := cprime + (a*b)

            '
tmp = data.frame(
  COGNITION = scale(test2$g),
  CRP = scale(test2$CRP),
  TB = scale(test2$global.wbv),
  SEX = test2$sex,
  AGE = scale(test2$st_age),
  ICV = scale(test2$Standardised_ICV),
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

TB <- table2 %>% filter(label == "TBIDE")

#### GM

model1_GM <-
  '
            # direct effect (a path)
            COGNITION ~ cprime*CRP + b*GM + AGE + SEX
            # mediator (b path)
            GM  ~ a*CRP + AGE + SEX + SITE + EDITS + BATCH + ICV
            #(c prime path)
            #COGNITION ~ b*GM
            # indirect effect (a*b)
            GMIDE := a*b
            # total effect
            c := cprime + (a*b)

            '
tmp = data.frame(
  COGNITION = scale(test2$g),
  CRP = scale(test2$CRP),
  GM = scale(test2$global.total.gm),
  SEX = test2$sex,
  AGE = scale(test2$st_age),
  ICV = scale(test2$Standardised_ICV),
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

GM <- table2 %>% filter(label == "GMIDE")

##model1_WM <-
'
            # direct effect (a path)
            COGNITION ~ cprime*CRP + b*WM + AGE + SEX
            # mediator (b path)
            WM  ~ a*CRP + AGE + SEX + SITE + EDITS + BATCH + ICV
            #(c prime path)
            #COGNITION ~ b*WM
            # indirect effect (a*b)
            WMIDE := a*b
            # total effect
            c := cprime + (a*b)

            '
tmp = data.frame(
  COGNITION = scale(test2$g),
  CRP = scale(test2$CRP),
  WM = scale(test2$global.cerebral.wm),
  SEX = test2$sex,
  AGE = scale(test2$st_age),
  ICV = scale(test2$Standardised_ICV),
  BATCH = test2$batch,
  SITE = test2$site,
  EDITS = test2$edited
)

Single_SEM <- sem(model1_WM, data = tmp, missing = "fiml")

summary(
  Single_SEM,
  rsq = T,
  standardized = TRUE,
  fit.measures = TRUE
)

##############outputting MEDIATION results into a  table#################
table <- parameterEstimates(Single_SEM, standardized = TRUE)

table2 <- table %>% filter(label == "WMIDE" | label == "c" | label == "cprime")

table2 %<>% select(label, est, se, pvalue, ci.lower, ci.upper)

WM <- table2 %>% filter(label == "WMIDE")

###
model1_gFA <-
  '
            # direct effect (a path)
            COGNITION ~ cprime*CRP + b*gFA + AGE + SEX
            # mediator (b path)
            gFA  ~ a*CRP + AGE + SEX + SITE 
            #(c prime path)
            #COGNITION ~ b*gFA
            # indirect effect (a*b)
            gFAIDE := a*b
            # total effect
            c := cprime + (a*b)

            '
tmp = data.frame(
  COGNITION = scale(test2$g),
  CRP = scale(test2$CRP),
  gFA = scale(test2$gFA),
  SEX = test2$sex,
  AGE = scale(test2$st_age),
  SITE = test2$site

)

Single_SEM <- sem(model1_gFA, data = tmp, missing = "fiml")

summary(
  Single_SEM,
  rsq = T,
  standardized = TRUE,
  fit.measures = TRUE
)

##############outputting MEDIATION results into a  table#################
table <- parameterEstimates(Single_SEM, standardized = TRUE)

table2 <- table %>% filter(label == "gFAIDE" | label == "c" | label == "cprime")

table2 %<>% select(label, est, se, pvalue, ci.lower, ci.upper)

gFA <- table2 %>% filter(label == "gFAIDE")
####

#### gMD

model1_gMD <-
  '
            # direct effect (a path)
            COGNITION ~ cprime*CRP + b*gMD + AGE + SEX
            # mediator (b path)
            gMD  ~ a*CRP + AGE + SEX + SITE 
            #(c prime path)
            #COGNITION ~ b*gMD
            # indirect effect (a*b)
            gMDIDE := a*b
            # total effect
            c := cprime + (a*b)

            '
tmp = data.frame(
  COGNITION = scale(test2$g),
  CRP = scale(test2$CRP),
  gMD = scale(test2$gMD),
  SEX = test2$sex,
  AGE = scale(test2$st_age),
  SITE = test2$site

)

Single_SEM <- sem(model1_gMD, data = tmp, missing = "fiml")

summary(
  Single_SEM,
  rsq = T,
  standardized = TRUE,
  fit.measures = TRUE
)

##############outputting MEDIATION results into a  table#################
table <- parameterEstimates(Single_SEM, standardized = TRUE)

table2 <- table %>% filter(label == "gMDIDE" | label == "c" | label == "cprime")

table2 %<>% select(label, est, se, pvalue, ci.lower, ci.upper)

gMD <- table2 %>% filter(label == "gMDIDE")

###

Single_SEM_neuroimaging_ab <- rbind(TB, GM, WM, gFA, gMD)

Single_SEM_neuroimaging_ab %<>% mutate(DNAm =
                                         "CRP")

x <- ggplot(Single_SEM_neuroimaging_ab,
            aes(y=est,
                x=reorder(label, -est),
                color = label,
                group = DNAm
                )) +
  
  geom_point(position = position_dodge(width = 0.9),
             size = 2.0,
             stroke = 0.9) +
  
  geom_errorbar(
    aes(
      ymin = ci.lower,
      ymax = ci.upper
    ),
    position = position_dodge(0.9),
    width = 0.2,
    colour = "darkgrey",
    alpha = 0.9,
    size = 0.8
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(axis.text.x = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(size = 12)) +
  
  coord_flip() +
  xlab("") +
  ylab("") +
  scale_y_continuous(limits = c(-0.05, 0.05))

#Visualise

x + 
  facet_wrap(vars(DNAm),  scales = "free_y", ncol =1) +
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
    )) +
  scale_shape_manual(values = c(16,
                                1)) +
  scale_colour_manual(values = c("darkgrey",
                                 "cornflowerblue",
                                  "#404788FF",
                                 "#CC0A4E",
                                 "#FF7F94"))


