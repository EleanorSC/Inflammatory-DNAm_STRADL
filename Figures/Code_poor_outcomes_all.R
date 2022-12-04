##### Figure 3 - which DNAm signatures associate with + brain and cognitive outcomes
# ----------------------------#
# CREATE SUPPLEMENTARY TABLE
# ----------------------------#
MMP12_test <- plot_global_methylation %>% filter(DNAm == "MMP12")


Poor_brain_health <- plot_global_methylation %>% 
  select(brain_metric, DNAm, estimate, std.error, r2, p.value, pFDR, modality)

Poor_brain_health_lifestyle <- plot_global_methylation_lifestyle %>%
  select(brain_metric, DNAm, estimate, std.error, r2, p.value, pFDR, modality) %>%
  rename(
    brain_metric_lifestyle = brain_metric,
    DNAm_lifestyle = DNAm,
    estimate_lifestyle = estimate,
    r2_lifestyle = r2,
    p.value_lifestyle = p.value,
    pFDR_lifestyle = pFDR,
    std.error_lifestyle = std.error,
    modality_lifestyle = modality
  )



###
## ----------------------------#
# Combining both sets of models to see percentage attenuation
## ----------------------------#

table_new <- cbind(Poor_brain_health, Poor_brain_health_lifestyle)


table_new %<>% mutate(estimate2 = abs(estimate),
                      estimate_lifestyle2 = abs(estimate_lifestyle),
                      difference =
                        case_when(
                          (estimate2 > estimate_lifestyle2) ~ (estimate2 - estimate_lifestyle2),
                          (estimate2 < estimate_lifestyle2) ~ (estimate_lifestyle2 - estimate2),
                          TRUE ~ 0),
                      percentage_increase_decrease = ((difference/estimate2)*100),
                      increase_or_decrease = 
                        case_when(
                          (estimate2 > estimate_lifestyle2) ~ "Attenuation",
                          (estimate2 < estimate_lifestyle2) ~ "Increase",
                          TRUE ~ "misc" )
)%>% 
  
  filter(brain_metric == "global grey matter" & estimate < 0 |
           brain_metric == "global white matter" & estimate < 0 |
           brain_metric == "total brain volume" & estimate < 0 |
           brain_metric == "global cortical thickness" & estimate < 0 |
           brain_metric == "global cortical volume" & estimate < 0 |
           brain_metric == "global cortical surface area" & estimate < 0 |
           brain_metric == "gFA" & estimate < 0 |
           brain_metric == "global subcortical volume" & estimate < 0 |
           brain_metric == "WMH" & estimate > 0 |
           brain_metric == "gMD" & estimate > 0 |
           brain_metric == "g" & estimate < 0 |
           brain_metric == "gf" & estimate < 0 |
           brain_metric == "relative brain age" & estimate > 0 |
           brain_metric == "APOE" & estimate > 0 |
           brain_metric =="processing speed" & estimate < 0 |
           brain_metric =="executive function" & estimate < 0 |
           brain_metric =="vocabulary" & estimate < 0|
           brain_metric =="verbal declarative memory" & estimate < 0|
           brain_metric =="matrix reasoning" & estimate < 0) %>% 
  
  
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
  
  filter(pFDR < 0.05) %>%
  
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
    increase_or_decrease,
    percentage_increase_decrease
  )    


table_new <- subset(table_new, select = -c(modality)) 

table_new %<>% group_by(brain_metric, DNAm) %>%
  arrange(desc(abs(estimate)),
          by_group = TRUE)
# ----------------------------#
# write to .csv
# ----------------------------#
write.csv(table_new, "poor_methylation_models.csv")

# bits for writeup results
range(abs(table_new$estimate))
unique(table_new$DNAm)

table2 <- table_new %>% filter(pFDR_lifestyle <0.05) %>% arrange(estimate_lifestyle)
unique(table2$DNAm)
range(abs(table2$estimate_lifestyle))

MMP12_poor <- table_new %>% filter(DNAm == "MMP12" & pFDR <0.05)

table2 <- table_new %>% filter(increase_or_decrease == "Attenuation")

table_new_att <- table_new %>% filter(increase_or_decrease == "Attenuation")
mean(table_new_att$percentage_increase_decrease)
range(abs(table_new$estimate_lifestyle))

table_new_att <- table_new %>% filter(increase_or_decrease == "Increase")
mean(table_new_att$percentage_increase_decrease)

table_new_att <- table_new %>% filter(pFDR_lifestyle <0.05 & increase_or_decrease == "Increase" & brain_metric == "relative brain age")




# ----------------------------#
# / END CREATE SUPPLEMENTARY TABLE
# ----------------------------#


Poor_brain_health2<- plot_global_methylation %>% 
  filter(brain_metric == "global grey matter" & estimate < 0 |
         brain_metric == "global white matter" & estimate < 0 |
         brain_metric == "total brain volume" & estimate < 0 |
         brain_metric == "global cortical thickness" & estimate < 0 |
         brain_metric == "global cortical volume" & estimate < 0 |
         brain_metric == "global cortical surface area" & estimate < 0 |
         brain_metric == "gFA" & estimate < 0 |
         brain_metric == "global subcortical volume" & estimate < 0 |
         brain_metric == "WMH" & estimate > 0 |
         brain_metric == "gMD" & estimate > 0 |
         brain_metric == "g" & estimate < 0 |
         brain_metric == "gf" & estimate < 0 |
         brain_metric == "relative brain age" & estimate > 0 |
         brain_metric == "APOE" & estimate > 0 |
         brain_metric =="processing speed" & estimate < 0 |
         brain_metric =="executive function" & estimate < 0 |
         brain_metric =="vocabulary" & estimate < 0|
         brain_metric =="verbal declarative memory" & estimate < 0|
         brain_metric =="matrix reasoning" & estimate < 0) 

### Which DNAm proxy has the most FDR significant hits?
Poor_brain_health2 %<>%
  mutate(number_significant =
           case_when(significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Poor_brain_health <-
  aggregate(
    Poor_brain_health2$number_significant,
    by = list(DNAm = Poor_brain_health2$DNAm),
    FUN = sum
  )

Top_hits_Poor_brain_health %<>% arrange(desc(x)) %>% rename(p = x)

# Run FDR
Poor_brain_health2 %<>%
  mutate(number_significant_FDR =
           case_when(FDR_significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Poor_brain_health_2 <-
  aggregate(
    Poor_brain_health2$number_significant_FDR,
    by = list(DNAm = Poor_brain_health2$DNAm),
    FUN = sum
  ) %>% 
  rename(pFDR = x)


Poor_brain <- merge(Top_hits_Poor_brain_health,
                    Top_hits_Poor_brain_health_2,
                    by = "DNAm") %>%
  mutate(direction =
           "Poor brain health outcome")


####
Poor_brain_sigs2 <- Poor_brain %>% filter(pFDR > 1) %>% arrange(desc(pFDR))

##
#CRP         MMP12       PIGR        IGFBP4      FGF.21      SKR3        PRSS2      
#SERPIND1    THBS2       VEGFA       ICAM5       MMP.1_olink RARRES2     TGF.alpha  
#ACY1        AFM         CCL17       CCL18       HGF         MMP1        MMP9       
#SIGLEC1  
#
Poor_brain_sigs <- Poor_brain %>% filter(p > 0) %>%
  mutate(DNAm = stringr::str_replace(DNAm, "_olink", ".x")) 

ggplot(Poor_brain_sigs,
       aes(
         y = reorder(DNAm, p),
         x = p
       )) +
  
  ##
  geom_col(aes(
    y = reorder(DNAm, p),
    x = p,
    fill = reorder(DNAm, p)
  )) +
  
  geom_col(aes(
    y = reorder(DNAm, p),
    x = pFDR,
    colour = c("grey"),
    alpha = 0.95
  )) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",
    
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
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(
      size = 7,
      face = "bold",
      family = "sans"
    ),
    axis.title.x = element_text(
      size = 7,
      face = "bold",
      family = "sans"
    )
  ) +
  
  labs(x = "",
       y = "") +
  
  facet_wrap( ~ direction) +
  scale_fill_manual(values = c(
    colorspace::sequential_hcl(102, 
                               palette = "SunsetDark") 
  )) +
  scale_x_reverse ()+
  scale_colour_manual(values = c("#808080"))
