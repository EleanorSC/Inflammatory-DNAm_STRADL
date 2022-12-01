##### Figure 3 - which DNAm signatures associate with + brain and cognitive outcomes
# ----------------------------#
# CREATE SUPPLEMENTARY TABLE
# ----------------------------#
Good_brain_health <- plot_global_methylation %>% 
  select(brain_metric, DNAm, estimate, std.error, r2, p.value, pFDR, modality)

Good_brain_health_lifestyle <- plot_global_methylation_lifestyle %>%
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

table_new <- cbind(Good_brain_health, Good_brain_health_lifestyle)


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

  filter(brain_metric == "global grey matter" & estimate > 0 |
           brain_metric == "global white matter" & estimate > 0 |
           brain_metric == "total brain volume" & estimate > 0 |
           brain_metric == "global cortical thickness" & estimate > 0 |
           brain_metric == "global cortical volume" & estimate > 0 |
           brain_metric == "global cortical surface area" & estimate > 0 |
           brain_metric == "gFA" & estimate > 0 |
           brain_metric == "global subcortical volume" & estimate > 0 |
           brain_metric == "WMH" & estimate < 0 |
           brain_metric == "gMD" & estimate < 0 |
           brain_metric == "g" & estimate > 0 |
           brain_metric == "gf" & estimate > 0 |
           brain_metric == "relative brain age" & estimate < 0 |
           brain_metric == "APOE" & estimate < 0 |
           brain_metric =="processing speed" & estimate > 0 |
           brain_metric =="executive function" & estimate > 0 |
           brain_metric =="vocabulary" & estimate > 0|
           brain_metric =="verbal declarative memory" & estimate > 0|
           brain_metric =="matrix reasoning" & estimate > 0) %>% 
  
  group_by(DNAm, brain_metric) %>%
  
  arrange(desc(estimate),
          .by_group = TRUE) %>%
  
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
# ----------------------------#
# write to .csv
# ----------------------------#
write.csv(table_new, "favourable_methylation_models.csv")

table_new_att <- table_new %>% filter(increase_or_decrease == "Increase")

mean(table_new_att$percentage_increase_decrease)


# ----------------------------#
# / END CREATE SUPPLEMENTARY TABLE
# ----------------------------#

Good_brain_health <- plot_global_methylation_lifestyle %>% 
  filter(brain_metric == "global grey matter" & estimate > 0 |
           brain_metric == "global white matter" & estimate > 0 |
           brain_metric == "total brain volume" & estimate > 0 |
           brain_metric == "global cortical thickness" & estimate > 0 |
           brain_metric == "global cortical volume" & estimate > 0 |
           brain_metric == "global cortical surface area" & estimate > 0 |
           brain_metric == "gFA" & estimate > 0 |
           brain_metric == "global subcortical volume" & estimate > 0 |
           brain_metric == "WMH" & estimate < 0 |
           brain_metric == "gMD" & estimate < 0 |
           brain_metric == "g" & estimate > 0 |
           brain_metric == "gf" & estimate > 0 |
           brain_metric == "relative brain age" & estimate < 0 |
           brain_metric == "APOE" & estimate < 0 |
           brain_metric =="processing speed" & estimate > 0 |
           brain_metric =="executive function" & estimate > 0 |
           brain_metric =="vocabulary" & estimate > 0|
           brain_metric =="verbal declarative memory" & estimate > 0|
           brain_metric =="matrix reasoning" & estimate > 0) %>% 
 # filter(pFDR < 0.05) %>%
  group_by(DNAm, metric)%>%
  arrange(estimate, by_group = TRUE)


### Which DNAm proxy has the most FDR significant hits?
Good_brain_health %<>%
  mutate(number_significant =
           case_when(significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Good_brain_health <-
  aggregate(
    Good_brain_health$number_significant,
    by = list(DNAm = Good_brain_health$DNAm),
    FUN = sum
  )

Top_hits_Good_brain_health %<>% arrange(desc(x)) %>% rename(p = x)

# Run FDR
Good_brain_health %<>%
  mutate(number_significant_FDR =
           case_when(FDR_significance == "No" ~ 0,
                     TRUE ~ 1))

# Group by sum using dplyr
Top_hits_Good_brain_health_2 <-
  aggregate(
    Good_brain_health$number_significant_FDR,
    by = list(DNAm = Good_brain_health$DNAm),
    FUN = sum
  ) %>% 
  rename(pFDR = x)


Good_brain <- merge(Top_hits_Good_brain_health,
                    Top_hits_Good_brain_health_2,
                    by = "DNAm") %>%
  mutate(direction =
           "Favourable brain health outcome")


####
Good_brain_sigs2 <- Good_brain %>% filter(pFDR > 1) 

Good_brain_sigs <- Good_brain %>% filter(p > 0) %>%
  mutate(DNAm = stringr::str_replace(DNAm, "_olink", ".x")) 

ggplot(Good_brain_sigs,
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
    colorspace::sequential_hcl(46, 
                               palette = "Dark Mint") 
  )) +
  scale_x_reverse ()+
  scale_colour_manual(values = c("#808080")) #+
  #scale_x_continuous(limits = c(0, 13))