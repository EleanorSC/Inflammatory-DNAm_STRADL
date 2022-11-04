# PheWAS plot attempt


CRP <- plot2 %>% filter(DNAm == "CRP")

ggplot(CRP, aes(x = metric, 
                y = -log(p.value)
                )
       ) +
  
  geom_point(aes(col = brain_var,
                 size = -(estimate))
             ) +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_blank(),
    panel.grid.minor = element_line(colour = "grey", linetype = "dashed"),
    axis.ticks = element_blank()
  ) +
  
  labs(color = "Category",
       size = "Effect size",
       x = "neuroimaging metric",
       y = "log(p-value)") +
  
  ggrepel::geom_text_repel(
    data = . %>% mutate(label = ifelse(p.value < 0.05, as.character(metric), "")),
    aes(label = label),
    size = 4.1,
    box.padding = unit(0.7, "lines")
  ) +
  
  geom_hline(
    yintercept = -log(0.01),
    color = "red",
    size = 1,
    alpha = 0.5
  ) +
  
  geom_hline(
    yintercept = -log(0.05),
    color = "darkgrey",
    size = 1,
    alpha = 0.5
  ) 