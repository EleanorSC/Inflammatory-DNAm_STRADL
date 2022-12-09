### Protein histograms


## ---------------------------
##
## Script Purpose: Histogram plots of DNAm signatures n = 778
##
##                
##
##                
##
##  
## ----------------------------

df <- PROTEOMICS_DNAm_DATA %>% select(
  stradl_ID,
  ACY1_proteomics    ,
  ADAMTS13_proteomics  ,
  ADIPOQ_proteomics,
  AFM_proteomics       ,
  B2M_proteomics      ,
  BCAM_proteomics     ,
  BMP1_proteomics   ,
  C4A.C4B_proteomics ,
  C5_proteomics      ,
  C9_proteomics       ,
  CCL17_proteomics     ,
  CCL18_proteomics   ,
  CCL21_proteomics    ,
  CCL22_proteomics  ,
  CCL25_proteomics  ,
  CD163_proteomics ,
  CD209_proteomics   ,
  CD48_proteomics      ,
  CD5L_proteomics      ,
  CHIT1_proteomics     ,
  CLEC11A_proteomics  ,
  CLEC11A.1_proteomics,
  CNTN4_proteomics  ,
  CRP_proteomics     ,
  CXCL10_proteomics  ,
  CXCL11_proteomics    ,
  EDA_proteomics       ,
  ENPP7_proteomics     ,
  ESM1_proteomics     ,
  F7_proteomics       ,
  FAP_proteomics    ,
  FCER2_proteomics   ,
  FCGR3B_proteomics  ,
  GHR_proteomics       ,
  GZMA_proteomics ,
  GNLY_proteomics      ,
  GP1BA_proteomics     ,
  HGFAC_proteomics    ,
  ICAM5_proteomics  ,
  IDUA_proteomics    ,
  IGFBP1_proteomics  ,
  IGFBP4_proteomics    ,
  IL19_proteomics      ,
  INSR_proteomics      ,
  LGALS3BP_proteomics ,
  LGALS4_proteomics   ,
  LTA.LTB_proteomics,
  LTF_proteomics   ,
  LY9_proteomics     ,
  LYZ_proteomics     ,
  MIA_proteomics       ,
  MMP1_proteomics     ,
  MMP12_proteomics    ,
  MMP2_proteomics   ,
  MMP9_proteomics    ,
  MPL_proteomics     ,
  MPO_proteomics       ,
  MRC2_proteomics      ,
  MST1_proteomics      ,
  NCAM1_proteomics    ,
  NOTCH1_proteomics   ,
  NTRK3_proteomics  ,
  OMD_proteomics     ,
  PAPPA_proteomics   ,
  PIGR_proteomics      ,
  PRSS2_proteomics     ,
  RARRES2_proteomics   ,
  RETN_proteomics     ,
  S100A9_proteomics   ,
  SELE_proteomics   ,
  SELL_proteomics   ,
  SEMA3E_proteomics  ,
  SERPINA3_proteomics  ,
  SERPIND1_proteomics ,
  SHBG_proteomics      ,
  SLITRK5_proteomics  ,
  SPOCK2_proteomics   ,
  STC1_proteomics   ,
  THBS2_proteomics   ,
  TNFRSF17_proteomics  ,
  TNFRSF1B_proteomics ,
  TPSB2_proteomics     ,
  VCAM1_proteomics    ,
  WFIKKN2_proteomics,
  
  # LBC trained
  
  CRTAM_proteomics    ,
  EZR_proteomics      ,
  NMNAT1_proteomics   ,
  SMPD1_proteomics    ,
  CCL11_proteomics  ,
  CXCL9_proteomics    ,
  HGF_proteomics      ,
  OSM_proteomics      ,
  VEGFA_proteomics )


#Convert into long format

data_long <- reshape2::melt(df,
                            # ID variables - all the variables to keep but not split apart on
                            id.vars=c("stradl_ID"),
                            # The source columns
                            measure.vars=c("ACY1_proteomics"    ,
                                           "ADAMTS13_proteomics"  ,
                                           "ADIPOQ_proteomics",
                                           "AFM_proteomics"       ,
                                           "B2M_proteomics"      ,
                                           "BCAM_proteomics"     ,
                                           "BMP1_proteomics"   ,
                                           "C4A.C4B_proteomics" ,
                                           "C5_proteomics"      ,
                                           "C9_proteomics"       ,
                                           "CCL17_proteomics"     ,
                                           "CCL18_proteomics"   ,
                                           "CCL21_proteomics"    ,
                                           "CCL22_proteomics"  ,
                                           "CCL25_proteomics"  ,
                                           "CD163_proteomics" ,
                                           "CD209_proteomics"   ,
                                           "CD48_proteomics"      ,
                                           "CD5L_proteomics"      ,
                                           "CHIT1_proteomics"     ,
                                           "CLEC11A_proteomics"  ,
                                           "CLEC11A.1_proteomics",
                                           "CNTN4_proteomics"  ,
                                           "CRP_proteomics"     ,
                                           "CXCL10_proteomics"  ,
                                           "CXCL11_proteomics"    ,
                                           "EDA_proteomics"       ,
                                           "ENPP7_proteomics"     ,
                                           "ESM1_proteomics"     ,
                                           "F7_proteomics"       ,
                                           "FAP_proteomics"    ,
                                           "FCER2_proteomics"   ,
                                           "FCGR3B_proteomics"  ,
                                           "GHR_proteomics"       ,
                                           "GZMA_proteomics" ,
                                           "GNLY_proteomics"      ,
                                           "GP1BA_proteomics"     ,
                                           "HGFAC_proteomics"    ,
                                           "ICAM5_proteomics"  ,
                                           "IDUA_proteomics"    ,
                                           "IGFBP1_proteomics"  ,
                                           "IGFBP4_proteomics"    ,
                                           "IL19_proteomics"      ,
                                           "INSR_proteomics"      ,
                                           "LGALS3BP_proteomics" ,
                                           "LGALS4_proteomics"   ,
                                           "LTA.LTB_proteomics",
                                           "LTF_proteomics"   ,
                                           "LY9_proteomics"     ,
                                           "LYZ_proteomics"     ,
                                           "MIA_proteomics"       ,
                                           "MMP1_proteomics"     ,
                                           "MMP12_proteomics"    ,
                                           "MMP2_proteomics"   ,
                                           "MMP9_proteomics"    ,
                                           "MPL_proteomics"     ,
                                           "MPO_proteomics"       ,
                                           "MRC2_proteomics"      ,
                                           "MST1_proteomics"      ,
                                           "NCAM1_proteomics"    ,
                                           "NOTCH1_proteomics"   ,
                                           "NTRK3_proteomics"  ,
                                           "OMD_proteomics"     ,
                                           "PAPPA_proteomics"   ,
                                           "PIGR_proteomics"      ,
                                           "PRSS2_proteomics"     ,
                                           "RARRES2_proteomics"   ,
                                           "RETN_proteomics"     ,
                                           "S100A9_proteomics"   ,
                                           "SELE_proteomics"   ,
                                           "SELL_proteomics"   ,
                                           "SEMA3E_proteomics"  ,
                                           "SERPINA3_proteomics"  ,
                                           "SERPIND1_proteomics" ,
                                           "SHBG_proteomics"      ,
                                           "SLITRK5_proteomics"  ,
                                           "SPOCK2_proteomics"   ,
                                           "STC1_proteomics"   ,
                                           "THBS2_proteomics"   ,
                                           "TNFRSF17_proteomics"  ,
                                           "TNFRSF1B_proteomics" ,
                                           "TPSB2_proteomics"     ,
                                           "VCAM1_proteomics"    ,
                                           "WFIKKN2_proteomics",
                                           
                                           # LBC trained
                                           
                                           "CRTAM_proteomics"    ,
                                           "EZR_proteomics"      ,
                                           "NMNAT1_proteomics"   ,
                                           "SMPD1_proteomics"    ,
                                           "CCL11_proteomics"  ,
                                           "CXCL9_proteomics"    ,
                                           "HGF_proteomics"      ,
                                           "OSM_proteomics"      ,
                                           "VEGFA_proteomics" ),
                            # Name of the destination column that will identify the original
                            # column that the score came from
                            variable.name="DNAm",
                            value.name="score")


# Calculate the mean of each
## Mean and SD
msd <- na.omit(data_long) %>% 
  group_by(DNAm) %>% 
  summarise(mean= mean(score), 
            sd=sd(score),
            n = n()
  )

## ----------------------------# 
# PLOT 1
## ----------------------------#

ggplot(data_long, 
       aes(score, fill = DNAm)
) + 
  
  geom_histogram(alpha = 0.5, 
                 aes(y = ..density..), 
                 #  binwidth = 30,
                 position = 'identity') +
  
  geom_vline(data = msd, 
             aes(xintercept = mean, 
                 color = DNAm), 
             size = 1,
             alpha = 0.6
  ) +
  
  geom_vline(data = msd, 
             aes(xintercept = mean - sd, 
                 color = DNAm), 
             size = 0.5,
             alpha = 0.6,
             linetype="dotdash") +
  
  geom_vline(data = msd, 
             aes(xintercept = mean + sd, 
                 color = DNAm), 
             size = 0.5,
             alpha = 0.6,
             linetype="dotdash") +
  
  facet_wrap(~DNAm,
             scales = "free"
             #nrow = 2
  ) +
  
  labs (x = "test score",
        y = "frequency") +
  
  viridis::scale_fill_viridis(discrete = TRUE,
                              option = "F",
                              direction = -1) +
  
  viridis::scale_colour_viridis(discrete = TRUE,
                                option = "F",
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
    
    legend.position="none",
    
    
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