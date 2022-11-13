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
  CRTAM           ,
  EZR             ,
  FcRL2           ,
  G.CSF           ,
  GDF.8           ,
  GZMA_olink      ,
  N.CDase         ,
  NEP             ,
  NMNAT1          ,
  NTRK3_olink    ,
  SIGLEC1        ,
  SKR3           ,
  SMPD1          ,
  CCL11          ,
  CD6            ,
  CXCL10_olink   ,
  CXCL11_olink   ,
  CXCL9          ,
  EN.RAGE        ,
  FGF.21         ,
  HGF            ,
  MMP.1_olink    ,
  OSM            ,
  TGF.alpha      ,
  VEGFA          ,
  CCL21          ,
  MMP9           ,
  MPO            ,
  NTRK3          ,
  TNFRSF17       ,
  MIA            ,
  CCL25          ,
  IGFBP1         ,
  LTF            ,
  BCAM           ,
  EDA            ,
  C5             ,
  GHR            ,
  IGFBP4         ,
  CLEC11A        ,
  VCAM1          ,
  LGALS4         ,
  CD209          ,
  IL19           ,
  CXCL11         ,
  MRC2           ,
  CCL18          ,
  RETN           ,
  C9             ,
  RARRES2        ,
  TNFRSF1B       ,
  IDUA           ,
  ADAMTS13       ,
  F7             ,
  GNLY           ,
  PIGR           ,
  WFIKKN2        ,
  FCER2          ,
  CD48           ,
  CD5L           ,
  CNTN4          ,
  FCGR3B         ,
  SERPIND1       ,
  LY9            ,
  THBS2          ,
  ACY1           ,
  BMP1           ,
  TPSB2          ,
  GZMA           ,
  INSR           ,
  SELE           ,
  MPL            ,
  B2M            ,
  `LTA|LTB`      ,
  CCL22          ,
  CCL17          ,
  ADIPOQ         ,
  CHIT1          ,
  HGFAC          ,
  ESM1           ,
  CXCL10         ,
  PAPPA          ,
  SERPINA3       ,
  MMP2           ,
  CRP            ,
  MST1           ,
  ENPP7          ,
  `C4A|C4B`        ,
  MMP12          ,
  NCAM1          ,
  CLEC11A.1      ,
  SLITRK5        ,
  AFM            ,
  SELL           ,
  LYZ            ,
  MMP1           ,
  SHBG           ,
  STC1           ,
  GP1BA          ,
  LGALS3BP      ,
  CD163         ,
  FAP           ,
  PRSS2         ,
  NOTCH1        ,
  ICAM5         ,
  S100A9        ,
  OMD           ,
  SEMA3E        ,
  SPOCK2
)

#Convert into long format

data_long <- reshape2::melt(df,
                            # ID variables - all the variables to keep but not split apart on
                            id.vars=c("stradl_ID"),
                            # The source columns
                            measure.vars=c("CRTAM"           ,
                                           "EZR"             ,
                                           "FcRL2"           ,
                                           "G.CSF"           ,
                                           "GDF.8"           ,
                                           "GZMA_olink"      ,
                                           "N.CDase"         ,
                                           "NEP"             ,
                                           "NMNAT1"          ,
                                           "NTRK3_olink"    ,
                                           "SIGLEC1"        ,
                                           "SKR3"           ,
                                           "SMPD1"          ,
                                           "CCL11"          ,
                                           "CD6"            ,
                                           "CXCL10_olink"   ,
                                           "CXCL11_olink"   ,
                                           "CXCL9"          ,
                                           "EN.RAGE"        ,
                                           "FGF.21"         ,
                                           "HGF"            ,
                                           "MMP.1_olink"    ,
                                           "OSM"            ,
                                           "TGF.alpha"      ,
                                           "VEGFA"          ,
                                           "CCL21"          ,
                                           "MMP9"           ,
                                           "MPO"            ,
                                           "NTRK3"          ,
                                           "TNFRSF17"       ,
                                           "MIA"            ,
                                           "CCL25"          ,
                                           "IGFBP1"         ,
                                           "LTF"            ,
                                           "BCAM"           ,
                                           "EDA"            ,
                                           "C5"             ,
                                           "GHR"            ,
                                           "IGFBP4"         ,
                                           "CLEC11A"        ,
                                           "VCAM1"          ,
                                           "LGALS4"         ,
                                           "CD209"          ,
                                           "IL19"           ,
                                           "CXCL11"         ,
                                           "MRC2"           ,
                                           "CCL18"          ,
                                           "RETN"           ,
                                           "C9"             ,
                                           "RARRES2"        ,
                                           "TNFRSF1B"       ,
                                           "IDUA"           ,
                                           "ADAMTS13"       ,
                                           "F7"             ,
                                           "GNLY"           ,
                                           "PIGR"           ,
                                           "WFIKKN2"        ,
                                           "FCER2"          ,
                                           "CD48"           ,
                                           "CD5L"           ,
                                           "CNTN4"          ,
                                           "FCGR3B"         ,
                                           "SERPIND1"       ,
                                           "LY9"            ,
                                           "THBS2"          ,
                                           "ACY1"           ,
                                           "BMP1"           ,
                                           "TPSB2"          ,
                                           "GZMA"           ,
                                           "INSR"           ,
                                           "SELE"           ,
                                           "MPL"            ,
                                           "B2M"            ,
                                           "LTA|LTB"      ,
                                           "CCL22"          ,
                                           "CCL17"          ,
                                           "ADIPOQ"         ,
                                           "CHIT1"          ,
                                           "HGFAC"          ,
                                           "ESM1"           ,
                                           "CXCL10"         ,
                                           "PAPPA"          ,
                                           "SERPINA3"       ,
                                           "MMP2"           ,
                                           "CRP"            ,
                                           "MST1"           ,
                                           "ENPP7"          ,
                                           "C4A|C4B"        ,
                                           "MMP12"          ,
                                           "NCAM1"          ,
                                           "CLEC11A.1"      ,
                                           "SLITRK5"        ,
                                           "AFM"            ,
                                           "SELL"           ,
                                           "LYZ"            ,
                                           "MMP1"           ,
                                           "SHBG"           ,
                                           "STC1"           ,
                                           "GP1BA"          ,
                                           "LGALS3BP"      ,
                                           "CD163"         ,
                                           "FAP"           ,
                                           "PRSS2"         ,
                                           "NOTCH1"        ,
                                           "ICAM5"         ,
                                           "S100A9"        ,
                                           "OMD"           ,
                                           "SEMA3E"        ,
                                           "SPOCK2"),
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
