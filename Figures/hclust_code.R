### Correlation for DNAm signatures, n =778

PROTEOMICS_DNAm_DATA_CORR <-
  subset(
    PROTEOMICS_DNAm_DATA,
    select = c(
      "CRTAM"           ,
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
      "LTA.LTB"      ,
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
      "C4A.C4B"        ,
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
      "SPOCK2"
    )
  )


### clusters
# Dissimilarity matrix
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend)


df <- scale(na.omit(PROTEOMICS_DNAm_DATA_CORR))


df <- M
d <- dist(M, method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )

# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1)

# Compute with agnes
hc2 <- agnes(df, method = "complete")

# Agglomerative coefficient
hc2$ac

#factoextra::fviz_nbclust(M, FUN = hcut, method = "wss")

# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# function to compute coefficient
ac <- function(x) {
  agnes(df, method = x)$ac
}

map_dbl(m, ac)

hc3 <- agnes(df, method = "ward")
pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 

# Ward's method
hc5 <- hclust(d, method = "ward.D2" )

# Cut tree into 9 groups
sub_grp <- cutree(hc5, k = 9)

# Number of members in each cluster
table(sub_grp)

plot(hc5, cex = 0.6)
rect.hclust(hc5, k = 9, border = 2:5)


fviz_cluster(list(data = df, cluster = sub_grp))
############
# CORRPLOT
############
M <-cor(na.omit(PROTEOMICS_DNAm_DATA_CORR), use = "pairwise.complete.obs")

corrplot(M, 
         order = "hclust",
         method = "color",
         number.digits = 2, 
         # addCoef.col="black", 
         number.cex = 0.5,
         tl.col = "black",
         tl.cex = 0.4, 
         cl.cex = 0.4,
         addrect = 10,
         col = colorRampPalette(c("#739dff",
                                  "#faefde",
                                  "#f55156"))(100)
)