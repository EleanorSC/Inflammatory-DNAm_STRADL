## ---------------------------
##
## Script Purpose: Code for correlation matrix of DNAm signatures 
##                 
##                 
##                 
##                 
##
##                
##
##  
## ---------------------------
## ---------------------------

methylation_proxies_KORA_LBC <-
  subset(
    KORA_LBC_DNAm,
    select = c(
      "ACY1",
      "ADAMTS13",
      "ADIPOQ",
      "AFM",
      "B2M",
      "BCAM",
      "BMP1",
      "C4A|C4B",
      "C5",
      "C9",
      "CCL17",
      "CCL18",
      "CCL21",
      "CCL22",
      "CCL25",
      "CD163",
      "CD209",
      "CD48",
      "CD5L",
      "CHIT1",
      "CLEC11A",
      "CLEC11A.1",
      "CNTN4",
      "CRP",
      "CXCL10.y",
      "CXCL11.y",
      "EDA",
      "ENPP7",
      "ESM1",
      "F7",
      "FAP",
      "FCER2",
      "FCGR3B",
      "GHR",
      "GNLY",
      "GP1BA",
      "GZMA.y",
      "HGFAC",
      "ICAM5",
      "IDUA",
      "IGFBP1",
      "IGFBP4",
      "IL19",
      "INSR",
      "LGALS3BP",
      "LGALS4",
      "LTA|LTB",
      "LTF",
      "LY9",
      "LYZ",
      "MIA",
      "MMP1",
      "MMP12",
      "MMP2",
      "MMP9",
      "MPL",
      "MPO",
      "MRC2",
      "MST1",
      "NCAM1",
      "NOTCH1",
      "NTRK3.y",
      "OMD",
      "PAPPA",
      "PIGR",
      "PRSS2",
      "RARRES2",
      "RETN",
      "S100A9",
      "SELE",
      "SELL",
      "SEMA3E",
      "SERPINA3",
      "SERPIND1",
      "SHBG",
      "SLITRK5",
      "SPOCK2",
      "STC1",
      "THBS2",
      "TNFRSF17",
      "TNFRSF1B",
      "TPSB2",
      "VCAM1",
      "WFIKKN2",
      
      # LBC trained proxies
      ##### WHICH STRADL KORA AND STRADL LBC match?
      # GZMA r = 0.71, MMP.1 r = 0.46, CXCL10 r = 0.35, NTRK3 r = 0.26, and CXCL11 r = 0.09
      # x = LBC, y = KORA
      
      "CCL11",
      "CD6",
      "CRTAM",
      "CXCL10.x",
      "CXCL11.x",
      "CXCL9",
      "EN.RAGE",
      "EZR",
      "FcRL2",
      "FGF.21",
      "G.CSF",
      "GDF.8",
      "GZMA.x",
      "HGF",
      "MMP.1",
      "N.CDase",
      "NEP",
      "NMNAT1",
      "NTRK3.x",
      "OSM",
      "SIGLEC1",
      "SKR3",
      "SMPD1",
      "TGF.alpha",
      "VEGFA"
    )
  )


#colnames(methylation_proxies_KORA_LBC) <- paste0(colnames(methylation_proxies_KORA_LBC),'_DNAm')


M <-cor(methylation_proxies_KORA_LBC, use = "pairwise.complete.obs")

install.packages("corrplot")
library(corrplot)
#library(wesanderson)
#library(RColorBrewer)
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(M)
head(p.mat[, 1:5])


corrplot(M, method = "color",
         #type = "upper",
         order = "hclust",
         p.mat = p.mat, 
         sig.level = 0.05, 
         insig = "blank",
         cl.pos = "b", 
         tl.col = "black", 
         tl.cex = 0.4, 
         cl.cex = 0.4, 
         addCoef.col = "black", 
         number.digits = 2, 
         number.cex = 0.01, 
         col = colorRampPalette(c("turquoise4","bisque", "#DC143C"))(100))

#col = wesanderson::wes_palette("Zissou1", 84, type = "continuous"))
#col = colorRampPalette(c("darkred","white","darkcyan"))(100)