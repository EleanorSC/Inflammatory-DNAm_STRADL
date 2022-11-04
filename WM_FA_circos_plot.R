## ---------------------------
##
## Script Purpose: Circos plot showing significant associations of WM FA with DNAm-signatures
##
## n.b "plot2" is a dataframe created in "iterate_through_two_lists_function.R"
##
## ---------------------------
plot3 <- plot2 %>% filter(Hemisphere == "Left")

plot3 <- plot3 %>% filter(significance == "Yes")

plot3 <- plot3 %>% select(DNAm, metric, estimate)
  
plot3 %<>% 
  dplyr::rename(
    rowname = DNAm,
    key = metric,
    value = estimate)

plot3 <- subset(plot3, select = -c(Hemisphere))

plot3$rowname <- as.character(plot3$rowname)

#install.packages("circlize")
#library(circlize)
# parameters
circos.clear()
circos.par(
  gap.after = c(rep(5, nrow(plot3)-1), 15, rep(5, ncol(plot3)-1), 15),
  start.degree = 90, gap.degree = 0.1, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE
           )
par(mar = rep(0, 4))


#install.packages("colorspace")
pal <- colorspace::sequential_hcl(88, palette = "SunsetDark")
col_fun = function(x) ifelse(x < -0.1 | x > 0.1, pal, # "#ECECEC"
                             "#00000000"
                             )


# Base plot
chordDiagram(
  x = plot3, 
  big.gap = 30,
#  order = c(    "SEMA3E",    "PIGR",      "SERPIND1",  "VEGFA",     "MMP12",     "MMP.1",     "NTRK3.x",   "SKR3",      "CCL11",    
#                "CRP",       "NCAM1",     "NTRK3.y",   "CCL17",     "CCL18",     "FGF.21",    "ICAM5",     "PRSS2",     "SLITRK5",  
#                "WFIKKN2",   "HGF",       "BCAM",      "IGFBP4",    "CNTN4",     "NOTCH1",    "OMD",       "TGF.alpha", "INSR",     
#                "SPOCK2",    "AFM",       "CXCL11.x",  "GDF.8",     "GZMA.y",    "SIGLEC1",   "ESM1",      "MMP9",      "EDA",      
#                "ACY1",      "C4A|C4B",   "CCL22",     "ENPP7",     "GZMA.x",    "IL19",      "NEP",       "SELL",      "LYZ",      
#                "MMP1",      "CD209",     "CD5L",      "CD6",       "SHBG",      "THBS2",     "CD48",      "EN.RAGE",   "F7",       
#                "LTA|LTB",   "MMP2",      "SELE",      "STC1",      "VCAM1",     "LTF",       "OSM",       "SERPINA3",  "ADAMTS13", 
#                "CRTAM",     "GNLY",      "LGALS3BP",  "LGALS4",    "MIA" ,      "NMNAT1",    "RARRES2",   "ADIPOQ",    "BMP1",     
#                "C9",        "CCL25",     "CLEC11A",   "CLEC11A.1", "CXCL9",     "FcRL2",     "G.CSF",     "MPO",       "TPSB2",   
#                "B2M",       "C5",        "CCL21",     "CD163",     "CHIT1",     "CXCL10.x",  "CXCL10.y",  "CXCL11.y",  "EZR",      
#                "FAP",       "FCER2",     "FCGR3B",    "GHR",       "GP1BA",     "HGFAC",     "IDUA",      "IGFBP1",    "LY9",      
#                "MPL",       "MRC2",      "MST1",      "N.CDase",   "PAPPA",     "RETN",      "S100A9",    "SMPD1",     "TNFRSF17", 
#                "TNFRSF1B", )
  grid.col = pal,
  col = col_fun,
  transparency = 0.25,
 symmetric = TRUE,
  #directional = 1,
  #direction.type = c("arrows", "diffHeight"), 
  #diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.decreasing = TRUE
  #link.largest.ontop = TRUE
)

# Add text and axis
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 2.5, 
      labels = sector.index, 
      #adj = c(0, degree(5)), 
      niceFacing = TRUE,
      facing = "clockwise",
      cex = 0.5
    )
  
    # Attempt to separate labels 
  #circos.axis(h = "top", 
  #            labels.cex = 0.5, 
  #            major.tick.percentage = 0.2, 
  #            sector.index = get.cell.meta.data("sector.index")
  #            
  #            #track.index = 2
  #            )
  #
  }
)


### Now try but per WM tract e.g CR


