# This script makes a Venn diagram

library(RColorBrewer)
 myCol <- brewer.pal(3, "Pastel2")
 library(RColorBrewer)
 myCol <- brewer.pal(4, "Pastel2")

venn.diagram(x = list(y, s, x, z), category.names = c("LDWT" , "SDWT" , "SDCLOCK", "LDCLOCK"), filename = 'IntersectedGenes.png', output=TRUE, imagetype="png", height = 480 , width = 480 , resolution = 300, compression = "lzw", lwd = 2, lty = 'blank', fill = myCol, cex = .6, fontface = "bold", fontfamily = "sans", cat.cex = 0.6, cat.fontface = "bold", cat.default.pos = "outer", cat.pos = c(0, 0, 0, 0), cat.dist = c(0.055, 0.055, 0.085, 0.085), cat.fontfamily = "sans")
#################################################################################
