


# libraries --------------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("SNPRelate")

library(dartR) # version:    2.4.3
library(tidyverse) #         1.3.1
library(devtools) #          2.4.3
library(RColorBrewer) #      1.1-3
library(patchwork) #         1.1.1
library(scales) #            1.2.0
library(gridExtra) #         2.3
library(grid) #              4.1.2
library(raster) #            3.5-2
library(webshot) # version   0.5.4
library(otuSummary) #        0.1.1
library(sjPlot) #            2.8.11
library(spaa) #              0.2.2


# functions --------------------------------------------------------------------
source("dobry_figure_functions.R")

# setup ------------------------------------------------------------------------
## load ------------------------------------------------------------------------
snps <- read.csv('snp_varanas_acanthurus.csv')

snp_matrix <- as.matrix(snps[,-1])
rownames(snp_matrix) <- snps$X

indmeta <- read.csv('indmeta_varanas_acanthurus.csv')
locmeta <- read.csv('locmeta_varanas_acanthurus.csv')

# make genlight object ---------------------------------------------------------
data <- t(snp_matrix)
loci <- rownames(snp_matrix)
individuals <- colnames(snp_matrix)

# genlight
gl <- new("genlight", data, ploidy = 2, loc.names = loci, 
          ind.names = individuals)
gl <- gl.compliance.check(gl, verbose = 1)


# add ind meta
gl@other$ind.metrics <- gl@other$ind.metrics %>% 
  left_join(indmeta)


# add loci meta
table(gl@loc.names == locmeta$AlleleID)
gl@other$loc.metrics <- gl@other$loc.metrics %>% 
  dplyr::select(-`array(NA, nLoc(x))`) %>%
  cbind(locmeta) %>% 
  relocate(AlleleID)

table(gl@loc.names == gl@other$loc.metrics$AlleleID)

# add latlon (from meta)
gl@other$latlon <- data.frame(lat = gl@other$ind.metrics$lat,
                              lon = gl@other$ind.metrics$lon,
                              row.names = gl@other$ind.metrics$id)

# asign pop
pop(gl) <- gl@other$ind.metrics$pop
vac.pops <- gl

## filter ----------------------------------------------------------------------

vac.gl1 <- gl.filter.callrate(vac.pops, method = "loc", threshold = 0.80)
vac.gl1 <- gl.filter.callrate(vac.gl1, method = "ind", threshold = 0.75) 
vac.gl1 <- gl.filter.reproducibility(vac.gl1, t = 0.99)
vac.gl1 <- gl.filter.monomorphs(vac.gl1, v = 0)

## colours ---------------------------------------------------------------------
popcolors <- c(East = "#0571b0", North = "#404040", 
               South = "#5e3c99", West = "#e66101")

karyocolors <- c("MM" = "red", "MA" = "black", "AA" = "blue")

## pcoa ------------------------------------------------------------------------
system.time(pc <- gl.pcoa(vac.gl1, nfactors = 5))

## ibd data --------------------------------------------------------------------
# pop
indexXY <- complete.cases(vac.gl1@other$latlon) # missing locations
vac.gl17 <- gl.keep.ind(vac.gl1,ind.list = indNames(vac.gl1)[indexXY])

latAdjust <- runif(nInd(vac.gl17))/10000 # unique individual locations
vac.gl17@other$latlon$lat <- vac.gl17@other$latlon$lat+ latAdjust

# karyo
vac.gl18 <- gl.drop.pop(vac.gl17, as.pop = "karyo", pop.list = "")
pop(vac.gl18) <- vac.gl18@other$ind.metrics$karyo

# figure 3 ---------------------------------------------------------------------
popi <- levels(vac.gl1@pop)
vac.edit <- vac.gl1
# all East population assumed to be MM
vac.edit@other$ind.metrics$karyo[pop(vac.gl1)=="East"] <- "MM"


pList2 <- list()
pList3 <- list()
for(i in 1:4){
  #pcoa
  vac.popi <- gl.keep.pop(vac.edit, pop.list = popi[i], as.pop = "pop")
  vacp <- gl.filter.monomorphs(vac.popi)
  pop(vacp) <- vacp@other$ind.metrics$karyo
  vacpK <- gl.keep.pop(vacp, as.pop = "karyo", pop.list = names(karyocolors))
  pca <- gl.pcoa(vacpK, nfactors = 5)
  pList2[[i]] <- plot.pcoa(pca, vacpK, xpc = 1, ypc = 2, 
                           popColors = karyocolors)
  
  plotfiles <- list.files(tempdir())[grep("Plot", list.files(tempdir()))]
  unlink(paste0(tempdir(), "/", plotfiles))
  
  #ibd 
  fig3ibd <-  gl.ibd(vacpK, distance = "propShared", paircols = "pop",
                     Dgeo_trans = "Dgeo", save2tmp = TRUE)
  
  pList3[[i]]<- ibd_fig3_fun(popibd = fig3ibd, popgl = vacpK, karyocolors)
  
}

fig3 <- grid.arrange(pList2[[2]],pList3[[2]], pList2[[1]], pList3[[1]],
                     pList2[[3]],pList3[[3]], pList2[[4]], pList3[[4]],
                     nrow=4)

ggsave("./figure_3.png", fig3,
       units = "cm", width = 15*1.5, height = 21*1.5,
       dpi = 600)

# ggsave("./figure_3.tiff", fig3,
#        units = "cm", width = 15*1.5, height = 21*1.5,
#        dpi = 600, device = "tiff")


# figure 4 ---------------------------------------------------------------------
## pcoa figures ----------------------------------------------------------------
p1 <- plot.pcoa(pc, vac.gl1, xpc = 1, ypc = 2, popColors = popcolors)
p2 <- plot.pcoa(pc, vac.gl1, xpc = 3, ypc = 4, popColors = popcolors)

## ibd figures -----------------------------------------------------------------

# pop
popibd <- gl.ibd(vac.gl17, distance = "propShared",
                 paircols = "pop", Dgeo_trans = "Dgeo")

pop_ibd_plotjd <- ibd_fig4_fun(popibd = popibd, popgl = vac.gl17,
                       index = "pop", popCol = popcolors)

# karyo
karyoibd <- gl.ibd(vac.gl18, distance = "propShared",
                   paircols = "pop", Dgeo_trans = "Dgeo")

karyo_ibd_plot2jd <- ibd_fig4_fun(popibd = karyoibd, popgl = vac.gl18, 
                        index = "karyo", popCol = karyocolors)

## reviewer figures ------------------------------------------------------------

rplotBetween <- review.analysis(gl = vac.gl18, analysis = "between") %>% 
  jd.theme.classic(leg.position = "bottom") 


rplotWithin <- review.analysis(gl = vac.gl18, analysis = "within") %>% 
  jd.theme.classic(leg.position = "bottom") 


## FINAL figures ---------------------------------------------------------------
fig4 <- grid.arrange(p1, p2, 
                     pop_ibd_plotjd,
                     karyo_ibd_plot2jd,
                     rplotWithin,
                     rplotBetween,
                     nrow=3)


ggsave("./figure_4.png", fig4,
       units = "cm", width = 15*1.25, height = 21*1.25,
       dpi = 600)

# ggsave("./figure_4.tiff", fig4,
#        units = "cm", width = 15*1.25, height = 21*1.25,
#        dpi = 600, device = "tiff")


