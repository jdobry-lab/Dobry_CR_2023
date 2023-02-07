
# figure functions -------------------------------------------------------------
jd.theme <- function(p){
 pp <- p + theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "none",
      axis.ticks = element_line(colour = "grey70", size = 0.2),
      axis.text = element_text(color = "black", size = 12),
      axis.title.x = element_text(face = "bold", color = "black", size = 12),
      axis.title.y = element_text(face = "bold", color = "black", size = 12),
      panel.grid.major = element_line(colour = "grey70", size = 0.2),
      panel.grid.minor = element_blank()
    )
 return(pp) 
}

jd.theme.classic <- function(p, leg.position = "none"){
  pp <- p + theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = leg.position,
      axis.ticks = element_line(colour = "grey70", size = 0.2),
      axis.text = element_text(color = "black", size = 12),
      axis.title.x = element_text(face = "bold", color = "black", size = 12),
      axis.title.y = element_text(face = "bold", color = "black", size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  return(pp) 
}


plot.pcoa <- function(pc, gl, xpc = 1, ypc = 2, popColors = popcolors){

unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

  gl.pcoa.plot(
  pc,
  gl,
  scale = FALSE,
  ellipse = FALSE,
  plevel = 0.95,
  pop.labels = "none",
  interactive = FALSE,
  as.pop = NULL,
  hadjust = 1.5,
  vadjust = 1,
  xaxis = xpc,
  yaxis = ypc,
  zaxis = NULL,
  pt.size = 5,
  pt.colors = popColors,
  pt.shapes = NULL,
  label.size = 1,
  axis.label.size = 1.5,
  save2tmp = TRUE,
  verbose = NULL
)

p <- gl.print.reports(1)

p <- jd.theme(p)

return(p)

}


# figure 4
ibd_fig4_fun <- function(popibd, popgl, index = "pop", popCol = karyocolors){
  
  # data
  df1a <- matrixConvert(popibd$Dgen, c("id1","id2", "gen"))
  df1b <- matrixConvert(popibd$Dgeo, c("id1","id2", "geo"))
  df1 <- left_join(df1a, df1b)
  
  popdfmeta <- popgl@other$ind.metrics[c("id", index)]
  
  df2 <- df1 %>% left_join(popdfmeta, by = c("id1" = "id")) %>% 
    left_join(popdfmeta, by = c("id2" = "id"))
  colnames(df2) <- c("id1", "id2", "gen", "geo", "pop.x", "pop.y")
  df3<-  df2 %>% 
    rowwise() %>% 
    mutate(pairs = paste(sort(c(pop.x, pop.y)), collapse = "-")) %>% 
    separate(pairs, sep = "-",into = c("p1","p2"))
  df3$pops <- factor(paste(df3$p1, df3$p2))
  
  # model grob
  m <- lm(popibd$Dgen ~ popibd$Dgeo)
  b0 <- round(m$coefficients[1],3)
  b1 <- signif(m$coefficients[2], 2)
  r2 <- round(popibd$mantel$statistic^2,2)
  mExpress <- paste0(b0 ,
                     " + ", b1)
  pval <- popibd$mantel$signif
  
  grob_m <- grobTree(textGrob(bquote(paste(italic("y = "), .(mExpress),italic(" x,  R"^2)," = ", .(r2),
                                           italic(", p = "), .(pval))),
                              x= 0.15,  
                              y= 0.05,
                              hjust=0, rot = 0,
                              gp=gpar(col="black", fontsize=9, fontface="italic",
                                      backgroun = "white")))
  
  # plot
  ggplot(df3, aes(x = geo/1000, y = gen)) + #put col in geom_point
    geom_smooth(method = "lm", formula = y~x)+
    geom_point(aes(colour = p1, fill = p2), alpha = 0.9, size = 3, pch = 21, stroke = 2)+
    #geom_jitter(aes(shape = karyo), size = 2, width = 50, height = 0.02) +
    theme_classic()+
    scale_fill_manual(values = popCol) +
    scale_colour_manual(values = popCol) +
    annotation_custom(grob_m)+
    xlab("Distance (km)") +
    ylab("1-(Proportion of Shared Alleles)") -> pop_ibd_plot; pop_ibd_plot
  ibd_plotjd <- jd.theme.classic(pop_ibd_plot)
  
  return(ibd_plotjd)
}

# figure 3

ibd_fig3_fun <- function(popibd = fig3ibd, popgl = vacpK, karyocolors){

  df1a <- matrixConvert(popibd$Dgen, c("id1","id2", "gen"))
  df1b <- matrixConvert(popibd$Dgeo, c("id1","id2", "geo"))
  df1 <- left_join(df1a, df1b)
  
popdfmeta <- popgl@other$ind.metrics[c("id", "karyo")]

df2 <- df1 %>% left_join(popdfmeta, by = c("id1" = "id")) %>% 
  left_join(popdfmeta, by = c("id2" = "id"))

df3 <- df2 %>% 
  rowwise() %>% 
  mutate(pairs = paste(sort(c(karyo.x, karyo.y)), collapse = "-")) %>% 
  separate(pairs, sep = "-",into = c("k1","k2"))


m <- lm(popibd$Dgen ~ popibd$Dgeo)
b0 <- round(m$coefficients[1],2)
b1 <- signif(m$coefficients[2], 2)
r2 <- round(popibd$mantel$statistic^2,2)
mExpress <- paste0(b0 ,
                   " + ", b1)
pval <- popibd$mantel$signif

grob_m <- grobTree(textGrob(bquote(paste(italic("y = "), .(mExpress),italic(" x,  R"^2)," = ", .(r2),
                                         italic(", p = "), .(pval))),
                            x= 0.17,  
                            y= 0.05,
                            hjust=0, rot = 0,
                            gp=gpar(col="black", fontsize=11, fontface="italic",
                                    backgroun = "white")))

#levels(df3$karyo.shape)<-c(1,1,1, 2, 1 ,1)
ggplot(df3, aes(x = geo/1000, y = gen)) + #put col in geom_point
  geom_smooth(method = "lm")+
  geom_point(aes(colour = k1, fill = k2), alpha = 0.8, 
             size = 3, pch =21, stroke =2)+
  #geom_jitter(aes(shape = karyo), size = 2, width = 50, height = 0.02) +
  theme_classic()+
  scale_fill_manual(values = karyocolors) +
  scale_colour_manual(values = karyocolors) +
  annotation_custom(grob_m) +
  xlab("Distance (km)") +
  ylab("1-(Proportion of Shared Alleles)") -> karyo_ibd_fig3



  karyo_ibd_fig3jd  <-  jd.theme.classic(karyo_ibd_fig3)
  return(karyo_ibd_fig3jd)
}


# reviewer ibd karyo
review.analysis <- function(gl, analysis = "between"){
  
  if(analysis == "within"){ ndist <- 51; diffCol = "#FFC33B"}
  if(analysis == "between"){ ndist <- 1050; diffCol = "#FF5AAF"}
  
  karyoibe <- gl.ibd(gl, distance = "propShared",
                     paircols = "pop", Dgeo_trans = "Dgeo")
  
  mat <- karyoibd$Dgeo/1000
  prop <- karyoibd$Dgen
  
  dfDist <- matrixConvert(mat, c("id1", "id2", "geo"))
  dfGen <- matrixConvert(prop, c("id1", "id2", "gen"))
  jason2 <- gl
  dfkaryo <- jason2@other$ind.metrics[,c("id", "karyo")]
  dfpop <- jason2@other$ind.metrics[,c("id", "pop")]
  
  
  
  df <- left_join(dfGen, dfDist) %>% 
    left_join(dfkaryo, by = c("id1" = "id")) %>% 
    left_join(dfkaryo, by = c("id2" = "id")) %>% 
    left_join(dfpop, by = c("id1" = "id")) %>% 
    left_join(dfpop, by = c("id2" = "id")) %>% 
    mutate(pair = paste(karyo.x,"-", karyo.y)) %>% 
    filter(geo <= ndist) %>%        # change based on question between or within
  rowwise() %>% 
    mutate(karyo = ifelse(karyo.x == karyo.y, "same", "different"),
           diffpop = pop.x == pop.y) 

  indexPost <- c(df$id1, df$id2)
  

  m <- lm(gen ~ geo, data = df)
  summary(m)
  
  m2 <- lm(gen ~  geo + karyo, data = df)
  summary(m2)
  
 # tab_model(m, m2, show.aic = TRUE, file = "temp.html")
 #  webshot("temp.html",file = paste0("model_karyo_ibd_", ndist, ".png") , 
 #  vwidth = 441, vheight = 351)
  
  aic <- AIC(m, m2)
  aic$delta <- aic$AIC - min(aic$AIC)
  
  x <- 1:ndist
  
  y1 <- predict(m, newdata = data.frame(geo = x), interval = 'confidence')
  y2 <- predict(m2, newdata = data.frame(geo = x, karyo = "same"))
  y3 <- predict(m2, newdata = data.frame(geo = x, karyo = "different"))
  
  predmodel <- data.frame(geo = rep(x,2), gen = c(y2, y3), 
                          karyo = rep(c("same", "different"), each = ndist))
  
  sm <- summary(m2)
  
  b0 <- round(m2$coefficients[1],3)
  b1 <- signif(m2$coefficients[2], 1)
  b2 <- round(m2$coefficients[3], 3)
  
  r2 <- round(sm$r.squared,2)
  tb2 <- round(sm$coefficients[3,3], 2)
  
  mExpress <- paste0(b0 ,
                     " + ", b1, " + ", b2)
  pval <- round(sm$coefficients[3,4],2)
  
  grob_m <- grobTree(textGrob(bquote(paste(italic("y = "), "β"[0], 
                                           " + ",  "β"[1], italic(" x"),"  + ",  
                                           "β"[2],italic(" same karyo"))),
                              x= 0.30,  
                              y= 0.12,
                              hjust=0, rot = 0,
                              gp=gpar(col="black", fontsize=10, fontface="italic",
                                      backgroun = "white")))
  
  
  grob_m2 <- grobTree(textGrob(bquote(paste("β"[2]," = ", .(b2), " (", italic("t = "),.(tb2),
                                            italic(", p = "), .(pval), ")")),
                               x= 0.30,  
                               y= 0.05,
                               hjust=0, rot = 0,
                               gp=gpar(col="black", fontsize=8, fontface="italic",
                                       backgroun = "white")))
  
  
  
  reviewerPlot <-ggplot(df, aes(geo, gen, fill = karyo)) +
    geom_point(alpha = 0.7, size = 4, pch = 21) +
    theme_classic() +
    geom_line(data = filter(predmodel, geo < ndist), 
              aes(x = geo, y = gen, colour = karyo),size = 1, lty = 5)+
    xlab("Distance (km)")+
    ylab("1-(Proportion of Shared Alleles)")+
    scale_fill_manual(values = c(different = diffCol, same = "black"))+
    scale_colour_manual(values = c(different = diffCol, same = "black"))+
    # theme(legend.position = "none") +
    geom_text(aes(x = 39, y = 0.05, label = c("")), cex = 5)+
    theme(legend.position = "bottom")+
    labs(fill = "Karyotype", colour = "Karyotype")+
    annotation_custom(grob_m)+
    annotation_custom(grob_m2)
  
  
  return(reviewerPlot)
  
}
