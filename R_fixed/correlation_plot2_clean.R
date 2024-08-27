#-----------------------------------------------------------------------------------
# script to get heatmap of correlations of two data
# Use data: samp1("chr_cpg_loc","beta") , samp2                                                                    
# Produce correlation plot and Pearson's coefficient
# Jingru Yu
# Author: Jingru Yu
# jingruy@stanford.edu                                                                 
# 
# 2023-07-12 UTC
# Differenve between v1 and v2: grp.x and grp.y (,], now they are the same
#------------------------------------------------------------------------------------

# get cleaned data (2 cols: chr_cpg_loc, beta)

# load function
CorrelatePlot2 <- function(samp1, samp2, assay_xr, assay, covnum){
  
  # samp_compare <- merge(samp1, samp2, by=c("chr","start","end"))
  # samp_compare <- samp_compare[,c("chr","start","end","beta.x","beta.y")]
  
  samp_compare <- merge(samp1, samp2, by="chr_cpg_loc")
  samp_compare <- samp_compare[,c("chr_cpg_loc","beta.x","beta.y")]
  
  rm(samp1, samp2);gc()
  
  # heatmap
  res1 <- cor.test(samp_compare$beta.x, samp_compare$beta.y,  
                   method = "spearman", use = "complete.obs",exact = FALSE)
  
  res2 <- cor.test(samp_compare$beta.x, samp_compare$beta.y,  
                   method = "pearson", use = "complete.obs")
  
  #data do not come from a bivariate normal distribution, rank-test
  #Warning in cor.test.default(x, y, method = "spearman",exact = FALSE)
  #"pearson"
  rho <- c(res1$estimate,res2$estimate)
  
  samp_compare$grp.x <- 0
  samp_compare$grp.y <- 0
  
  for(i in 1:25){
    samp_compare$grp.x <- ifelse(samp_compare$beta.x > (i-1)/25 & samp_compare$beta.x <=(i/25),
                                 yes=i/25, no=samp_compare$grp.x)
    samp_compare$grp.y <- ifelse(samp_compare$beta.y >= (i-1)/25 & samp_compare$beta.y <(i/25),
                                 yes=i/25, no=samp_compare$grp.y)
  }

  samp_compare$grp.x <- ifelse(samp_compare$beta.x ==0, yes=0.04, no=samp_compare$grp.x)
  samp_compare$grp.y <- ifelse(samp_compare$beta.y ==1, yes=1.00, no=samp_compare$grp.y)
  
  samp_compare$grp <- paste0(samp_compare$grp.x,"_",samp_compare$grp.y)
  
  freq <- as.data.frame(table(samp_compare$grp, useNA="ifany"))
  colnames(freq)[1] <- "grp"
  
  samp_compare2 <- merge(samp_compare,freq,by="grp",all.x=TRUE)
  samp_compare2$Freq <- ifelse(is.na(samp_compare2$Freq), yes=0, no=samp_compare2$Freq)
  
  vec25 <-seq(0,1,by=0.04)
  vec25 <- vec25[-1]
  
  grp.x <- rep(vec25, each = 25) # replicate each integer in vec 9 times
  grp.y <- rep(vec25, times = 25) #circulate 5 times
  vec25m <- cbind(as.data.frame(grp.x), as.data.frame(grp.y))
  vec25m <- unique(vec25m)
  
  samp_uniq <- unique(samp_compare2[,c("grp.x","grp.y","Freq")])
  totaln <- sum(samp_uniq[, "Freq"], na.rm = TRUE)
  samp_uniq$frac <- round(samp_uniq$Freq/totaln, digits=6)
  
  samp_uniq.order <- merge(vec25m,samp_uniq,by=c("grp.x","grp.y"),all.x=TRUE)
  
  summary(samp_uniq.order$frac)
  samp_uniq.order$frac <- ifelse(samp_uniq.order$frac<=0.00001, 
                                 yes=0.00001, no=samp_uniq.order$frac)
  samp_uniq.order$frac <- ifelse(samp_uniq.order$frac>=0.10, 
                                 yes=0.10, no=samp_uniq.order$frac)
  
  if(assay =="XRBS"|assay =="RRBS"){assay2=paste0(assay,", rep. 1&2")} else {assay2=assay}
  #samp_uniq.order[is.na(samp_uniq.order)] <- 0
  #samp_uniq.order <- na.omit(samp_uniq.order)
  #make plots
  
  colorn <- c("#000071", "#0042F9","#00D4FD", "#46FF97", "#E8FD3A", "#FF5516", "#7A0606")
  
  suppressWarnings({
  p1 <- ggplot(samp_uniq.order, aes(grp.x, grp.y)) + 
    geom_tile(aes(fill = frac)) +
    scale_fill_gradientn(colors = colorn, na.value = "#000073",
                         limits = c(0.00001, 0.10), trans = "log", 
                         breaks = c(0.00001, 0.001, 0.10),
                         labels = c("0.001%", "0.1%", "10%"), 
                         name = "CpG density", 
                         guide = guide_colourbar(
                           direction = "horizontal", even.steps =FALSE,
                           frame.colour = "black", 
                           ticks.colour = "black", # you can also remove the ticks with NA
                           barwidth=10)) +
    scale_x_continuous(limits= c(0, 1),
                       breaks= c(0.02, 0.98),
                       labels = c("0","1"),
                       name = assay_xr, expand = c(0, 0)) +
    scale_y_continuous(limits= c(0, 1),
                       breaks= c(0.02, 0.98),
                       labels = c("0","1"),
                       name = assay2, expand = c(0, 0)) +
    theme_bw() + 
    theme(legend.position="bottom",
          legend.text = element_text(color="black",size=8),
          panel.grid=element_blank(), 
          panel.border=element_blank(),
          axis.title.y=element_text(size=14),
          axis.title.x=element_text(size=14),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14)
    ) 
  
  p2 <- p1 + annotate("segment", x = 0, y = 0.02, xend = 0.98, yend = 0.02) +
    annotate("segment", x = 0.02, y = 0, xend = 0.02, yend = 0.98) +
    annotate("segment", x = 0.98, y = 0, xend = 0.98, yend = 0.98) +
    annotate("segment", x = 0, y = 0.98, xend = 0.98, yend = 0.98)
  })
  
  print(rho)
  rm(p1,samp_compare,samp_compare2);gc()
  return(p2)
  
}
