#-----------------------------------------------------------------------------------
# script to sample.TSNE.InterBlocks based on RealsampXR_array_inter,  
# to get TSNE plots 
# Use data: sample_beta.RData, ref_hg38_betas (TCGA_hg38_betas_filter,CNStumor_hg38_betas_filter),
# anno_ref (anno_TCGA.sub1, anno_CNStumor3)
# Produce files: sample.tnse.pdf, sample_mcf_Tsne.pdf, sample_Tsne.csv
# Author: Jingru Yu
# jingruy@stanford.edu                                                                 
# 
# 2023-12-17 UTC
# Update: 2024-02-29
#------------------------------------------------------------------------------------ 

library(Rtsne)
library(RSpectra)
library(ggrepel)

sample.TSNE.InterBlocks <- 
  function(pathtemp, TFs, ref_hg38_betas, anno_ref, Output, filetype, samp_tp){
  ## Loop over samples in the loop
  
  for(i in 1:length(TFs)){
    
    load(file=file.path(pathtemp,Output,paste0(TFs[i],"_beta.RData")))
    colnames(samp_inter_beta) <- c("CG",TFs[i])
    
    betas_m <- merge(ref_hg38_betas, samp_inter_beta, by="CG") 
    #30k CGs for TCGA, and 16k for CNStumor with intersection
    #chk_na <- sum(is.na(betas_m[,2470])) / nrow(betas_m) #3%
    
    hg38_betas_mat <- data.matrix(betas_m, rownames.force = NA) #matrix is faster
    
    betas_trans <- t(hg38_betas_mat)
    betas_trans <- betas_trans[-1,] #twice if CG and position
    ## Remove any NA
    betas_trans <- betas_trans[,!colSums(!is.finite(betas_trans))]
    
    ## get final annotation file
    anno3 <- data.frame(TCGA_ID = TFs[i],
                        Sentrix_ID = TFs[i],
                        material = "Frozen",
                        TCGA_Project = ifelse(
                          substr(sub(".*_(.*)", "\\1", TFs[i]),1,3)=="sub",
                          yes="Purified",no="Original"),
                        TCGA_Project_subgrp = ifelse(
                          substr(sub(".*_(.*)", "\\1", TFs[i]),1,3)=="sub",
                          yes="Purified",no="Original"),
                        TCGA_Project_mcf = ifelse(
                          substr(sub(".*_(.*)", "\\1", TFs[i]),1,3)=="sub",
                          yes="Purified",no="Original"),
                        shape_grp=2,
                        purity.abs=0)
    anno_final <- rbind(anno_ref, anno3)
    
    
    if ( ncol(betas_trans) > 200) {
      # calculate first 94 PCs
      pca <- prcomp_svds(betas_trans,k=94)
    } else {
      pca <- prcomp_svds(betas_trans,k=floor(ncol(betas_trans)/3*2))
    }
    
    # set.seed(6924)
    # set.seed(724538)
    # set.seed(31245)
    # set.seed(12453)
    set.seed(123)
    
    # calculate tSNE
    res <- Rtsne(pca$x, pca=FALSE, max_iter=2500,theta=0.0,verbose=T)
    #remove seed
    rm(.Random.seed, envir=globalenv())
    
    shape_class <- factor(anno_final$shape_grp, order=T,
                          levels=c(1,2), labels=c("TCGA/CNS","New sample"))
    
    if(samp_tp=="TCGA60k"){
      
      meth_class <- as.factor(anno_final$TCGA_Project_subgrp)
      levels(meth_class) #43
      
      tsne_plot.df <- 
        data.frame(x = res$Y[,1], y = res$Y[,2],
                   cluster = factor(meth_class, labels=as.vector(levels(meth_class))),
                   shape.f = factor(shape_class, labels=as.vector(levels(shape_class))),
                   Sentrix_ID = rownames(pca$x))
      
      #save plot data
      write.csv(tsne_plot.df, file=file.path(pathtemp,Output,
        paste0(TFs[i],"_",filetype,"Tsne.csv")),row.names = FALSE)
      
      #different color pattern        
      meth3 <- 
        factor(tsne_plot.df$cluster, order=T,
               levels=c("ACC","BLCA","BRCA_B","BRCA_L","CESC_AD",
                        "CESC_SCC","COAD", "DLBC","ESCA_AD","ESCA_SCC",
                        "GBM","HNSC","KICH","KIRC","KIRP",
                        "KIRP_CIMP","LAML","LGG", "CHOL","LIHC",
                        "Low_fraction","LUAD","LUSC","MESO","OV",
                        "PAAD","PNET","PRAD","READ","SKCM",
                        "STAD","TGCT_NS","TGCT_S","THCA","THYM",
                        "UCEC_AD", "UCEC_Other","UCS", "UVM","UVM_MP",
                        "CTRL (BLOOD)","CTRL (MUS)","CTRL (REA)","CTRL (HighT)","CTRL (N)",
                        "Original","Purified")) #"BRCA_Other"
      #"THCA_Classical","THCA_CpGislandMeth","THCA_Follicular",
      hc.norm = hclust(dist(tsne_plot.df[,c(1,2)]))
      tsne_plot.df$hclust = meth3
      hc.norm.cent = tsne_plot.df %>% group_by(hclust) %>% dplyr::select(x, y) %>% 
        summarise(x=median(x), y=median(y))
      
      col_code <- c(
        "Original"="#000000","Purified"="#000000","Low_fraction"="#a8a3a3",
        "ACC"="#0b4068","BLCA"="#ff6dfa","BRCA_B"="#c95f7b","BRCA_L"="#6e1e33",
        "CESC_AD"="#68ab97","CESC_SCC"="#764e90","COAD"="#ad6ef1",
        "DLBC"="#e4703a","ESCA_AD"="#7af870","ESCA_SCC"="#11a0aa","GBM"="#44aeff",
        "HNSC"="#55c056","KICH"="#b90f38","KIRC"="#9d5723","KIRP"="#fe576c",
        "KIRP_CIMP"="#94946f","LAML"="#fe576c","LGG"="#006ead","CHOL"="#fa8fcc",
        "LIHC"="#b37dec","LUAD"="#b890e0","LUSC"="#27780d","MESO"="#aa4a1c",
        "OV"="#0581be","PAAD"="#734ba9","PNET"="#a37c92","PRAD"="#f5811d",
        "READ"="#f0c507","SKCM"="#63632f","STAD"="#f79246","TGCT_NS"="#cd4d00",
        "TGCT_S"="#ddca11","THCA"="#87037a",
        "THYM"="#3c1b46","UCEC_Other"="#f50a41","UCEC_AD"="#43dad7",
        "UCS"="#6da7d9","UVM"="#ecbf00","UVM_MP"="#8d7200",
        "BRCA"="#c299a4","UCEC"="#0af597",
        "CTRL (BLOOD)"="#67686c", "CTRL (MUS)"="#7c7e82", "CTRL (REA)"="#cacbcd",
        "CTRL (HighT)"="#e8a9cc", "CTRL (N)"="#c57d95")
      #"BRCA_Other"="#c299a4"
      #"THCA_Classical"="#87037a",
      #"THCA_CpGislandMeth"="#677864","THCA_Follicular"="#dbb3d7",
      pdf(file=file.path(pathtemp,Output,
                         paste0(TFs[i],"_",filetype,"Tsne.pdf")), 
          width = 10, height = 10)
      
      print(labeled_tsne <-
              ggplot(data=tsne_plot.df, aes(x=x, y=y, color = hclust, label=TRUE))+
              xlab("t-SNE 1") + ylab("t-SNE 2") +
              geom_point(aes(shape=shape.f)) + 
              scale_shape_manual(values=c(20,2)) + 
              theme_bw() +
              geom_label_repel(aes(label = hclust), data=hc.norm.cent, 
                               label.size = 0.02, max.overlaps = Inf) + 
              scale_colour_manual(values = col_code, labels=meth3) +
              ggtitle(paste0("Sample, ", TFs[i])) +
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(colour = "black", fill=NA),
                    legend.position="none"))
      #labeled_tsne
      dev.off()  
    } else {
      
      meth_class <- as.factor(anno_final$TCGA_Project_subgrp)
      levels(meth_class) #43
      
      tsne_plot.df <- 
        data.frame(x = res$Y[,1], y = res$Y[,2],
                   cluster = factor(meth_class, labels=as.vector(levels(meth_class))),
                   shape.f = factor(shape_class, labels=as.vector(levels(shape_class))),
                   Sentrix_ID = rownames(pca$x))
      
      #save plot data
      write.csv(tsne_plot.df, 
                file=file.path(pathtemp,Output,
                paste0(TFs[i],"_",filetype,"Tsne.csv")),row.names = FALSE)
      
      meth3 <- factor(tsne_plot.df$cluster, order=T,
                      levels=c("Original","Purified",
                               "A IDH", "A IDH, HG", "ANA PA",  "ATRT, MYC", "ATRT, SHH",
                               "ATRT, TYR", "CHGL", "CHORDM", "CN", "CNS NB, FOXR2",
                               "CONTR, ADENOPIT", "CONTR, CEBM", "CONTR, HEMI", "CONTR, HYPTHAL", "CONTR, INFLAM",
                               "CONTR, PINEAL", "CONTR, PONS", "CONTR, REACT","CONTR, WM",  "CPH, ADM",
                               "CPH, PAP","DLGNT", "DMG, K27", "EFT, CIC","ENB, A",
                               "ENB, B",  "EPN, MPE", "EPN, PF A", "EPN, PF B", "EPN, RELA",
                               "EPN, SPINE","EPN, YAP","ETMR", "EWS","GBM, G34",
                               "GBM, MES", "GBM, MID","GBM, MYCN", "GBM, RTK I","GBM, RTK II",
                               "GBM, RTK III", "HGNET, BCOR", "HGNET, MN1", "HMB", "IHG",
                               "LGG, DIG/DIA","LGG, DNT","LGG, GG", "LGG, MYB","LGG, PA MID",
                               "LGG, PA PF", "LGG, PA/GG ST", "LGG, RGNT", "LGG, SEGA", "LIPN",
                               "LYMPHO", "MB, G3", "MB, G4", "MB, SHH CHL AD",
                               "MB, SHH INF", "MB, WNT", "MELAN",   "MELCYT",  "MNG",
                               "O IDH",   "PGG, nC", "PIN T,  PB A","PIN T,  PB B","PIN T, PPT",
                               "PITAD, ACTH", "PITAD, FSH LH", "PITAD, PRL","PITAD, STH DNS A", "PITAD, STH DNS B",
                               "PITAD, STH SPA", "PITAD, TSH", "PITUI", "PLASMA",  "PLEX, AD",
                               "PLEX, PED A", "PLEX, PED B", "PTPR, A", "PTPR, B", "PXA",
                               "RETB", "SCHW", "SCHW, MEL", "SFT HMPC","SUBEPN, PF",
                               "SUBEPN, SPINE", "SUBEPN, ST"))
      
      hc.norm = hclust(dist(tsne_plot.df[,c(1,2)]))
      tsne_plot.df$hclust = meth3
      hc.norm.cent = tsne_plot.df %>% group_by(hclust) %>% select(x, y) %>% 
        summarise(x=median(x), y=median(y))
      
      col_code <- c("Original"="#000000","Purified"="#000000",
                    "GBM, G34"="#55c056",  "DMG, K27"="#a3f187",  "ATRT, SHH"="#009cea",
                    "CONTR, CEBM"="#bdbdbd", "GBM, MYCN"="#48e948", "LGG, PA PF"="#9d73c9" ,     
                    "CONTR, REACT"="#d7d7d7", "LGG, PA MID"="#b890e0",  "LGG, RGNT"="#f79246",
                    "MB, WNT"="#6da7d9", "ATRT, MYC"="#006ead", "ATRT, TYR"="#44aeff",
                    "LGG, PA/GG ST"="#764e90",    "LGG, SEGA"="#ad6ef1", "MB, G4"="#9dc3f5",   
                    "MB, SHH INF"="#8bc0ff",  "MB, SHH CHL AD"="#2e80c6",   "MB, G3"="#0581be",   
                    "EPN, RELA"="#d84c4b", "PIN T,  PB A"="#cbf5af",     "SUBEPN, PF"="#dd052e",      
                    "PTPR, B"="#bfff9c",   "SUBEPN, ST"="#ff4156","EPN, YAP"="#fe576c", 
                    "EPN, PF A"="#e30149", "EPN, PF B"="#df3854", "EFT, CIC"="#955c9d", 
                    "ETMR"="#1693d6",      "CNS NB, FOXR2"="#518ccb",    "LYMPHO"="#3c1b46",   
                    "HGNET, BCOR"="#739aca",  "GBM, MES"="#538533",  "GBM, RTK II"="#7af870",     
                    "GBM, RTK I"="#68c62d","LGG, MYB"="#9159de",  "CONTR, HEMI"="#939393",     
                    "HGNET, MN1"="#ab94c1","GBM, MID"="#27780d",  "HMB"="#ed5cfc",      
                    "EPN, MPE"="#b90f38",  "SCHW"="#ffe1b2",      "O IDH"="#ffe336",    
                    "A IDH, HG"="#ddca11", "MNG"="#d566d7","A IDH"="#ecbf00",    
                    "LGG, GG"="#e4703a",   "CN"="#b85400", "SUBEPN, SPINE"="#fa3b44",   
                    "PIN T,  PB B"="#b6dc95", "PXA"="#7c46c9","ANA PA"="#734ba9",   
                    "CONTR, INFLAM"="#dbdbdb", "PIN T, PPT"="#b4f79b","CPH, ADM"="#72d4c8", 
                    "CPH, PAP"="#43dad7",  "ENB, A"="#aa4a1c",    "PITUI"="#58d8db", 
                    "CONTR, PINEAL"="#696969",    "PGG, nC"="#f5811d",   "LGG, DNT"="#cb5216",
                    "CHGL"="#b37dec",  "MELAN"="#0b4068",     "PLEX, AD"="#7b3e0a", 
                    "ENB, B"="#ff8339",    "LIPN"="#de7b2c",      "EPN, SPINE"="#ff044c",      
                    "PTPR, A"="#aede8a",   "SFT HMPC"="#ff6dfa",  "PLEX, PED A"="#9d5723",  
                    "GBM, RTK III"="#78e435",     "IHG"="#d0bbd8", "MELCYT"="#102a46",   
                    "DLGNT"="#d1723a",     "PITAD, ACTH"="#137d75",      "PITAD, STH DNS B"="#14c5c9",
                    "PITAD, PRL"="#00b1a9","PITAD, FSH LH"="#44a199",    "PLEX, PED B"="#754027",     
                    "EWS"="#ce69ec","SCHW, MEL"="#f7ca8a", "CONTR, ADENOPIT"="#d3d3d3", 
                    "LGG, DIG/DIA"="#cd4d00",     "PITAD, STH SPA"="#41ccbe",   "PITAD, STH DNS A"="#01c3bd",
                    "CONTR, WM"="#a8a8a8", "PLASMA"="#632f63",    "CHORDM"="#fa8cff",   
                    "RETB"="#ff9755", "CONTR, PONS"="#545454", "CONTR, HYPTHAL"="#7e7e7e",  
                    "PITAD, TSH"="#12d7c7")
      
      pdf(file=file.path(pathtemp,Output, 
                         paste0(TFs[i],"_",filetype,"Tsne.pdf")), 
          width = 12, height = 12) 
      
      print(labeled_tsne <-
              ggplot(data=tsne_plot.df, aes(x=x, y=y, color = hclust, label=TRUE))+
              xlab("t-SNE 1") + ylab("t-SNE 2") +
              geom_point(aes(shape=shape.f), size=1) + 
              scale_shape_manual(values=c(20,2)) + 
              theme_bw() +
              geom_label_repel(aes(label = hclust), data=hc.norm.cent, 
                               label.size = 0.015, max.overlaps = Inf) + 
              scale_colour_manual(values = col_code, labels=meth3) +
              ggtitle(paste0("Sample, ", TFs[i])) +
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(colour = "black", fill=NA),
                    legend.position="none"))
      #labeled_tsne
      dev.off() 
      
      ###########################
      meth_class_mcf <- as.factor(anno_final$TCGA_Project)
      levels(meth_class) #43
      
      tsne_plot.df <- 
        data.frame(x = res$Y[,1], y = res$Y[,2],
                   cluster = factor(meth_class_mcf, labels=as.vector(levels(meth_class_mcf))),
                   shape.f = factor(shape_class, labels=as.vector(levels(shape_class))),
                   Sentrix_ID = rownames(pca$x))
      
      #save plot data
      write.csv(tsne_plot.df, 
                file=file.path(pathtemp,Output,
                paste0(TFs[i],"_",filetype,"mcf_Tsne.csv")),row.names = FALSE)
      
      meth3_mcf <- factor(tsne_plot.df$cluster, order=T,
                          levels=c("MCF IDH GLM", "MCF ENB","MCF MB G3G4", 
                                   "MCF MB SHH", "MCF GBM", "MCF PA","MCF PLEX T", "MCF SCHW", 
                                   "ANA PA", "ATRT, MYC", "ATRT, SHH", "ATRT, TYR",
                                   "CHGL", "CHORDM", "CN", "CNS NB, FOXR2",
                                   "CONTR, ADENOPIT", "CONTR, CEBM", "CONTR, HEMI", "CONTR, HYPTHAL", "CONTR, INFLAM",
                                   "CONTR, PINEAL", "CONTR, PONS", "CONTR, REACT", "CONTR, WM", 
                                   "CPH, ADM", "CPH, PAP","DLGNT", "DMG, K27", "EFT, CIC",
                                   "EPN, MPE", "EPN, PF A", "EPN, PF B", "EPN, RELA",
                                   "EPN, SPINE","EPN, YAP","ETMR", "EWS","GBM, G34",
                                   "HGNET, BCOR", "HGNET, MN1", "HMB", "IHG", #40
                                   "LGG, DIG/DIA","LGG, DNT","LGG, GG", "LGG, MYB",
                                   "LGG, RGNT", "LGG, SEGA", "LIPN", "LYMPHO",
                                   "MB, WNT", "MELAN",   "MELCYT",  "MNG",
                                   "PGG, nC", "PIN T,  PB A","PIN T,  PB B","PIN T, PPT",
                                   "PITAD, ACTH", "PITAD, FSH LH", "PITAD, PRL","PITAD, STH DNS A", "PITAD, STH DNS B",
                                   "PITAD, STH SPA", "PITAD, TSH", "PITUI", "PLASMA",
                                   "PTPR, A", "PTPR, B", "PXA",
                                   "RETB", "SFT HMPC","SUBEPN, PF",
                                   "SUBEPN, SPINE", "SUBEPN, ST", "Original","Purified")) #78
      
      hc.norm = hclust(dist(tsne_plot.df[,c(1,2)]))
      tsne_plot.df$hclust = meth3_mcf
      hc.norm.cent = tsne_plot.df %>% group_by(hclust) %>% dplyr::select(x, y) %>% 
        summarise(x=median(x), y=median(y))
      
      col_code_mcf <- c("Original"="#000000","Purified"="#000000",
                        "MCF GBM"="#55c056",  "DMG, K27"="#a3f187", 
                        "ATRT, MYC"="#006ead", "ATRT, SHH"="#009cea","ATRT, TYR"="#44aeff",
                        "CONTR, CEBM"="#bdbdbd",  "CONTR, REACT"="#d7d7d7",  "LGG, RGNT"="#f79246",
                        "MB, WNT"="#6da7d9", "MCF PA"="#764e90",    "LGG, SEGA"="#ad6ef1",    
                        "MCF MB SHH"="#8bc0ff",   "MCF MB G3G4"="#0581be",  "MCF IDH GLM"="#ecbf00", 
                        "EPN, RELA"="#d84c4b", "PIN T,  PB A"="#cbf5af",     "SUBEPN, PF"="#dd052e",      
                        "PTPR, B"="#bfff9c",   "SUBEPN, ST"="#ff4156","EPN, YAP"="#fe576c", 
                        "EPN, PF A"="#e30149", "EPN, PF B"="#df3854", "EFT, CIC"="#955c9d", 
                        "ETMR"="#1693d6",      "CNS NB, FOXR2"="#518ccb",    "LYMPHO"="#3c1b46",   
                        "HGNET, BCOR"="#739aca", "LGG, MYB"="#9159de",  "CONTR, HEMI"="#939393",     
                        "HGNET, MN1"="#ab94c1", "HMB"="#ed5cfc", "EPN, MPE"="#b90f38",  
                        "MNG"="#d566d7",  "GBM, G34"="#27780d",
                        "LGG, GG"="#e4703a",   "CN"="#b85400", "SUBEPN, SPINE"="#fa3b44",   
                        "PIN T,  PB B"="#b6dc95", "PXA"="#7c46c9","ANA PA"="#734ba9",   
                        "CONTR, INFLAM"="#dbdbdb", "PIN T, PPT"="#b4f79b","CPH, ADM"="#72d4c8", 
                        "CPH, PAP"="#43dad7",  "MCF ENB"="#aa4a1c",    "PITUI"="#58d8db", 
                        "CONTR, PINEAL"="#696969",    "PGG, nC"="#f5811d",   "LGG, DNT"="#cb5216",
                        "CHGL"="#b37dec",  "MELAN"="#0b4068", "LIPN"="#de7b2c",     
                        "EPN, SPINE"="#ff044c",  "PTPR, A"="#aede8a",   "SFT HMPC"="#ff6dfa",  
                        "MCF PLEX T"="#9d5723",  "IHG"="#d0bbd8", "MELCYT"="#102a46",   
                        "DLGNT"="#d1723a",     "PITAD, ACTH"="#137d75",      "PITAD, STH DNS B"="#14c5c9",
                        "PITAD, PRL"="#00b1a9","PITAD, FSH LH"="#44a199",   "EWS"="#ce69ec",   
                        "MCF SCHW"="#f7ca8a", "CONTR, ADENOPIT"="#d3d3d3", "PITAD, TSH"="#12d7c7",
                        "LGG, DIG/DIA"="#cd4d00",     "PITAD, STH SPA"="#41ccbe",   "PITAD, STH DNS A"="#01c3bd",
                        "CONTR, WM"="#a8a8a8", "PLASMA"="#632f63",    "CHORDM"="#fa8cff",   
                        "RETB"="#ff9755", "CONTR, PONS"="#545454", "CONTR, HYPTHAL"="#7e7e7e"  
      )
      
      
      pdf(file=file.path(pathtemp,Output, 
                         paste0(TFs[i],"_",filetype,"mcf_Tsne.pdf")), 
          width = 12, height = 12) 
      
      print(labeled_tsne <-
              ggplot(data=tsne_plot.df, aes(x=x, y=y, color = hclust, label=TRUE))+
              xlab("t-SNE 1") + ylab("t-SNE 2") +
              geom_point(aes(shape=shape.f), size=1) + 
              scale_shape_manual(values=c(20,2)) + 
              theme_bw() +
              geom_label_repel(aes(label = hclust), data=hc.norm.cent, 
                               label.size = 0.015, max.overlaps = Inf) + 
              scale_colour_manual(values = col_code_mcf, labels=meth3_mcf) +
              ggtitle(paste0("Sample, ", TFs[i])) +
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(colour = "black", fill=NA),
                    legend.position="none"))
      labeled_tsne
      dev.off() 
    }
    
    cg_num <- data.frame(BF_id=TFs[i],interFixedXR_CpG=dim(samp_inter_beta)[1],
                         impute_CpG=0)
    write.csv(cg_num, file.path(pathtemp, Output, 
                                paste0(TFs[i],"_",filetype,"cpg.csv")),row.names = FALSE)
    
    print(TFs[i])
  }
  
} 

