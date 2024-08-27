#-----------------------------------------------------------------------------------
# script to GetIntersectCG_fixList based on GetIntersectCG, intersect XR with blocks 
# to get new beta values 
# Use data: sample.txt, allBlocks_CpG (ref_blocks2k_3cpgs.CpG_loc_230603.bed), 
#   array list(tcga_60k, cns_32k)
# Produce data: sample_original_prep.RData, sample_original_beta.RData
# Author: Jingru Yu
# jingruy@stanford.edu                                                                 
# 
# 2023-08-01 UTC
#------------------------------------------------------------------------------------ 

# Function
GetIntersectCG_fixList <- function(TFs, TFs_after, sampXRfixed_interBlocks, 
                           top_CG, samp_tp, maxcov, orig, pathtemp, interFiles){
  
  #1) load sample
  samp <- fread(
    file.path(pathtemp, paste0(TFs,TFs_after,".txt")),
    header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  colnames(samp) <- c("chr","start","end","meth.cnt","total") #14M
  samp <- samp %>% filter(total>=5)
  
  #2) merge sample with all CpGs in blocks (optional: ArrayInterBlock)
  # to get new beta values of blocks
  
  samp_block <- merge(samp, sampXRfixed_interBlocks, by=c("chr","start","end"),all.y=T)
  #count NA to 0
  samp_block[is.na(samp_block)] <- 0
  #sum cpgs in the same intersected scMarker
  samp_block_sum <- samp_block[!duplicated(samp_block[,"chr_loc"]),] %>% 
    group_by(chr_marker_loc) %>%
    summarise_at(vars(meth.cnt,total),sum) %>%
    mutate(beta=meth.cnt/total) 
  #pay attention, deduplicate by one variable is much faster
  
  samp_block_sum <- samp_block_sum %>% #1.6M
    filter(total!=0) %>%
    mutate(chr_marker = sapply(str_split(chr_marker_loc, "_",  n = 3), `[`, 1),
           start_marker = as.numeric(sapply(str_split(chr_marker_loc, "_",  n = 3), `[`, 2)),
           end_marker = as.numeric(sapply(str_split(chr_marker_loc, "_",  n = 3), `[`, 3)))
  
  
  # merge individual CpG locations by marker_loc, keep all CpGs in the markers
  beta1 <- merge(sampXRfixed_interBlocks[,1:6],samp_block_sum,
                 by=c("chr_marker","start_marker","end_marker"), all.y=T)
  
  beta1 <- beta1 %>%
    mutate(chr=ifelse(is.na(chr), yes=chr_marker, no=chr),
           start=ifelse(is.na(start), yes=start_marker, no=start),
           end=ifelse(is.na(end), yes=end_marker, no=end))
  
  #3) merge individual CpGs from sample with re-calculated beta values in the blocks
  # to include CpGs that are not in the blocks
  
  beta2 <- merge(beta1,samp,by=c("chr","start","end"),all.x=T,all.y=T) #16M
  
  beta2 <- beta2 %>%
    mutate(chr_marker=ifelse(is.na(chr_marker), yes=chr, no=chr_marker),
           start_marker=ifelse(is.na(start_marker), yes=start, no=start_marker),
           end_marker=ifelse(is.na(end_marker), yes=end, no=end_marker),
           chr_cpg_loc = paste0(chr,"_",start,"_",end),
           chr_marker_loc = paste0(chr_marker,"_",start_marker,"_",end_marker),
           
           meth.cnt.x=ifelse(is.na(meth.cnt.x), yes=meth.cnt.y, no=meth.cnt.x),
           total.x=ifelse(is.na(total.x), yes=total.y, no=total.x),
           beta=ifelse(is.na(beta), yes=meth.cnt.y/total.y, no=beta))
  
  #4) filter by the max cov in each block
  
  beta_maxcov <- beta2 %>%
    group_by(chr_marker_loc) %>%
    summarise(max = max(total.x, na.rm=TRUE))
  beta_maxcov <- beta_maxcov %>% filter(max>=maxcov)
  
  beta_maxcov_filter <- 
    beta2[(beta2$chr_marker_loc %in% beta_maxcov$chr_marker_loc),] #10M
  
  # save(beta_maxcov_filter, 
  #      file=file.path(pathtemp,interFiles,paste0(TFs, "-interblocks_allvars.RData")))
  
  #5) filter CpGs that are included in the array probe list, TCGA/CNS classifier
  
  if(samp_tp=="TCGA60k"){
    # merge CGs by CpG_loc
    beta_merge_tp <- merge(top_CG,beta_maxcov_filter,
                           by.x=c("chr_cg","start_cg","end_cg"),
                           by.y=c("chr","start","end"),all.x=T)
    
    samp_inter_beta <- beta_merge_tp %>% 
      filter(!is.na(CG), !is.na(chr_marker_loc)) %>% dplyr::select(CG,beta)
    
    save(samp_inter_beta, 
         file=file.path(pathtemp,interFiles,paste0(sub("_[^_]+$", "", TFs),
         "-",maxcov,"cov.",samp_tp,"_",orig,"_beta.RData")))
    
    samp_inter_marker <- beta_merge_tp %>% 
      filter(!is.na(CG), !is.na(chr_marker_loc)) %>% 
      dplyr::select(CG,chr_marker_loc,beta)
    samp_inter_beta_prep <- beta_merge_tp[(beta_merge_tp$chr_marker_loc %in% 
                                             samp_inter_marker$chr_marker_loc),]
    
    save(samp_inter_beta_prep, 
         file=file.path(pathtemp,interFiles,paste0(sub("_[^_]+$", "", TFs),
         "-",maxcov,"cov.",samp_tp,"_",orig,"_prep.RData")))
    
  } else{
    if(samp_tp=="CNS32k"){
      # merge CGs by CpG_loc
      beta_merge_tp <- 
        merge(top_CG,beta_maxcov_filter,by.x=c("chr_cg","start_cg","end_cg"),
              by.y=c("chr","start","end"),all.x=T)
      
      samp_inter_beta <- beta_merge_tp %>% 
        filter(!is.na(CG), !is.na(chr_marker_loc)) %>% dplyr::select(CG,beta)
      save(samp_inter_beta, 
           file=file.path(pathtemp,interFiles,paste0(sub("_[^_]+$", "", TFs),
                          "-",maxcov,"cov.",samp_tp,"_",orig,"_beta.RData")))
      
      samp_inter_marker <- beta_merge_tp %>% 
        filter(!is.na(CG), !is.na(chr_marker_loc)) %>% 
        dplyr::select(CG,chr_marker_loc,beta)
      samp_inter_beta_prep <- beta_merge_tp[(beta_merge_tp$chr_marker_loc %in% 
                                               samp_inter_marker$chr_marker_loc),]
      
      save(samp_inter_beta_prep, 
           file=file.path(pathtemp,interFiles,paste0(sub("_[^_]+$", "", TFs),
           "-",maxcov,"cov.",samp_tp,"_",orig,"_prep.RData")))
      
    } else { 
      beta_maxcov_filter_uniq <- 
        beta_maxcov_filter[!duplicated(beta_maxcov_filter$chr_cpg_loc),] %>%
        filter(!is.na(total.y)) %>% 
        dplyr::select(chr_cpg_loc,beta,meth.cnt.y, total.y)
      
      
      save(beta_maxcov_filter_uniq, 
           file=file.path(pathtemp,interFiles,paste0(sub("_[^_]+$", "", TFs), 
           "-",maxcov,"cov.blocks_noExtraCpG_",orig,"_prep.RData")))
    }
  }
  print(TFs)
  rm(samp,samp_block,beta1,beta2,beta_maxcov,beta_maxcov_filter);gc()
}




