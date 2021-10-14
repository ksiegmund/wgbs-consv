
read_pwd <- function(x){
  dat <- fread(x, header = F)
  colnames(dat) <- c("chr", "start", "end", "pwd", "beta.x", "beta.y")
  #dat <- dat %>% filter(chr=="chr19")
  dat <- as.data.frame(dat)
  return(dat)
}

# Annotation function
annot <- function(dat, region){
  wgpos <- GRanges(dat)
  dm_annotated <- annotate_regions(
    regions = wgpos,
    annotations = region,
    ignore.strand = TRUE,
    quiet = FALSE)
  
  dat_annot <- data.frame(dm_annotated)
  
  return(dat_annot)
}

# make chromosome region -- 100bp bins
mkchrreg <- function(howlong, whichchr){
  a <- rep(0, howlong)
  b <- rep(0, howlong)
  i <- 1
  for(i in 1:howlong){
    a[i] <- (i-1)*100+1
    b[i] <- i*100
    i+1
  }
  
  table1 <- cbind(a, b)
  colnames(table1) <- c("start", "end")
  
  chrreg <- as.data.frame(table1) %>% mutate(seqnames=whichchr, regid=row_number()) %>% 
    dplyr::select(seqnames, start, end, regid)
  
  return(chrreg)
}

# mkdata for plot
mkdata <- function(dat, whichchr, minvalue){
  
  wkdata <- dat %>% filter(ave_score>minvalue) %>% 
    rowwise() %>% 
    mutate(n_sites=max(c_across(starts_with("n_sites_")))) %>% 
    filter(n_sites>=5) %>% ungroup %>% 
    mutate(N1N2=mean(ave_pwd_N1N2),
           JAJB=mean(ave_pwd_JAJB),
           IAIB=mean(ave_pwd_IAIB),
           DADB=mean(ave_pwd_DADB),
           HAHB=mean(ave_pwd_HAHB),
           MAMB=mean(ave_pwd_MAMB)) %>% 
    distinct(N1N2, JAJB, IAIB, DADB, HAHB, MAMB) %>% 
    pivot_longer(N1N2:MAMB)
  
  nongene <- c(0.11473, 0.21769, 0.11609, 0.11547, 0.13374, 0.10755)
  nogenedat <- as.data.frame(nongene)
  
  comb <- cbind(wkdata, nogenedat)
  comb$seqnames <- whichchr
  
  return(comb)
  
}

# mktable for each cell line

mkdata1 <- function(dat, whichscore, minvalue){
  
  index <- which(colnames(dat)==whichscore)
  
  wktmp <- dat %>% filter(dat[,index]>minvalue) %>% 
    rowwise() %>% 
    mutate(n_sites=max(c_across(starts_with("n_sites_")))) %>% 
    filter(n_sites>=5) %>% ungroup %>% 
    dplyr::select(seqnames, start, end, n_sites, whichscore, starts_with("ave_pwd_"))
  
  return(wktmp)
  
}