seq2mfe <- function(sequence, window = 51){
  # as input takes DNA sequence ("DNAString") and calculates secondary structure (mfe in kcal/mol)
  # for all possible sequences of length equal to window
  # returns data.frame object 
  
  require(Biostrings)
  require(stringr)

  dir.create("temp")
  on.exit(unlink("temp", recursive = T))
  
  my_ranges <- data.frame(start = 1:(length(sequence) - window +1), 
                          end = window:length(sequence),
                          mfe = NA,
                          mid = NA)
  
  for(i in 1:nrow(my_ranges)){
    my_seq <- sequence[my_ranges[i,1]:my_ranges[i,2]]
    writeXStringSet(RNAStringSet(list(seq1 = my_seq)),
                    filepath = "temp/seqs.fa")
    system("RNAfold < temp/seqs.fa > temp/out.txt")
    my_lines <- readLines("temp/out.txt")[3]
    mfe <- as.numeric(str_sub(my_lines, 
                              start = window+3,
                              end = str_length(my_lines)-1))
    my_ranges[i,3] <- mfe
    my_ranges[i,4] <- ((my_ranges[i,1]):(my_ranges[i,2]))[round(length(((my_ranges[i,1]):(my_ranges[i,2])))/2)]
  }
  unlink("temp", recursive = T)
  unlink(x = c("seq1_ss.ps"))
  
  return(my_ranges)
}



free2bind <- function(sequence1, sequence2 = "UCCUCC", window = 8, loc_free2bind = "/home/piotr/Programs/free2bind"){
  # takes DNA sequence (("DNAString", sequence1) and looks for Shine-Dalgarno sequences (UCCUCC) in all
  # possible windows of length = window using free2bind
  # as a result returns a data.frame object
  
  require(Biostrings)
  require(stringr)
  
  sequence2 <- RNAString(sequence2)
  
  dir.create("temp")
  on.exit(unlink("temp", recursive = T))
  
  my_ranges <- data.frame(start = 1:(length(sequence1) - window +1), 
                          end = window:length(sequence1),
                          deltaG = NA)
  
  for(i in 1:nrow(my_ranges)){
    temp_sequence1 <- sequence1[(my_ranges[i,1]):(my_ranges[i,2])]
    temp_rna <- RNAString(temp_sequence1)
    options( warn = -1 )
    a <- system(paste("cd ", loc_free2bind, "; ./free_align.pl ", temp_rna, " ", sequence2, sep = ""), intern = T)[8]
    deltaG <- as.numeric(str_sub(string = a, start = str_locate(a, "=")[1]+2, end = str_length(a)))
    my_ranges[i,3] <- deltaG
  }
  unlink("temp", recursive = T)
  return(my_ranges)
}





## with trinucleotydes

## with trinucleotydes
codon.cov <- function(transcripts, reads, ribo_shift = -7, length = "all"){
  # takes reads - GAlignments object 
  # and DNAStringSet object as input 
  # and calculates pausing scores (p_score) at each position
  ### ---- # ---- # ---- # ---- # ---- # ----# ---- # ---- ### 
  require(dplyr)
  print("subseting reads")
  reads <- tbl_df(reads[strand(reads) == "+"])
  if((length == "all")[1]){
  }else{reads <- reads %>% filter(width %in% length)}
  
  reads$start <- reads$end - 51 + ribo_shift
  
  # 2. calculating pausing scores for each nucleotyde
  reads_by_pos <- reads %>% group_by(seqnames, start) %>% summarise(occ = length(seqnames))
  #colnames(reads_by_pos) <- c("seqnames", "start", "occ")
  reads_by_pos <- split(reads_by_pos, reads_by_pos$seqnames)
  my_seqnames <- names(reads_by_pos)
  stopifnot(all(names(reads_by_pos) %in% names(transcripts)))
  
  print("calculating p_scores")
  # 3. subset reads_by_pos to the cds region and calculare p_score
  reads_by_pos <- lapply(seq_along(reads_by_pos), function(x){
    my_df <- reads_by_pos[[x]]
    my_name <- names(reads_by_pos[x])
    my_width <- width(transcripts[my_name])
    
    my_df <- merge(my_df, data_frame(seqnames = as.character(my_df$seqnames[1]), start = 1:my_width),
                   all = T)
    my_df <- filter(my_df, start %in% 1:my_width)
    my_df[is.na(my_df)] <- 0
    my_df <- as_data_frame(my_df)
    my_df <- mutate(my_df, p_score = occ / (sum(occ)/nrow(my_df)), 
                    mean_reads = sum(occ)/nrow(my_df))
    my_codons <- DNAs2codons(transcripts[my_name])
    stopifnot(nrow(my_df) == nrow(my_codons[[1]])*3)
    my_df$codon <- rep(my_codons[[1]]$seq, each = 3)
    
    return(my_df)
  })
  names(reads_by_pos) <- my_seqnames
  return(reads_by_pos)
}


DNAs2codons <- function(transcripts){
  # x = DNAStringSet
  my_list <- list()
  for(i in 1:length(transcripts)){
    element <- DNAStringSet(transcripts[[i]],
                            start = seq(from = 1, by = 3,  to = length(transcripts[[i]])-2),
                            width = 3)
    element <- as.data.frame(element)
    colnames(element) <- "seq"
    element$codon <- 1:nrow(element)
    
    my_list <- c(my_list, list(element))
  }
  names(my_list) <- names(transcripts)
  return(my_list)
}

se <- function(x) sqrt(var(x)/length(x))



mfe.ccf <- function(x, cutoff = 100, window = 50, bkg = FALSE){
  # input: 
  # list of data_frames with p_scores calculated for each codon and mfe at nucleotyde resolution
  # output:
  # returns list of data.frames
  
  if(as.numeric(x[1,4]) < 1){
    my_df <- data_frame(ccf = character(), 
                        codon = character(),
                        pos = character(),
                        gene = character(),
                        gene_pos = character())    
  }else{
    selected_rows <- which(x$p_score > cutoff)
    selected_rows <- selected_rows[selected_rows %in% (52):(nrow(x)-51)]
    
    if(bkg == T){
      set.seed(500)
      selected_rows <- sample((52):(nrow(x)-51), size = 10*length(selected_rows), replace = T)
    }else{
      
    }
    
    
    my_df <- data_frame(ccf = character(), 
                        codon = character(),
                        pos = character(),
                        gene = character(),
                        gene_pos = character())
    
    if(length(selected_rows) < 1){
      # do nothing
    }else{
      for(i in selected_rows){
        my_region <- x[((i-window):(i+window)),]
        my_df_temp <- data_frame(ccf =  as.numeric(ccf(my_region$p_score, -my_region$mfe, 
                                                       lag.max = window, plot = F)$acf),
                                 mfe = my_region$mfe,
                                 pos = -window:window,
                                 gene = as.character(my_region$seqnames[1]),
                                 gene_pos = as.character(i))
        my_df <- rbind(my_df, my_df_temp)
      }
    }
  }
  return(my_df)
}




f2b.ccf <- function(x, cutoff = 100, window = 50, bkg = FALSE){
  # input: 
  # list of data_frames with p_scores calculated for each codon  and deltaG at nucleotyde resolution
  # output:
  # returns list of data.frames
  
  if(as.numeric(x[1,4]) < 1){
    my_df <- data_frame(ccf = character(), 
                        deltaG = character(),
                        pos = character(),
                        gene = character(),
                        gene_pos = character())    
  }else{
    selected_rows <- which(x$p_score >= cutoff)
    selected_rows <- selected_rows[selected_rows %in% (51+1):(nrow(x)-51)]
    
    if(bkg == T){
      set.seed(500)
      selected_rows <- sample((51+1):(nrow(x)-51), size = length(selected_rows)*5, replace = F)
    }else{
      
    }
    
    my_df <- data_frame(ccf = character(), 
                        deltaG = character(),
                        pos = character(),
                        gene = character(),
                        gene_pos = character())
    
    if(length(selected_rows) < 1){
      # do nothing
    }else{
      for(i in selected_rows){
        
        my_region <- x[((i-window):(i+window)),]
        if(sum(my_region$p_score) == 0){
          
        }else{
          my_df_temp <- data_frame(ccf =  as.numeric(ccf(my_region$p_score, my_region$deltaG, 
                                                         lag.max = window, plot = F)$acf),
                                   deltaG = my_region$deltaG,
                                   pos = -window:window,
                                   gene = as.character(my_region$seqnames[1]),
                                   gene_pos = as.character(i))
          my_df <- rbind(my_df, my_df_temp)
          
        }
      }
    }
  }
  return(my_df)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
