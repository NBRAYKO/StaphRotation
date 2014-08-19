function(str_name,ref_name="N315", cmd_path,out_path ){
  require(ape)
  #### Check that the inputs are all DNA sequences
  if (unique(names(table(as.character(read.FASTA(paste0(out_path, str_name,".fasta")))))!= c("a", "c", "g" ,"t"))) stop("Str not a DNA sequence")
  if (unique(names(table(as.character(read.FASTA(paste0(out_path, ref_name,".fasta")))))!= c("a", "c", "g" ,"t"))) stop("Ref not a DNA sequence")
  #### 1) run NUCmer w following options #####
  # --mum Use anchor matches that are unique in both the reference and query;
  return_filename=paste0("nucmer_",ref_name,"_",str_name)
  system(paste0(cmd_path,"nucmer --mum --prefix=", return_filename," ",
                out_path, ref_name, ".fasta ",
                out_path, str_name,".fasta"))
  
  # 1.1) filter NUCmer w following options------
  # -q  Query alignment using length*identity weighted LIS. For each query, leave only the alignments which form the longest consistent set for the query
  # -r  Reference alignment using length*identity weighted LIS. For each reference, leave only the alignments which form the longest consistent set for the reference.
  # -u float  Set the minimum alignment uniqueness, i.e. percent of the alignment matching to unique reference AND query sequence [0, 100], (default 0)
  # -o float  Set the maximum alignment overlap for -r and -q options as a percent of the alignment length [0, 100], (default 75)
  
  system(paste0("/data1/home/rpetit/bin/MUMmer/delta-filter -q -r -o 0 ",
                return_filename,".delta > ",
                return_filename,"_f",".delta"))
  
  ret_filename=paste0(return_filename,"_f")
  return(paste0(ret_filename))
}

function(ret_filename,cmd_path){
  #### 2) generate summary of the alignments ######
  # -g nly display alignments included in the Longest Ascending Subset, i.e. the global alignment. Recommened to be used in conjunction with the -r or -q options. Does not support circular sequences
  # -l includes seq length
  # -T tab-delimited
  #-q sort by query length (since we are interested in the qry strain as teh starting block of ancestraal reconstruciton further down the pipeline)
  system(paste0(cmd_path, "show-coords  -g -H -l -T -q ", 
                ret_filename,".delta > ",ret_filename,".coords"))
  
  # 2.1 )load aligns-----
  aligns_df=read.table(paste0(ret_filename,".coords"),blank.lines.skip=TRUE)
  names(aligns_df)= c("ref_pos","ref_end", "str_pos",  "str_end", "ref_al_len", "str_al_len", "pct_match", "ref_len", "str_len" , "str_name","ref_name")
  ref_fasta_name=levels(aligns_df$ref_name)[1]
  str_fasta_name=levels(aligns_df$str_name)[1]
  
  #### 3) generate summary of SNPs ######
  #-q sort by query strain.
  #-H drop header, make it tab-delimited (-T), include sequence length
  system(paste0(cmd_path, "show-snps -H -T -q -l ", 
                ret_filename,".delta > ",ret_filename,".snps"))
  
  # 3.1) load snps-----
  snps_df=read.table(paste0(ret_filename,".snps"),blank.lines.skip=TRUE,)
  snps_df=snps_df[,1:6]
  names(snps_df)= c("ref_pos","ref_sub", "str_sub",  "str_pos", "len", "dist")
  
  #### 4) show aligns ######
  system(paste0(cmd_path,"show-aligns ",
                ret_filename,".delta > ",
                ret_filename,".aligns '", str_fasta_name, "' '",ref_fasta_name,"'"))
  
  return(list(
    snps_df=snps_df,
    aligns_df=aligns_df,
    aligns_filename=paste0(ret_filename,".aligns")
  ))
}


function(snps_df, aligns_df){
  ###########
  #5)create matrix of alignment where columns are the aligned positions of strain vs reference, and snp=0 for mum, 1 for subs, 2 for indel------
  mum_snps=matrix(nrow=max(aligns_df$str_len), ncol=3)
  colnames(mum_snps)=c("str_pos", "ref_pos", "snp")
  mum_snps[,"snp"]="non-align"
  
  ###########
  #6) Create matrix [mum_snps] aligning the qry to the reference strain, and indicating SNPs -------
  #loop over the aligns_df list of alignments and the snps_df list of SNPs
  
  #for each alignmnet in the aligns_df list (going down the qry strain):
  
  for (a in 1:nrow(aligns_df)){
    #get start and end of alignment
    temp_align_s=aligns_df[a,]$str_pos
    temp_align_e=aligns_df[a,]$str_end
    #fill out the snp status variable to default to "mum"
    mum_snps[temp_align_s:temp_align_e,"snp"]="mum"
    
    #fill out the qry strain positions in matrix
    mum_snps[temp_align_s:temp_align_e,"str_pos"]=temp_align_s:temp_align_e
    # cut up the SNPs df
    snps_df_temp =  snps_df[which(snps_df$str_pos>=temp_align_s & snps_df$str_pos<= temp_align_e),]
    if (nrow(snps_df_temp) ==0) {
      warning(paste("Ambiguous alignment between qry strain pos",temp_align_s,temp_align_e,", no mapping produced" )) 
    } else print(paste("Aligning qry strain pos",temp_align_s,"-", temp_align_e ))
    
    #make string variable describing type of SNP
    snps_df_temp$snp = paste0(snps_df_temp$str_sub,">",snps_df_temp$ref_sub)
    #fill positions up to first snp of that alignment:
    mum_snps[temp_align_s:(snps_df_temp[1,]$str_pos),"str_pos"]= temp_align_s:(snps_df_temp[1,]$str_pos)
    diff= temp_align_s-aligns_df[a,]$ref_pos
    mum_snps[temp_align_s:(snps_df_temp[1,]$str_pos),"ref_pos"]= temp_align_s:(snps_df_temp[1,]$str_pos)-diff
    mum_snps[!is.na(mum_snps)&mum_snps<=0]=NA
    
    #for each SNP in the alignment:
    for (i in 1:nrow(snps_df_temp) ) {   
      #if the SNP is a substitution:
      if (snps_df_temp[i,]$ref_sub!="." & snps_df_temp[i,]$str_sub!="."  ) {  
        #fill in positions for the qry strain for basepairs between current and previous sub/indel
        temp_diff=snps_df_temp[i,]$str_pos - snps_df_temp[i,]$ref_pos
        #get start position of mum gap for str
        temp_s=snps_df_temp[i,]$str_pos
        #get the end. if the end of alignment is reached, the temp_end of the mum is temp_align_e
        if (i==nrow(snps_df_temp)) temp_e = temp_align_e else temp_e=snps_df_temp[i+1,]$str_pos 
        mum_snps[temp_s:temp_e,"str_pos"]=temp_s:temp_e
        #the ref position is the query - the shift between the frames
        mum_snps[temp_s:temp_e,"ref_pos"]=(temp_s:temp_e)-temp_diff 
        #mark as snp
        #mum_snps[temp_s,"snp"]=snps_df_temp[i,"snp"]
        mum_snps[temp_s,"snp"]="sub"
      }
      
      #if SNP is a deletion in qry string (insertion)
      if (snps_df_temp[i,]$ref_sub!="." & snps_df_temp[i,]$str_sub=="."  ) {
        #fill in positions for the qry strain for basepairs between current and previous sub/indel
        temp_diff=snps_df_temp[i,]$str_pos - snps_df_temp[i,]$ref_pos
        temp_s=snps_df_temp[i,]$str_pos+1
        #if the end of alignment is reached, the temp_end of the mum is temp_align_e
        if (i==nrow(snps_df_temp)) temp_e = temp_align_e+1 else temp_e=snps_df_temp[i+1,]$str_pos+1
        mum_snps[temp_s:temp_e,"str_pos"]=temp_s:temp_e
        #the ref position is the query - the shift between the frames
        mum_snps[temp_s:temp_e,"ref_pos"]=(temp_s:temp_e)-temp_diff     
        mum_snps[(temp_s-1),"snp"]="ins_s"
        mum_snps[(temp_s),"snp"]="ins_e"
      }  
      
      if (snps_df_temp[i,]$ref_sub=="." & snps_df_temp[i,]$str_sub!=".") {
        #fill in positions for the qry strain for basepairs between current and previous sub/indel
        temp_diff=snps_df_temp[i,]$str_pos - snps_df_temp[i,]$ref_pos
        #get start position of mum gap for str
        temp_s=snps_df_temp[i,]$str_pos
        #get the end. if the end of alignment is reached, the temp_end of the mum is temp_align_e
        if (i==nrow(snps_df_temp)) temp_e = temp_align_e else temp_e=snps_df_temp[i+1,]$str_pos  
        mum_snps[temp_s:temp_e,"str_pos"]=temp_s:temp_e
        #the ref position is the query - the shift between the frames
        mum_snps[(temp_s):temp_e,"ref_pos"]=c( NA,(temp_s+1):(temp_e)-temp_diff )
        #change the snp variable to mark as deletion in reference string
        #mum_snps[temp_s,"snp"]=snps_df_temp[i,"snp"]
        mum_snps[temp_s,"snp"]="del"
      }
    }
  }
  return(mum_snps)
}
