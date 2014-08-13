
rm(list=ls())
#USE MUMMER for pairwise alignmnet of 2 strains, check out table output and map from sry strain to reference----

str_name="TCH1516"
ref_name="N315"
system("cd /data1/home/nbrayko/completed_genomes")


#### 1) run NUCmer w following options #####
# --mum Use anchor matches that are unique in both the reference and query;
system(paste0("/data1/home/rpetit/bin/MUMmer/nucmer --mum --prefix=ref_qry /data1/home/nbrayko/completed_genomes/", 
              ref_name, ".fasta /data1/home/nbrayko/completed_genomes/",
              str_name,".fasta"))

# 1.1) filter NUCmer w following options------
# -q  Query alignment using length*identity weighted LIS. For each query, leave only the alignments which form the longest consistent set for the query
# -r	Reference alignment using length*identity weighted LIS. For each reference, leave only the alignments which form the longest consistent set for the reference.
# -u float	Set the minimum alignment uniqueness, i.e. percent of the alignment matching to unique reference AND query sequence [0, 100], (default 0)
# -o float  Set the maximum alignment overlap for -r and -q options as a percent of the alignment length [0, 100], (default 75)

system("/data1/home/rpetit/bin/MUMmer/delta-filter -q -r -o 0 ref_qry.delta > ref_qry_f.delta")

#### 2) generate summary of the alignments ######
# -g nly display alignments included in the Longest Ascending Subset, i.e. the global alignment. Recommened to be used in conjunction with the -r or -q options. Does not support circular sequences
# -l includes seq length
# -T tab-delimited
#-q sort by query length (since we are interested in the qry strain as teh starting block of ancestraal reconstruciton further down the pipeline)
system("/data1/home/rpetit/bin/MUMmer/show-coords  -g -H -l -T -q ref_qry_f.delta > ref_qry.coords")

# 2.1 )load aligns-----
aligns_df=read.table("ref_qry.coords",blank.lines.skip=TRUE)
names(aligns_df)= c("ref_pos","ref_end", "str_pos",  "str_end", "ref_al_len", "str_al_len", "pct_match", "ref_len", "str_len" , "str_name","ref_name")
ref_fasta_name=levels(aligns_df$ref_name)[1]
str_fasta_name=levels(aligns_df$str_name)[1]


#### 3) generate summary of SNPs ######
#-q sort by query strain.
#-H drop header, make it tab-delimited (-T), include sequence length
system("/data1/home/rpetit/bin/MUMmer/show-snps -H -T -q -l ref_qry_f.delta > ref_qry.snps")

# 3.1) load snps-----
snps_df=read.table("ref_qry.snps",blank.lines.skip=TRUE,)
snps_df=snps_df[,1:6]
names(snps_df)= c("ref_pos","ref_sub", "str_sub",  "str_pos", "len", "dist")

# [str_pos] position of the SNP in the reference sequence . 
#   NB! For indels, this position refers to the 1-based position of the first character before the indel, e.g. for an indel at the very beginning of a sequence this would report 0. For indels on the reverse strand, this position refers to the forward-strand position of the first character before indel on the reverse-strand, e.g. for an indel at the very end of a reverse complemented sequence this would report 1. 
# [str_sub and ref_sub] character or gap at this position in the reference [SUB] character or gap at this position in the query 
# [ref_sub] position of the SNP in the query sequence 
# [len] distance from this SNP to the !nearest! mismatch (end of alignment, indel, SNP, etc) in the same alignment 
# [dist] distance from this SNP to the nearest sequence end 


#### 4) show aligns ######
system(paste0("/data1/home/rpetit/bin/MUMmer/show-aligns ref_qry_f.delta > ref_qry.aligns '", str_fasta_name, "' '",ref_fasta_name,"'"))

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
  mum_snps[temp_align_s:temp_align_e,"snp"]="non-align"
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
        mum_snps[temp_s,"snp"]=snps_df_temp[i,"snp"]
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
        mum_snps[temp_s,"snp"]=snps_df_temp[i,"snp"]
      }
    }
  }
} 


# 
# 
# #try mummer
# system("/data1/home/rpetit/bin/MUMmer/run-mummer3  /data1/home/nbrayko/completed_genomes/TCH1516.fasta /data1/home/nbrayko/completed_genomes/N315.fasta  bla")
# #parse the mummer output file. Mummer returns alignment info for reverse sequence, split at every >;
# system("csplit -s -f part bla.errorsgaps  '/^>/' {*}")
# 
# #just grab the first one----
# mum_df=read.table("part01",skip=1,blank.lines.skip=TRUE)
# 
# 
# 
# 
# #note: the reference strain in mummer (1st col of output) is strain X in this script, so we call it str_pos
# names(mum_df)= c("str_pos", "ref_pos", "mum_len", "errs", "str_gap", "ref_gap", "n_diff")
# 
# #get the stop positions in the ref
# mum_df$str_end=mum_df$str_pos+mum_df$mum_len-1
# mum_df$ref_end=mum_df$ref_pos+mum_df$mum_len-1
# #get the length of gaps
# 
# #data frame of MUMs where first column (and the row index) is the position of X str and second is corresponding reference. 
# #the mismatches are NA
# mums=matrix(nrow=max(mum_df$str_end), ncol=2)
# colnames(mums)=c("str_pos", "ref_pos")
# 
# for (mrow in 1:nrow(mum_df) ){
#   mums[mum_df[mrow,]$str_pos:mum_df[mrow,]$str_end
#        ,"str_pos" ]=mum_df[mrow,]$str_pos:mum_df[mrow,]$str_end
#   
#   mums[mum_df[mrow,]$str_pos:mum_df[mrow,]$str_end
#        ,"ref_pos"  ]=mum_df[mrow,]$ref_pos:mum_df[mrow,]$ref_end
# } 
# 
# #data frame of gaps. of the str_gap and ref_gap are same lenght, it was just a substitution. if ref is longer, a deletion
# mum_df$str_gpos=mum_df$str_pos-as.numeric(as.character(mum_df$str_gap))
# mum_df$str_gend=mum_df$str_gpos+as.numeric(as.character(mum_df$str_gap))-1
# mum_df$ref_gpos=mum_df$ref_pos-as.numeric(as.character(mum_df$ref_gap))
# mum_df$ref_gend=mum_df$ref_gpos+as.numeric(as.character(mum_df$ref_gap))-1
# 
# 
# 
# 
# 
# gaps=matrix(nrow=max(mum_df$str_end), ncol=2)
# colnames(gaps)=c("str_pos", "ref_pos")
# 
# 
# #create an array the length of x strain
# 
# 
# 
# 
# # [enter password d******d]
# # 
# # 1  passwd
# # 2  dir
# # 3  mugsy
# # 4  MUGSY
# # 5  exit
# # 6  mugsy
# # 7  mugsyWGA
# # 8  ls
# # 9  completed_genomes
# # 10  cd completed_genomes
# # 11  ls
# # 12  mugsyWGS
# # 13  mugsyWGA
# # 14  cat N315.fasta MSSA476.fasta > align.fasta 
# # 15  mugsyWGA
# # 16  -seq align.fasta
# # 17  mugsyWGA -seq align.fasta
# # 18  mugsyWGA --help
# # 19  mugsyWGA --seq align.fasta
# # 20  ls
# # 21  head outfile
# # 22  head outfile.maf
# # 23  rm Eco-K12-NC_000913.fasta 
# # 24  rm Kp-NC_009648.fasta 
# # 25  rm PAO1.fasta 
# # 26  ls
# # 27  cat N315.fasta MSSA476.fasta MW2.fasta Mu50.fasta TCH1516.fasta > align.fasta 
# # 28  mugsyWGA --seq align.fasta
# # 29  mugsyWGA --seq align.fasta -f msf
# # 30  ls
# # 31  grep "gi|" outfile.maf | awk '{print ">" $2 "\n" $7}' > outfile.fasta
# 
# #create array of mums as long as the x strain; replace position with reference strain mummer position
# 
