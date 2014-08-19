test_that("mummer mapping is correct",{
  ####BUILD TEST DATASETS AND EXPECTED OUTPUTS OF align mummer####
  require(rentrez)
  require(testthat)
  source('./mummer align functions.R')
  
  #1. download the first 5000 bps of the N315 sequence from Rentrez-----
  test_raw=entrez_fetch(db="nuccore", id="NC_002745.2", rettype="fasta", from=0, to=5000)
  write(test_raw, "test_raw.fasta")
  test1=read.FASTA("test_raw.fasta") 
  
  #2. simulate 100 subs every 45 bps, delete first 500 ------
  #remove the first 500 characters; 
  names(test1) = "orig"
  subs=as.character(test1)$`orig`[501:5000]
  test1=as.character(test1)$`orig`
  #Create test sequence without DNA characters
  test_nondna=test1
  test_nondna[test_nondna=="a"]="W"
  # create function for making subs
  inv.base=function(x){
    if (x=="a") x="c"
    else if (x=="c") x="g"
    else if (x=="t") x="a"
    else if (x=="g") x="t"
    return(x)
  }
  #run function
  for (i in 1:100){
    subs[45*i-15]=inv.base(subs[45*i-15])
  }
  
  #3. add 4 dels every 1000 bps, each del is 3 bases long-----
  del_subs=subs
  for (i in 1:4){
    del_subs[i*1000]=""
    del_subs[i*1000-1]=""
    del_subs[i*1000+1]=""
  }
  
  #5. add insertion splicing the first 500 chars in middle------
  indel_subs=c(del_subs[1:500], test1[1:500],del_subs[501:4500]) 
  
  #6. save all as fasta files-----
  write(c(">gi|test\n",test1),"test1.fasta")
  test1=read.FASTA("test1.fasta")
  write(c(">gi|test_non\n",test_nondna),"test_nondna.fasta")
  test_nondna=read.FASTA("test_nondna.fasta")
  write(c(">gi|test_subs\n",subs),"test_subs.fasta")
  test_subs=read.FASTA("test_subs.fasta")
  write(c(">gi|test_dels\n",del_subs),"test_del_subs.fasta")
  test_del_subs=read.FASTA("test_del_subs.fasta")
  write(c(">gi|test_subs\n",indel_subs),"test_indel_subs.fasta")
  test_indel_subs=read.FASTA("test_indel_subs.fasta")
  
  #7. generate expected output mum_snps files for alignments:------ 
  #7.1 subs expected------
  mum_snps_test1_subs_expect=matrix(nrow=5000, ncol=3)
  colnames(mum_snps_test1_subs_expect)=c("str_pos", "ref_pos", "snp")
  mum_snps_test1_subs_expect[,"snp"]="non-align"
  mum_snps_test1_subs_expect[501:5000,"str_pos"]=501:5000
  mum_snps_test1_subs_expect[501:5000,"ref_pos"]=1:4500
  mum_snps_test1_subs_expect[501:5000,"snp"]="mum"
  mum_snps_test1_subs_expect[seq(530,5000,45),"snp"]="sub"
  #7.2.1 dels expected------
  mum_snps_test1_del_subs_expect= mum_snps_test1_subs_expect
  mum_snps_test1_del_subs_expect[c(1499:1501, 2499:2501,3499:3501,4499:4501),"ref_pos"]=NA
  mum_snps_test1_del_subs_expect[!is.na(mum_snps_test1_del_subs_expect[,"ref_pos"]),"ref_pos"]=c(rep(1:4488))
  mum_snps_test1_del_subs_expect[c(1499:1501, 2499:2501,3499:3501,4499:4501),"snp"]="del"
  #7.2.2 ins expected--------
  mum_snps_test1_ins_subs_expect= mum_snps_test1_subs_expect[1:4488,]
  mum_snps_test1_ins_subs_expect[,"str_pos"]=c(1:4488)
  mum_snps_test1_ins_subs_expect[,"ref_pos"]=c(501:1498,1502:2498,2502:3498,3502:4498,4502:5000)
  mum_snps_test1_ins_subs_expect[,"snp"]="mum"
  mum_snps_test1_ins_subs_expect[which(mum_snps_test1_ins_subs_expect[,"ref_pos"]%in%seq(530,5000,45)),"snp"]="sub"
  mum_snps_test1_ins_subs_expect[c(998,1995,2992,3989),"snp"]="ins_s"
  mum_snps_test1_ins_subs_expect[c(999, 1996,2993,3990),"snp"]="ins_e"
 
  #8. the actual tests--------
  system("cd /data1/home/nbrayko/completed_genomes/")
  cmd_path="/data1/home/rpetit/bin/MUMmer/"
  out_path="" 
  
  #8.1 identical sequences: mummer should throw its own default error------
  delta_filename=run.nucmer("test1", "test1" ,cmd_path,out_path)
  expect_that(run.summary(delta_filename,cmd_path), throws_error())
  
  #8.2 non-DNA: should throw error coded by the align mummer functions.R--------
  expect_that(run.nucmer("test1", "test_nondna" ,cmd_path,out_path), throws_error("Ref not a DNA sequence"))
  expect_that(run.nucmer( "test_nondna", "test1" ,cmd_path,out_path), throws_error("Seq not a DNA sequence"))
  
  #8.3 test1 agains subs------
  delta_filename=run.nucmer("test1", "test_subs" ,cmd_path,out_path)
  align_dfs=run.summary(delta_filename,cmd_path)
  mum_snps_test1_subs=build.mummat(align_dfs$snps_df,align_dfs$aligns_df) 
  ##expect output to match mock expected output
  expect_equal(mum_snps_test1_subs,mum_snps_test1_subs_expect)
  ##expect there to be 4501 unique values in firt column  (1:4500+NA)
  expect_that(length(unique(mum_snps_test1_subs[,1])),equals(4501))
  
  #8.4 test2 against dels--------
  delta_filename=run.nucmer("test1", "test_del_subs" ,cmd_path,out_path)
  align_dfs=run.summary(delta_filename,cmd_path)
  mum_snps_test1_del_subs=build.mummat(align_dfs$snps_df,align_dfs$aligns_df) 
  ##expect output to match mock expected output
  expect_equal(mum_snps_test1_del_subs,mum_snps_test1_del_subs_expect)
  ##expect there to be 4501 unique values in firt column  (1:4500+NA)
  expect_that(length(unique(mum_snps_test1_del_subs[,1])),equals(4501) )
  
  #8.5 test1 agains insertions---------
  delta_filename=run.nucmer("test_del_subs","test1" ,cmd_path,out_path)
  align_dfs=run.summary(delta_filename,cmd_path)
  mum_snps_test1_ins_subs=build.mummat(align_dfs$snps_df,align_dfs$aligns_df) 
  ##expect output to match mock expected output
  expect_equal(mum_snps_test1_ins_subs,mum_snps_test1_ins_subs_expect)
  ##expect there to be 4501 unique values in firt column  (1:4500-12 skips)
  expect_that(length(unique(mum_snps_test1_ins_subs[,1])),equals(4488) )
  
  #8.6 test1 agains indels---------
  delta_filename=run.nucmer("test_indel_subs" , "test1" ,cmd_path,out_path)
  align_dfs=run.summary(delta_filename,cmd_path)
  mum_snps_test1_indels=build.mummat(align_dfs$snps_df,align_dfs$aligns_df) 
  
}