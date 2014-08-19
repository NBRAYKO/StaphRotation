
source('./mummer align functions.R')
#USE MUMMER for pairwise alignmnet of 2 strains, check out table output and map from sry strain to reference----

str_name="Mu50"
ref_name="N315"
system("cd /data1/home/nbrayko/completed_genomes")
cmd_path="/data1/home/rpetit/bin/MUMmer/"
out_path="/data1/home/nbrayko/completed_genomes/" 

delta_filename=run.nucmer(str_name,ref_name,cmd_path,out_path)
align_dfs=run.summary(delta_filename,cmd_path)
mum_snps=build.mummat(align_dfs$snps_df,align_dfs$aligns_df)