<<<<<<< HEAD
rm(list=ls())
getwd()
=======
>>>>>>> dd81f1bdf4e2b72bf0f1acd60b551cafeafc2eb6

source('./mummer align functions.R')
#USE MUMMER for pairwise alignmnet of 2 strains, check out table output and map from sry strain to reference----

<<<<<<< HEAD
str_name="MW2"
=======
str_name="Mu50"
>>>>>>> dd81f1bdf4e2b72bf0f1acd60b551cafeafc2eb6
ref_name="N315"
system("cd /data1/home/nbrayko/completed_genomes")
cmd_path="/data1/home/rpetit/bin/MUMmer/"
out_path="/data1/home/nbrayko/completed_genomes/" 

<<<<<<< HEAD
#run mummer
delta_filename=run.nucmer(str_name,ref_name,cmd_path,out_path)
align_dfs=run.summary(delta_filename,cmd_path) 
#build matrix of positions of str against qry
mum_snps=build.mummat(align_dfs$snps_df,align_dfs$aligns_df)

#run test!!!

#run char acnestral recosntruction------- 
con <- dbConnect(MySQL(), dbname="snp_staph",user="snp_select", password = "selectonly", host="peregrine.genetics.emory.edu")
dbListTables(con)
#pull 10000 bases from 200 random strains
str_base_ref= get.align.sql(con,str=200)
#create phangorn object
str_phang=phang.obj(align_mat=str_base_ref)
#create ML and parsimony tree
tree=phang.tree(str_phang)
#reconstruct bases corresponding to site patterns
anc_list=phang.anc(tree,str_phang)
#generate list of bases for ML and parsimony trees
base_list_ML=phang.bases(anc=anc_list$ML,str_base_ref,str_phang)
base_list_pars=phang.bases(anc=anc_list$ACCTRAN,str_base_ref,str_phang)

#plot tree for random char

#lookup character (FIX THIS FUNCTION)
char.lookup(char=9480,base_list_var=base_list_ML)


#12) Plot as unrooted and return list of matrices and plots of positions in the alignment-----
plotAnc(tree$ParsTree, anc_list$MPR, 8,
        show.tip.label = FALSE, 
        type="fan")

plotAnc(tree_like, anc_ml, 4,
        show.tip.label = FALSE, 
        type="fan")

=======
delta_filename=run.nucmer(str_name,ref_name,cmd_path,out_path)
align_dfs=run.summary(delta_filename,cmd_path)
mum_snps=build.mummat(align_dfs$snps_df,align_dfs$aligns_df)
>>>>>>> dd81f1bdf4e2b72bf0f1acd60b551cafeafc2eb6
