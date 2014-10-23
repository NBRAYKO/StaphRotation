##COLLATING INPUT FOR `APE` PHYLOGENY PACKAGE##

####
#Load packages-----
library(RMySQL)
library(phangorn)
library(rphast)
# Set up a connection to your database management system-----

"%w/o%" <- function(x, y) x[!x %in% y]

get.align.sql =function(db=con, pos, str=200){
  #1) Get first pos=? bases in `reference` db-----
  if (is.na(pos)){
    ref_seq<-dbGetQuery(con,"SELECT * FROM `reference`")
  } else ref_seq<-dbGetQuery(con,paste0("SELECT * FROM `reference` WHERE `pos`<=",pos)) 
  
  
  #2) Get get the strain ids corresponding to the first 200 randomly chosen strains  strians----
  strain_id<-dbGetQuery(con,"SELECT DISTINCT `id` 
                        FROM `sample`")
  n=nrow(strain_id)  
  #if str=y is integer, randomly select y strains to include in query; otherwise, if a list is provided, queryt those strains
  if (is.numeric(str) & length(str==1)){
    strain_id =subset(strain_id, strain_id$id %in% sample(strain_id$id, 200, replace=FALSE))
    print(paste("Randomy drawing",str,"of", n, "strains from database"))
  } else if  (is.vector(str)) {
    strain_id =subset(strain_id, strain_id$id %in% c(str))
    print(paste("Selecting predefined list of",length(str),"of", n, "strains from database"))
  } else stop("Enter str=n for random slection of n strains, or enter list of numbers designating strain IDs to be drawn")
  str_id_str= paste0("str_",strain_id[,1])
  n=nrow(strain_id)
  #2.1) Query for substitutions corresponsing to first x pos over the set of y str -----
  if (is.na(pos)){
    subs<-dbGetQuery(con,"SELECT *  FROM `snp` ") 
  } else {
    subs<-dbGetQuery(con,paste("SELECT * 
                               FROM `snp` 
                               WHERE `ref_pos`<",pos)) #note: we don't care about snips, jsut about deletions?
  }
  
  #Query for the strain IDs corresponding to the subs in the strain_id list of queried strians
  subs_str<- dbGetQuery(con, sprintf("SELECT *
                                     FROM `straintosnp` 
                                     WHERE `snpid` IN (
                                     SELECT `id`
                                     FROM `snp` ) 
                                     AND `sid` in (%s)", paste(strain_id$id, collapse = ",")))
  
  #merge into index of deltions
  subs_coord = merge(subs, subs_str, by.x="id", by.y="snpid")
  
  #2.2) Query for deletions------
  #NB: exclude insertions. indel has to be longer than strain-> exclude them
  if (is.na(pos)){
    dels<-dbGetQuery(con,"SELECT * 
                     FROM `indel` 
                     WHERE `ref`>`alt`")
  } else dels<-dbGetQuery(con,paste("SELECT * 
                   FROM `indel` 
                                    WHERE `ref_pos`<",pos ,"AND `ref`>`alt`"))
  
  #calculate deletion length
  dels$len =nchar(dels$ref)-nchar(dels$alt)
  
  #   #RANDOM: query to see what the max deletion leght is, jsut out of cutiosity... 
  #   max_del_size=dbGetQuery(con,"SELECT  MAX(CHAR_LENGTH(`ref`)-CHAR_LENGTH(`alt`)) AS `diff`
  #                           FROM `indel` 
  #                           WHERE `ref`>`alt`")
  #   
  #Query for the strain IDs corresponding to the deletions in the strain_id list of queried strians
  dels_str= dbGetQuery(con,sprintf("SELECT *
                                   FROM `straintoindel` 
                                   WHERE `iid` IN (
                                   SELECT `id`  
                                   FROM `indel` 
                                   WHERE `ref`>`alt`)
                                   AND `sid` in (%s)", paste(strain_id$id, collapse = ",") ))
  
  #create object that has start and stop of each deletion for each strain
  dels_coord = merge(dels, dels_str, by.x="id", by.y="iid")
  dels_coord$start = dels_coord$ref_pos+1
  dels_coord$stop =  dels_coord$ref_pos+dels_coord$len
  
  #3) make test matrix of length+1x1000;------ 
  ##rows correspond to strain ID, cols to base pos
  str_pos =matrix(nrow=n,ncol=nrow(ref_seq), data=ref_seq$pos, byrow=T)
  #NB that the strain IDs do not correspond to the matrix index! creating a string so they are distinguished... 
  row.names(str_pos)=str_id_str 
  
  #4) Replace those IDs in the matrix with NA to mark deletions in those positions-----
  #could not do it with apply functions ugh :/
  for (str_num in dels_coord$sid ){
    repl_x=c(paste0("str_",str_num))
    repl_ystart=dels_coord[dels_coord$sid==str_num,]$start
    repl_ystop=dels_coord[dels_coord$sid==str_num,]$stop
    str_pos[repl_x,repl_ystart:repl_ystop]<-NA
  }      
  table(is.na(str_pos))
  
  #5) Replace positions with bases from reference in entire matrix if !is.na--------
  str_base=matrix(nrow=n,nrow(ref_seq), data=ref_seq$base, byrow=T)
  str_base[is.na(str_pos)]=NA
  row.names(str_base)=str_id_str 
  table(is.na(str_base))
  
  #6) Replace position with corresonding substitution for each sample strain---------
  str_base[c(paste0("str_",subs_coord$sid)),subs_coord$ref_pos]=subs_coord$alt_base
  
  #7) Translate to codon, translate gaps as "-" -----
  #Combine reference and base
  str_base_ref= rbind(rbind(t(as.matrix(ref_seq$base))), str_base)
  row.names(str_base_ref)[1]="REF_NC317"
  
  str_base_ref[str_base_ref=="1"]="A"
  str_base_ref[str_base_ref=="2"]="C"
  str_base_ref[str_base_ref=="3"]="G"
  str_base_ref[str_base_ref=="4"]="T"
  str_base_ref[is.na(str_base_ref)==TRUE]="-"
  
  return(str_base_ref)
  }

phang.obj=function(align_mat){
  #8) Import as phangorn object, make treez------
  #note: specify that there is a gap character!
  #note:phangorn compresses data so that only non-repeating columns are stored
  str_dbin=as.DNAbin(tolower(str_base_ref))
  str_phang=phyDat(str_dbin)
  names(str_dbin) <- names(str_phang)
  return(str_phang)
}

phang.tree=function(str_phang, ML=TRUE){
  #9) Build parsimony and likelihood trees, unrooted---- 
  tree_pars = pratchet(str_phang, trace=0) #pratchet implements the parsimony ratchet (Nixon, 1999) and is the prefered way to search for the best tree
  parsimony(tree_pars, str_phang) #returns parsimony score of optimal tree: number of base changes that have to be computed?
  
  #10) Use accelerated transofmration to get branch lengths-------
  #(pushes characters down the tree as far as possible)
  tree_pars = acctran(tree_pars, str_phang) #The funtion acctran (accelerated transformation) assigns edge length and internal nodes to the tree
  parsimony(tree_pars, str_phang) #returns parsimony score of optimal tree: number of base changes that have to be computed?
  #tree_nj = acctran(tree_nj, str_phang)
  if (ML==TRUE){
    #note: ml tree can only be calcualted afer lengths are known; see p.
    tree_like = pml(tree_pars,data=str_phang, k=4)
    tree_like <- optim.pml(tree_like, optGamma=TRUE,optBf=TRUE,optEdge=TRUE,optRate=TRUE)
    tree_like
    out_tree=list(MLTree=tree_like, ParsTree=tree_pars)
  }
  else out_tree=tree_pars 
  return(out_tree)
}

phang.anc=function(tree, str_phang){
  #11) Reconstruct ancestral bases for specific site pattern------
  #reconstruct the ancestral sequences usign 2 separate parsimony-based methods and maximum likelihood
  
  if (is.list(tree)){
    anc_acctran = ancestral.pars(tree$ParsTree, str_phang, "ACCTRAN") 
    anc_mpr = ancestral.pars(tree$ParsTree, str_phang, "MPR")
    anc_ml = ancestral.pml(tree$MLTree)
    out_anc=list(ACCTRAN=anc_acctran,MPR=anc_mpr,ML=anc_ml)
  } else{
    anc_acctran = ancestral.pars(tree, str_phang, "ACCTRAN") 
    anc_mpr = ancestral.pars(tree, str_phang, "MPR")
    out_anc=list(ACCTRAN=anc_acctran,MPR=anc_mpr)
  }
  
  return(out_anc)
}

phang.bases=function(anc, str_base_ref, str_phang){
  
  #get an array of lists contains all the positions in big alignment. Each list is a specific site pattern.
  pos_list=lapply(c(1:attributes(anc)$nr), function(z)
    which(attributes(anc)$index==z))
  #get frequency tables of the base changes for each site 
  pos_tab=lapply(c(1:attributes(anc)$nr), function(z)
    table(str_base_ref[,which(attributes(str_phang)$index==z)]))
  
  phang_bases=list(PosList=pos_list, PosTab=pos_tab)
  return(phang_bases)
}


char.lookup=function(char_ref,base_list_var){
  #for random character, look up if it belongs to a site pattern and reurn pattern number
 temp=lapply(base_list_var, function(z)
   as.character(char_ref)%in%z)
 site=which(unlist(temp)==TRUE)
 if (length(site) ==0) {
   print("Selected reference strain character does not belong to any character pattern in phylogeny")
   site=NA
   return(site)
  } 
     else{
   return(site)
  }
}

# 
# phang.mummap=function(char_qry, base_list, mum_mat){
#   #for random qry strain character, return the corresponding ref character and where it belongs
#   char_ref=as.numeric(mum_mat[char_qry,"ref_pos"])
#   snp_type=mum_mat[char_qry,"snp"]
#   if (is.na(char_ref)){
#     print("Selected querry strain character not in alignmnet")
#     out_list_empty=list(Qry_str=char_qry,RefStr=char_ref,SNP=snp_type,Site=NA,BaseTab=NA)
#     return(out_list_empty)
#   }
#   else{
#     site=char.lookup(char=char_ref,base_list_var=base_list[1])
#     out_list=list(Qry_str=char_qry,RefStr=char_ref,SNP=snp_type,Site=site,BaseTab=base_list[2][site])
#     return(out_list_empty)
#     if (is.na(site)) {    
#       print("")
#     }
#     else print("")
#     
#     base_tab=
#     
#   }
# }
# 
# 
#   
  
##print if it is in the list, if not say its 
## print to which pair this corresponds
#print if there has been a substitution




# 
# 
# ##MAP AGAINST QRY-----
# #call all the positions in the qry genome corresponding to the pattern sites in the reference alignments
# pos_map_in=lapply(c(1:attributes(anc_acctran)$nr), function(z)
#   ins=which(mum_snps[,"ref_pos"]%in%pos_list[[z]] ))
# 
# #call the bases in reference for which the qry genome had deletion
# pos_map_out=lapply(c(1:attributes(anc_acctran)$nr), function(z)
#   pos_list[[z]]%w/o%mum_snps[,"ref_pos"])      
# 
# #create array combining map_in and map_out
# pos_map_array=lapply(c(1:attributes(anc_acctran)$nr), function(z){
#   x=lapply(c(1: length(pos_list[[z]])), function(i) mum_snps[pos_list[[z]][i],"ref_pos"])
#   return(x)})
# 
# #call all the snps corresponding to the qry genome vs 
# pos_snp=lapply(c(1:attributes(anc_acctran)$nr), function(z)
#   mum_snps[pos_map[[z]],"snp"])
# 
# #13) get alignment of N315 w other full genomes that are on the server-----
# #see the code below. need to make those into system commands
# str_other=read.FASTA("/home/nbrayko/completed_genomes/outfile.fasta")
# #names are a little off...
# names(str_other)<- make.names(names(str_other))
# #test to see how they align to the reference imported earlier
# table(str_other[1,1:1000]==str_dbin[1:5,])
# 
