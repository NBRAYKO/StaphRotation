##COLLATING INPUT FOR `APE` PHYLOGENY PACKAGE##

####
#Load packages-----
library(RMySQL)
library(phangorn)
library(rphast)
# Set up a connection to your database management system-----
con <- dbConnect(MySQL(), dbname="snp_staph",user="snp_select", password = "selectonly", host="peregrine.genetics.emory.edu")
dbListTables(con)

#1) Get first 1000 bases form the reference-----
ref1000<-dbGetQuery(con,"SELECT * FROM `reference` WHERE `pos`<=10000")

#2) Get get the strain ids for the columns in the sample data----
strain_id<-dbGetQuery(con,"SELECT DISTINCT `id` 
                      FROM `sample`")
str_id_str= paste0("str_",strain_id[,1])
n=nrow(strain_id)


#2.1) Query for substitutions corresponsing to first 10000  positions over all ~3200 strains-----
subs<-dbGetQuery(con,"SELECT * 
                 FROM `snp` 
                 WHERE `ref_pos`<10000") #note: we don't care about snips, jsut about deletions?
#Query for the strain IDs corresponding to the subs
subs_str= dbGetQuery(con,"SELECT *
                     FROM `straintosnp` 
                     WHERE `snpid` IN (
                     SELECT `id`  
                     FROM `snp` 
                     WHERE `ref_pos`<10000)")

#merge into index of deltions
subs_coord = merge(subs, subs_str, by.x="id", by.y="snpid")

#2.2) Query for deletions fist 1000 positions over all ~3200 strains------
#NB: exclude insertions. indel has to be longer than strain-> exclude them
dels<-dbGetQuery(con,"SELECT * 
                 FROM `indel` 
                 WHERE `ref_pos`<10000 AND `ref`>`alt`")

#calculate deletion length
dels$len =nchar(dels$ref)-nchar(dels$alt)

#RANDOM: query to see what the max deletion leght is, jsut out of cutiosity... 
max_del_size=dbGetQuery(con,"SELECT  MAX(CHAR_LENGTH(`ref`)-CHAR_LENGTH(`alt`)) AS `diff`
                        FROM `indel` 
                        WHERE `ref`>`alt`")

#Query for the strain IDs corresponding to the deletions 
dels_str= dbGetQuery(con,"SELECT *
                     FROM `straintoindel` 
                     WHERE `iid` IN (
                     SELECT `id`  
                     FROM `indel` 
                     WHERE `ref_pos`<10000 AND `ref`>`alt`)")

#create object that has start and stop of each deletion for each strain
dels_coord = merge(dels, dels_str, by.x="id", by.y="iid")
dels_coord$start = dels_coord$ref_pos+1
dels_coord$stop =  dels_coord$ref_pos+dels_coord$len

#3) make test matrix of length+1x1000;------ 
##rows correspond to strain ID, cols to base pos
str_pos =matrix(nrow=n,ncol=nrow(ref1000), data=ref1000$pos, byrow=T)
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
str_base=matrix(nrow=n,nrow(ref1000), data=ref1000$base, byrow=T)
str_base[is.na(str_pos)]=NA
row.names(str_base)=str_id_str 
table(is.na(str_base))

#6) Replace position with corresonding substitution for each sample strain---------
str_base[c(paste0("str_",subs_coord$sid)),subs_coord$ref_pos]=subs_coord$alt_base

#7) Translate to codon, translate gaps as "-" -----
#Combine reference and base
str_base_ref= rbind(rbind(t(as.matrix(ref1000$base))), str_base)
row.names(str_base_ref)[1]="REF_NC317"

str_base_ref[str_base_ref=="1"]="A"
str_base_ref[str_base_ref=="2"]="C"
str_base_ref[str_base_ref=="3"]="G"
str_base_ref[str_base_ref=="4"]="T"
str_base_ref[is.na(str_base_ref)==TRUE]="-"


#8) Import as phangorn object, make treez------
#note: specify that there is a gap character!
#note:phangorn compresses data so that only non-repeating columns are stored
str_dbin=as.DNAbin(tolower(str_base_ref))
str_phang=phyDat(str_dbin)
names(str_dbin) <- names(str_phang)

#visualize distances
# dist=dist.dna(str_dbin, model="TN93")
# 
# temp <- t(as.matrix(dist))
# temp <- temp[,ncol(temp):1]
# par(mar=c(1,5,5,1))
# image(x=1:ncol(temp), y=1:nrow(temp), z=temp, col=rev(heat.colors(100)), xaxt="n", yaxt="n", xlab="",ylab="")
# axis(side=2, at=1:80, lab=rownames(str_dbin), las=2, cex.axis=.5)
# axis(side=3, at=1:80, lab=rownames(str_dbin), las=3, cex.axis=.5)

#9) Build parsimony, neighborhood join and likelihood trees, unrooted---- 
tree_pars = pratchet(str_phang, trace=0) #pratchet implements the parsimony ratchet (Nixon, 1999) and is the prefered way to search for the best tree
parsimony(tree_pars, str_phang) #returns parsimony score of optimal tree: number of base changes that have to be computed?
tree_nj=nj(dist)

#10) Use accelerated transofmration to get branch lengths-------
#(pushes characters down the tree as far as possible)
tree_pars = acctran(tree_pars, str_phang) #The funtion acctran (accelerated transformation) assigns edge length and internal nodes to the tree
parsimony(tree_pars, str_phang) #returns parsimony score of optimal tree: number of base changes that have to be computed?
tree_nj = acctran(tree_nj, str_phang)

#note: ml tree can only be calcualted afer lengths are known; see p.
tree_like = pml(tree_pars,data=str_phang, k=4)
tree_like <- optim.pml(tree_like, optGamma=TRUE,optBf=TRUE,optEdge=TRUE,optRate=TRUE)
tree_like

#11) Reconstruct ancestral bases for specific site pattern------
#reconstruct the ancestral sequences usign 2 separate parsimony-based methods
anc_acctran = ancestral.pars(tree_pars, str_phang, "ACCTRAN") 
anc_mpr = ancestral.pars(tree_pars, str_phang, "MPR")
anc_ml = ancestral.pml(tree_like)

#12) Plot as unrooted and return list of matrices and plots of positions in the alignment-----
plotAnc(tree_pars, anc_acctran, 4,
        show.tip.label = FALSE, 
        type="fan")

plotAnc(tree_like, anc_ml, 4,
        show.tip.label = FALSE, 
        type="fan")

"%w/o%" <- function(x, y) x[!x %in% y]
#get an array of lists contains all the positions in big alignment. Each list is a specific site pattern.
pos_list=lapply(c(1:attributes(anc_acctran)$nr), function(z)
  which(attributes(anc_acctran)$index==z))

#get frequency tables of the base changes for each site 
pos_tab=lapply(c(1:attributes(anc_acctran)$nr), function(z)
  table(str_base_ref[,which(attributes(str_phang)$index==z)]))

#call all the positions in the qry genome corresponding to the pattern sites in the reference alignments
pos_map_in=lapply(c(1:attributes(anc_acctran)$nr), function(z)
  ins=which(mum_snps[,"ref_pos"]%in%pos_list[[z]] ))

#call the bases in reference for which the qry genome had deletion
pos_map_out=lapply(c(1:attributes(anc_acctran)$nr), function(z)
  pos_list[[z]]%w/o%mum_snps[,"ref_pos"])      

#create array combining map_in and map_out
pos_map_array=lapply(c(1:attributes(anc_acctran)$nr), function(z){
    x=lapply(c(1: length(pos_list[[z]])), function(i) mum_snps[pos_list[[z]][i],"ref_pos"])
  return(x)})
                                  
                     
#call all the snps corresponding to the qry genome vs 
pos_snp=lapply(c(1:attributes(anc_acctran)$nr), function(z)
  mum_snps[pos_map[[z]],"snp"])

#13) get alignment of N315 w other full genomes that are on the server-----
#see the code below. need to make those into system commands
str_other=read.FASTA("/home/nbrayko/completed_genomes/outfile.fasta")
#names are a little off...
names(str_other)<- make.names(names(str_other))
#test to see how they align to the reference imported earlier
table(str_other[1,1:1000]==str_dbin[1:5,])

