## This script uses the gene expression matrix "ReadsPerGene.out.matrix". If used with other datasets
## the script should be modified accordingly. Ahmed Elewa elewa@ub.edu
## ReadsPerGene.out.matrix is produced by mapping 25 published P. waltl RNAseq libraries to genome
## assembly aPleWal.scaffolded.splitChromosomes.masked.fa. 
## The RNAseq libraries cover different body parts, life and regeneration stages (Elewa et al. 2017)

df = read.table('aPleWal.ReadsPerGene.out.matrix',header = TRUE,sep = '\t')
df = df[5:nrow(df),]

## Column names are numerical bu this table has necessary sample descriptions.
samplediscription = read.table('aPleWal.ReadsPerGene.out.matrix.samplesdisc.tab',header = TRUE,sep = '\t')


# Twenty five samples and 142,667 predicted gene models
dim(df) #[1] 142667     25

# Only 35,167 of the 142,667 predicted gene models encode a conserved protein. (i.e. gene has a UNIPROT ID)
dim(df[substr(rownames(df),1,4)!='gene',]) #[1] 35167    25

# Flag transposable elements in gene list (Optional)
TEstems =  as.list(read.table('TEstems_20220707.list',header = FALSE,sep = '\t'))
TEstatus = c()
for (i in 1:nrow(df)){
  TEstatus = c(TEstatus, unlist(strsplit(rownames(df[i,]),'.',1))[1]%in%unlist(TEstems))
}

# 26,854 conserved protein-coding genes do not encode transposable elements
dim(df[substr(rownames(df),1,4)!='gene' & TEstatus==FALSE,]) #[1] 26854    25

# Of the 26,854 genes, 3,479 are not expressed at all (putative pseudogenes)
dim(df[substr(rownames(df),1,4)!='gene' & rowSums(df)==0 & TEstatus==FALSE,]) #[1] 3479    25

# Therefore, 23,375 genes are Expressed Conserved Protein-coding genes that are not TEs and not putative pseudogenes (ECPCG)
dim(df[substr(rownames(df),1,4)!='gene' & rowSums(df)>0 & TEstatus==FALSE,]) #[1] 23375    25

# To work with this set of genes, pass it to a new variable. 
df_ECPCG = df[substr(rownames(df),1,4)!='gene' & rowSums(df)>0 & TEstatus==FALSE,]
dim(df_ECPCG) #[1] 23375    25

