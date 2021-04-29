
library(biomaRt)
s=read.table('Pancreas_Normal_FPKM_summary.txt',header = T,sep='\t')
ID=rownames(s)
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl") 
gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"), 
                  filters="ensembl_gene_id", values=ID, mart=human)
gene_coords$length=gene_coords$end_position - gene_coords$start_position

write.table(gene_coords,file = 'Gene_length.txt',row.names = F,quote = F,sep='\t')


#### Calculate FPKM from htseq counts

folder='/fast-data/BI/RUO_bioinfo/Pancreas_Normal_Dir/'
length=read.table('/fast-data/BI/RUO_bioinfo/Immune_profiling/Gene_length_keep.txt',header=T)

htseq_files=list.files(path=folder,pattern = '*htseq*')
names=sapply(htseq_files, function(x) strsplit(x, '.htseq')[[1]][1])

insertSize=0
cor=NULL
for (hf in htseq_files){
  htseq=read.table(gzfile(paste0(folder,hf)))
  fpkm_file=paste0(strsplit(hf,'.htseq')[[1]][1],'.FPKM.txt.gz')
  if (fpkm_file %in% list.files(path=folder)){
    fpkm=read.table(gzfile(paste0(folder,fpkm_file)))
    htseq$ID=sapply(as.character(htseq$V1),function(x) strsplit(x,'.',fixed=T)[[1]][1])
    htseq$Len=length[match(htseq$ID,length$ID),]$Length
    effLen=htseq$Len-insertSize+1
    idx=which(effLen>1)
    data=htseq[idx,]
    N=sum(data$V2)
    data$FPKM=exp( log(data$V2) + log(1e9) - log(effLen[idx]) - log(N))
    data$fpkm2=fpkm[match(data$V1,fpkm$V1),]$V2
    cor=c(cor,cor(data$FPKM,data$fpkm2))
    
  }
  
}
  
htseq=read.table(gzfile('../Pancreas_Normal_Dir/f748bf78-4dc1-47ad-8611-8186479d3e4b.htseq.counts.gz'))
sample3=read.table(gzfile('../Pancreas_Normal_Dir/dec7dc53-a480-493c-98ca-785db6a82c6d.FPKM.txt.gz'))

length=read.table('Gene_length_keep.txt',header=T)
htseq$ID=sapply(as.character(htseq$V1),function(x) strsplit(x,'.',fixed=T)[[1]][1])
htseq$Len=length[match(htseq$ID,length$ID),]$Length

effLen=htseq$Len-300+1
idx=which(effLen>1)
data=htseq[idx,]

N=sum(data$V2)
data$FPKM=exp( log(data$V2) + log(1e9) - log(effLen[idx]) - log(N))
data$fpkm2=sample3[match(data$V1,sample3$V1),]$V2
cor(data$FPKM,data$fpkm2)
plot(data$FPKM,data$fpkm2,xlim=c(0,1000),ylim=c(0,1000))

## kmeans Clutering based on xCell result
library("gplots")
result=read.table('Breast_tumor_xCell_result.txt',header=T,sep='\t')
result2=result[idx,]
kmean=kmeans(t(result2),3)

tres=as.data.frame(t(result))
tres=tres[order(kmean$cluster),]
tresnor=tres/rowSums(tres)

immu=as.data.frame(t(result2))
immu=immu[order(kmean$cluster),]

heatmap(as.matrix(tres),Rowv=NA)

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
RowSideColors=c(rep('red',table(kmean$cluster)[1]),rep('orange',table(kmean$cluster)[2]),rep("green",table(kmean$cluster)[3]))
heatmap.2(as.matrix(tres),Rowv=F,trace = "none", density.info = "none",symm=F,scale='row',
          col = my_palette,RowSideColors=RowSideColors)

heatmap.2(as.matrix(immu),Rowv=F,trace = "none", density.info = "none",symm=F,scale='none',
          col = my_palette)

## GSEA for immune cell signatures
library(GSVA)
library(NbClust)
library(readxl)
data=read_xlsx('C:/Users/lanqing.ying/Downloads/Immune escape supple tables.xlsx',sheet = 1,skip = 1)
gset.list=NULL
for (cell in unique(data$CellType)){
  list=list(cell=data[data$CellType==cell,]$ENSG)
  names(list)=cell
  gset.list=c(gset.list,list)
}

expr=read.table('Breast_tumor_FPKM_summary.txt',sep='\t',header=T)
gsva=gsva(as.matrix(expr), gset.list, method="ssgsea")
rownames(gsva)=names(gset.list)
colnames(gsva)=colnames(expr)
gsva=gsva/colMeans(gsva)
gsva=t(gsva)
kmean=kmeans(gsva,3)
res=NbClust(gsva,method = "kmeans")
gsva_c=gsva[order(kmean$cluster),]
heatmap.2(gsva,Rowv=F,trace = "none", density.info = "none",symm=F,scale='none',
          col = my_palette)

## Nbclust test to determine best number of clusters
res=NbClust(gsva,method = "kmeans")
res$Best.nc
table(res$Best.nc[1,])


## NSCLC data
tumor=read.table('NSCLC__FPKM_scaled_select.txt',header=T,sep='\t')
tumor=tumor[-which(is.na(tumor$GeneID)),]
tumor=tumor[-which(duplicated(tumor$GeneID)),]
rownames(tumor)=as.character(tumor$GeneName)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
par(mar=c(1,1,1,1))
pdf('test.pdf')
heatmap.2(as.matrix(tumor[,3:47]),Rowv=F,trace = "none", density.info = "none",symm=F,scale='none',
          col = my_palette,labCol=F)

dev.off()
