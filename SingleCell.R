library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

# Load the dataset
dir='Y:\\01-New Y Drive\\Bioinformatics\\Lanqing\\BioPharma\\Single Cell\\20018-01\\filtered_feature_bc_matrix'
dir='Y:\\01-New Y Drive\\Bioinformatics\\Lanqing\\BioPharma\\Single Cell\\1729-15\\filtered_feature_bc_matrix'
pbmc.data <- Read10X(data.dir = dir)

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "test", min.cells = 1, min.features = 1)

# The percentage of reads that map to the mitochondrial genome
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")


# QC metrics
head(pbmc@meta.data)
tail(pbmc@meta.data)

# Visualize QC metrics as a violin plot
group=sapply(rownames(pbmc@meta.data),function(x) strsplit(x,'-')[[1]][2])
pbmc@meta.data$groups=group
p=VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by='groups',combine = F)

p1=p[[1]] + labs(x='sample',title='nGene') + theme(legend.position = "none") 
p2=p[[2]] + labs(x='sample',title='nUMI') + theme(legend.position = "none")
p3=p[[3]] + labs(x='sample',title='% mito') + theme(legend.position = "none")

ggbld <- ggplot_build(p1)
p1_y=ggbld$layout$coord$labels(ggbld$layout$panel_params)[[1]]$y.major_source
p1_y=c(p1_y,0,5500,NA)

p1.1=p1 + geom_hline(yintercept = 5500, color='red',linetype="dashed",size=1) + 
  geom_hline(yintercept = 0, color='red',linetype="dashed",size=1) +
  scale_y_continuous(breaks=p1_y) + 
  geom_text(aes(x=1,y=6000),label='99%',color='red')

p2.1=p2 + geom_hline(yintercept = 110000, color='red',linetype="dashed",size=1) + 
  geom_hline(yintercept = 0, color='red',linetype="dashed",size=1) +
  scale_y_continuous(breaks=c(0,50000,100000,150000,200000,110000)) + 
  geom_text(aes(x=1,y=170000),label='99%',color='red')

p3.1=p3 + geom_hline(yintercept = 85, color='red',linetype="dashed",size=1) + 
  geom_hline(yintercept = 0, color='red',linetype="dashed",size=1) +
  scale_y_continuous(breaks=c(0,25,50,75,85)) + 
  geom_text(aes(x=1,y=100),label='95%',color='red')

figure=p1.1+p2.1+p3.1

### Cumulative percent
s1=pbmc@meta.data
df1=data.frame(nFeature=s1$nFeature_RNA[order(s1$nFeature_RNA)],cell=c(1:dim(s1)[1]))
df1$cumulative_percentage= 100*cumsum(df1$cell)/sum(df1$cell)

df2=data.frame(nFeature=s1$nCount_RNA[order(s1$nCount_RNA)],cell=c(1:dim(s1)[1]))
df2$cumulative_percentage= 100*cumsum(df2$cell)/sum(df2$cell)

df3=data.frame(nFeature=s1$percent.mt[order(s1$percent.mt)],cell=c(1:dim(s1)[1]))
df3$cumulative_percentage= 100*cumsum(df3$cell)/sum(df3$cell)

df=df3
ggplot(data=df, aes(y=nFeature,x=cumulative_percentage))+
  geom_point()+
  ylab('Number of genes')+
  geom_vline(xintercept = 95, color = "red")+
  geom_text(aes(x=95, label="95%", y=0),color='red')

df[df$nFeature>85 & df$nFeature<86,]

### subset data
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- subset(pbmc, subset = groups==1)

# normalize the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

head(pbmc[["RNA"]]@scale.data)[,1:5]

# Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.2)

pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE)

pbmc = RunTSNE(pbmc)
DimPlot(pbmc, reduction = "tsne", label = TRUE)

## Color by sample
barcode=names(Idents(pbmc))
pbmc@meta.data$sample = sapply(barcode, function(x) strsplit(x, '-')[[1]][2])
DimPlot(pbmc, reduction = "tsne", label = F,group.by='sample')

## subset by sample and print plots and tables
require(lattice)
for (i in levels(factor(pbmc@meta.data$sample))) {
  s=subset(pbmc,subset = sample==i)
  outImage=paste0('1729-20-aggre3-cluster-sample',i,'-cellType.pdf')
  pdf(outImage)
  print(DimPlot(s, reduction = "tsne", label = TRUE, group.by='cell_type'))
  dev.off()
}

for (i in levels(factor(pbmc@meta.data$sample))) {
  s=subset(pbmc,subset = sample==i)
  s.markers=FindAllMarkers(s, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.2)
  meanCounts=NULL
  for (j in 1:dim(s.markers)[1]){
    gene=s.markers$gene[j]
    expression=s[["RNA"]]@scale.data[gene,]
    clust=s.markers$cluster[j]
    mc = mean(expression[s$seurat_clusters==as.character(clust)])
    meanCounts=c(meanCounts,mc)
  }
  s.markers$MeanCounts=meanCounts
  outTable=paste0('1729-20-aggre2-sample',i,'-markers.txt')
  write.table(s.markers,file=outTable,quote = F, row.names = F)

}


# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# Take all cells in cluster 2, and find markers that separate cells in the 'g1' group (metadata variable 'group')
#markers <- FindMarkers(pbmc, ident.1 = "g1", group.by = 'group', subset.ident = "2")

sample1.markers <- FindMarkers(pbmc, ident.1 = "1", ident.2 = '2', group.by = 'sample',subset.ident = '1')
head(x = sample1.markers)


# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)
top_markders=pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene)


### Add mean counts for pbmc.markers
pbmc.markers$gene=as.character(pbmc.markers$gene)
meanCounts=NULL

for (i in 1:dim(pbmc.markers)[1]){
  gene=pbmc.markers$gene[i]
  expression=pbmc[["RNA"]]@scale.data[gene,]
  #clust=pbmc.markers$cluster[i]
  rowCounts=NULL
  for (clust in levels(pbmc.markers$cluster)){
    mc = mean(expression[pbmc$seurat_clusters==clust])
    rowCounts=c(rowCounts,mc)
  }
  meanCounts=rbind(meanCounts,rowCounts)
}

for (g in unique(pbmc.markers$gene)){
  expression=pbmc[["RNA"]]@scale.data[g,]
  rowCounts1=sapply(levels(DE_genes$cluster),function(x) mean(expression[pbmc$seurat_clusters==x]))
  meanCluster=rbind(meanCluster,rowCounts1)
  
  rowCounts2=sapply(levels(DE_genes$cell_type),function(x) mean(expression[pbmc$cell_type==x]))
  meanCellType=rbind(meanCellType,rowCounts2)
  
}
colnames(meanCounts)=paste0('Expression in Cluster ',levels(pbmc.markers$cluster))
rownames(meanCounts)=rownames(pbmc.markers)

pbmc_new_markers=cbind(pbmc.markers,meanCounts)
write.csv(pbmc_new_markers,file=paste0(project,'-DE-genes.csv'),row.names = F,quote = F)

### Reformat table

pbmc.markers$gene=as.character(pbmc.markers$gene)
df=NULL
for (g in unique(pbmc.markers$gene)){
  gtable=pbmc.markers[pbmc.markers$gene==g,]
  row=NULL
  names=NULL
  for (clust in levels(pbmc.markers$cluster)){
    names=c(names,paste0('MeanCounts_Cluster',clust),paste0('logFC_Cluster',clust),paste0('pval_Cluster',clust),paste0('pval_adj_Cluster',clust))
    if (clust %in% gtable$cluster){
      logFC=gtable[gtable$cluster==clust,]$avg_logFC
      pval=gtable[gtable$cluster==clust,]$p_val
      pval_adj=gtable[gtable$cluster==clust,]$p_val_adj}
    else {logFC=NA; pval=NA; pval_adj=NA}
      row=c(row,meanCounts,logFC,pval, pval_adj)
  }
  df=rbind(df,row)
}
colnames(df)=names
write.csv(df,file=paste0(project,'-table4.csv'),row.names = T,quote = F)

#############################################
### Do clustered heatmap
z <- DoHeatmap(pbmc, features = top10$gene) + NoLegend()
y <- z$data %>% tidyr::drop_na()
df <- y %>% group_by(Identity) %>% dplyr::select(Feature, Cell, Identity, Expression)%>% tidyr::spread(key = Feature, value = Expression)


df=df[order(df$Identity),]
m=apply(df[,3:ncol(df),], 1,as.numeric)
rownames(m)=colnames(df)[3:ncol(df)]

library(gplots)
tb=table(df$Identity)
colorbar=NULL
for (i in 1:length(tb)) {
  colorbar=c(colorbar,rep(colors()[i*5],tb[i]))
  
}

heatmap.2(m,scale = "none",Colv=FALSE, dendrogram='row',labCol = FALSE,col=redblue(100),
          trace = "none", density.info = "none", ColSideColors=colorbar)

#####################################################

# Annotate cluster by markers given
markers=list('Neutrophils'= c('Ly6g', 'Orm1'),'Macrophages'=c('Adgre1', 'Ms4a7'),'T-cells'=c('Cd3e','Cd8','Cd4'))

file='/fast-data/BI/RUO_bioinfo/singleCell/1729-20/Cell_markers.txt'
Markers=read.table(file,sep='\t',header=T)
markers=sapply(as.character(Markers$Markers), function(x) strsplit(x,', ')[[1]])
names(markers)=as.character(Markers$Cell_Type)

# For each cell type
cell_cluster=NULL
for ( n in 1:length(markers)){
  table=NULL
  for (m in markers[n][[1]]){
    sub=pbmc.markers[pbmc.markers$gene==m,]
    if (dim(sub)[1]==0){
      logFC=0
      cluster='NA'
      #print (m)
      table=rbind(table,c(m,cluster,logFC))
    }
    else{
      for (i in 1:nrow(sub)){
        logFC=sub$avg_logFC[i]
        cluster=as.character(sub$cluster)[i]
        table=rbind(table,c(m,cluster,logFC))
      }
    }
    
  }
  agg=aggregate(as.numeric(table[,3]),list(cluster=table[,2]),mean)
  agg$cell_type=names(markers[n])
  cell_cluster=rbind(cell_cluster,agg)
  
}

### Deal with situation when one cell type is assigned to multiple clusters
cluster_anno=NULL
for (cl in levels(factor(pbmc.markers$cluster))){
  data=cell_cluster[cell_cluster$cluster==cl,]
  if (dim(data)[1]==0){
    type='Others'
  }
  else{
    type=data[data$x==max(data$x),]$cell_type
  }
  cluster_anno=rbind(cluster_anno,c(cl,type))
}

## Color coding by cell names
new.cluster.ids=cluster_anno[,2]
names(new.cluster.ids) <- levels(pbmc)
#pbmc <- RenameIdents(pbmc, new.cluster.ids)

pbmc@meta.data$cell_type=factor(new.cluster.ids[pbmc@meta.data$seurat_clusters])

DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = 'cell_type')

###
ct=NULL
for (l in levels(pbmc)){
  ct=c(ct,pbmc@meta.data[pbmc@meta.data$seurat_clusters==l,]$cell_type[1])
}
names(ct)=levels(pbmc)
pbmc <- RenameIdents(pbmc, ct)

pbmc.markers$gene=as.character(pbmc.markers$gene)
####################################################
## Genes of interest analysis
# input: gene_lsit

gene1=read.table('Aging_genes.txt')
gene2=read.table('Inflammation_genes.txt')
gene_list=c(as.character(gene1$V1),as.character(gene2$V1))

require(gridExtra)
project='B-aggre8'
group=c(1,2,3)

DotPlot_2=function(gene){
df=NULL
for (n in group){
  s=subset(pbmc,subset = sample==n)
  pn=DotPlot(s,features = gene,group.by = 'cell_type')
  
  dfn=pn$data
  dfn$sample=n
  
  df=rbind(df,dfn)
}

df$sample=as.factor(df$sample)

plot=ggplot(df,aes(x=sample,y = id))+
  geom_point(aes(color=avg.exp,size=pct.exp))+
  ggtitle(gene)+
  ylab('Cell Type')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

return (plot)
}


for (i in 1:20){
  outImage=paste0(project,'-dotPlot-',i,'.pdf')
  pdf(outImage,height = 10)
  p1=DotPlot(pbmc,features = gene_list[(i*10-9):(i*10)],group.by = 'cell_type',split.by='sample',cols=colors())
  plot=ggplot(p1$data,aes(x=features.plot,y = id))+
    geom_point(aes(color=avg.exp,size=pct.exp))+
    scale_colour_gradient(low='green',high='red')+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  print(plot)

  dev.off()
}


for (n in group){
  s=subset(pbmc,subset = sample==n)

  for (i in 1:10){
    outImage=paste0(project,'-sample', n,'-vlnPlot-',i,'.pdf')
    pdf(outImage,width = 20,height=18)
    print(VlnPlot(s, features = gene_list[(i*20-19):(i*20)],group.by = 'cell_type', ncol=5))
    dev.off()
    
    outImage=paste0(project,'-sample', n,'-tSNEPlot-',i,'.pdf')
    pdf(outImage,width = 20,height=14)
    print(FeaturePlot(s, features = gene_list[(i*20-19):(i*20)], ncol=5))
    dev.off()
  }
  
}


DotPlot(pbmc, features = toupper(gene_list[1:10]))+ylab('Cluster')
VlnPlot(pbmc, features = toupper(gene_list[1:10]),log = F, ncol=5)
FeaturePlot(pbmc, features = toupper(gene_list[1:10]),ncol=5)

#### DE genes
DE_genes <- FindAllMarkers(pbmc, only.pos = F, min.pct = 0.1, logfc.threshold = 0.2)

# count up and down
up=DE_genes[DE_genes$avg_logFC>0,]
down=DE_genes[DE_genes$avg_logFC<0,]

tb=table(up$gene,up$cluster)>0
tb=as.data.frame(ifelse(tb==T,1,0))
tb$rowSum=rowSums(tb)

label=NULL
for (i in 1:dim(tb)[1]){
  label=c(label,paste(tb[i,1:dim(tb)[2]],collapse=''))
}
tb$label=label

### For genes only apprear in 1 group
tb2=tb[tb$rowSum==1,]

cell=NULL
tag=NULL
for (i in 1:(dim(tb2)[2])){
  cell=c(cell,rep(colnames(tb2)[i],dim(tb2)[1]))
  tag=c(tag,tb2[,i])
}
df=data.frame(cell=cell,tag=tag,label=rep(tb2$label,(length(tb2)-2)),
              rowSum=rep(tb2$rowSum,(length(tb2)-2)),
              gene=rep(rownames(tb2),(length(tb2)-2)))
df$tag=factor(df$tag)
df$Count=''
for (l in unique(df$label)){
  df$Count[which(df$label==l & df$tag==1)[1]]=table(df$label,df$tag)[,2][l]
}

ggplot(df,aes(fill=tag,x=cell,y=gene,group=label))+
  geom_bar(stat="identity")+
  #geom_text(aes(label=Count),vjust=-1)+
  scale_fill_manual(values=c('grey',"#4169E1"))+
  ylab('Number of DE genes')+
  xlab('Cell type')+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")


##### For genes appear in more than 1 group
tb3=tb[tb$rowSum!=1,]
tb3=tb3[order(tb3$rowSum),]
tb3=tb3[order(tb3$label),]

cell=NULL
tag=NULL
for (i in 1:(length(tb3)-2)){
  cell=c(cell,rep(colnames(tb3)[i],dim(tb3)[1]))
  tag=c(tag,tb3[,i])
}
df=data.frame(cell=cell,tag=tag,label=rep(tb3$label,(length(tb3)-2)),
              rowSum=rep(tb3$rowSum,(length(tb3)-2)),
              gene=rep(rownames(tb3),(length(tb3)-2)))
df$tag=factor(df$tag)

ggplot(df,aes(fill=tag,x=cell,y=label,group=label))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c('grey',"#FF6666"))+
  ylab(paste0('Number of DE genes (',dim(df)[1],')'))+
  xlab('Cell type')+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")

################################################
library(SingleR)
#obtain built-in reference data from the Human Primary Cell Atlas
hpca.se <- HumanPrimaryCellAtlasData()
pred.hesc <- SingleR(test = hESCs, ref = hpca.se, labels = hpca.se$label.main)
table(pred.hesc$labels)

dir="Y:\\BI\\Lanqing\\BioPharma\\Single Cell\\20018-01\\aggregated\\filtered_feature_bc_matrix"
singler = CreateSinglerSeuratObject(counts=dir,project.name = 'aggre',
                                    min.genes = 500, species = "Human", 
                                    normalize.gene.length = T, min.cells = 2, npca = 10,
                                    regress.out = "nUMI", reduce.seurat.object = T)

saveRDS(singler, file="Y:\\BI\\Lanqing\\BioPharma\\Single Cell\\20018-01\\aggregated\\singler.rds")

## Seperate groups
barcode=names(singler$other)
group = sapply(barcode, function(x) strsplit(x, '-')[[1]][2])
group = gsub('1','S03',group)
group = gsub('2','S04',group)                                                                                                

# singler$singler[[1]] is the annotations obtained by using ImmGen dataset as reference.
# singler$singler[[2]] is based on the Blueprint+ENCODE reference.
# use singler$singler[[i]]$about for meta-data on the reference.


out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,
                       singler$meta.data$xy, do.label = FALSE, do.letters = F,
                       labels=singler$meta.data$orig.ident,label.size = 4,
                       dot.size = 3)

out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,
                       singler$meta.data$xy,do.label=FALSE,
                       do.letters = F,labels = singler$singler[[1]]$SingleR.single$labels,
                       label.size = 4, dot.size = 3)

out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single.main,
                       singler$meta.data$xy, do.label = FALSE,
                       do.letters = F, dot.size = 2)

out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,
                       singler$meta.data$xy,do.label=F,
                       do.letters =F,labels=singler$other,
                       dot.size = 1.3)

out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,
                       singler$meta.data$xy,do.label=F,
                       do.letters =F,labels=group,
                       dot.size = 1.3)

out$p


kable(table(singler$singler[[1]]$SingleR.single$labels,singler$meta.data$orig.ident))

df = data.frame(x=singler$meta.data$xy[,1],
                y=singler$meta.data$xy[,2],
                t(as.matrix(singler$seurat@data[c('CD3E','CD4','CD8A',
                                                  'CCR7','GZMA','GNLY','MS4A1','CD14','CD34'),])))
df = melt(df,id.vars = c('x','y'))
ggplot(df,aes(x=x,y=y,color=value)) +
  geom_point(size=0.3)+scale_color_gradient(low="gray", high="blue") +
  facet_wrap(~variable,ncol=3) +theme_classic()+xlab('')+ylab('')


pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE)

# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)
top_markers = pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


###########################
library(SCENIC)
## database for mouse
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")

for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
}

load('Y:\\01-New Y Drive\\Bioinformatics\\Lanqing\\BioPharma\\Single Cell\\toyData\\sceMouseBrain.RData')
  
download.file("http://loom.linnarssonlab.org/clone/Previously%20Published/Cortex.loom", "Cortex.loom")
loomPath <- "Y:\\01-New Y Drive\\Bioinformatics\\Lanqing\\BioPharma\\Single Cell\\toyData\\Cortex.loom"

library(SCopeLoomR)
loom <- open_loom(loomPath, mode="r")
exprMat <- get_dgem(loom)
cellInfo <- get_cellAnnotation(loom)
close_loom(loom)

# or into a SingleCellExperiment object:
library(SingleCellExperiment)
sceMouseBrain <- load_as_sce(loomPath)
exprMat <- counts(sceMouseBrain)

cellInfo <- colData(sceMouseBrain)


###############################
library(monocle3)

cds <- load_cellranger_data(dir)

# Load the data
expression_matrix <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_expression.rds"))
cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_colData.rds"))
gene_annotation <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_rowData.rds"))

cell_metadata=pbmc@meta.data
expression_matrix=pbmc[["RNA"]]@counts
gene_annotation=data.frame(gene_short_name=rownames(pbmc[["RNA"]]@counts))
rownames(gene_annotation)=rownames(pbmc[["RNA"]]@counts)


# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# Pre-process the data
cds <- preprocess_cds(cds, num_dim = 10, use_genes = NULL)
cds <- align_cds(cds, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")

# Reduce dimensionality and visualize the results
cds <- reduce_dimension(cds,reduction_method='tSNE')
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
#cds <- order_cells(cds)
#plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "pseudotime")
plot_cells(cds, label_groups_by_cluster=T,  label_cell_groups=F,color_cells_by = "cell_type",
           label_leaves=FALSE,label_branch_points=FALSE,graph_label_size = 10)

# Cluster cells
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by="cao_cell_type")

cds <- reduce_dimension(cds, reduction_method="tSNE")
plot_cells(cds, reduction_method="tSNE", color_cells_by="cao_cell_type")

#### 
cds <- cluster_cells(cds,reduction_method = "UMAP")
cds <- learn_graph(cds)
cds <- order_cells(cds,reduction_method = "UMAP")
plot_cells(cds,color_cells_by = "pseudotime")


############ slingshot
# Normalize 
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
# Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE)

pbmc = RunTSNE(pbmc)
DimPlot(pbmc, reduction = "tsne", label = TRUE)

#######
library(slingshot)
library(scales)
library(ggplot2)

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

sds <- slingshot(Embeddings(pbmc, "tsne"), clusterLabels = pbmc$seurat_clusters)

sds <- slingshot(Embeddings(pbmc, "tsne"), clusterLabels = pbmc@meta.data$cell_type)

#cell_colors <- cell_pal(pbmc$cell_type, brewer_pal("qual", "Set2"))

colorGroup=pbmc$seurat_clusters
colorGroup=factor(pbmc@meta.data$cell_type)

cell_colors_clust <- cell_pal(colorGroup, hue_pal())
color=NULL
for (c in levels(colorGroup)){
  color=c(color,cell_colors_clust[c][1])
}

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

plot(reducedDim(sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(sds, lwd = 2, type = 'lineages', col = 'black')
legend('bottomright', inset=c(-0.3,-0), legend=levels(colorGroup),
       col=color, pch=19,cex=0.8)



curve=slingCurves(sds)[[2]]$s[slingCurves(sds)[[2]]$ord, ]
ggplot(data=as.data.frame(reducedDim(sds)),aes(x=tSNE_1,y=tSNE_2))+
  geom_point(aes(color=factor(pbmc@meta.data$cell_type)))+
  theme_minimal()+
  theme(legend.title = element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())

##################

file='/usr/local/lib/R/site-library/fgsea/extdata/boneMarrow.geneMarkers.celltype.gmt'
f=scan(file,character(),sep='\n')
Cell_Type=NULL
m=NULL
for (line in f){
  s=strsplit(line,'\t')[[1]]
  Cell_Type=c(Cell_Type,s[1])
  m=c(m,paste(s[3:length(s)], collapse = ', '))
  
}

Markers=data.frame(Cell_Type=Cell_Type,Markers=m)
write.table(Markers,file = 'boneMarrow.geneMarkers.celltype.txt',sep='\t',row.names = F, quote = F)

