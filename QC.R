library(Seurat)
library(cowplot)
library(stringr)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)
library(gridExtra)
library(rtracklayer)
library(GenomicRanges)
saved.seed <- 6789



gff_path <- "/media/robb-e532/mr.twelve/Downloads/ncbi_dataset_Genes(Annotation Features)/ncbi_dataset/data/GCF_015227675.2/genomic.gff"
gff_data <- import.gff(gff_path)

mitochondrial_genes <- subset(gff_data, seqnames == "NC_001665.2")

mitochondrial_genes_df <- as.data.frame(mitochondrial_genes) # converting the object into DF

head(mitochondrial_genes_df)
unique(seqnames(gff_data))

#storing all the mitochondrial genes in a df.
mito_gene_df <- mitochondrial_genes_df[mitochondrial_genes_df$type == "gene", ]

head(mito_gene_df)
#Extracting only gene column
mito_gene_df <- mito_gene_df[12]

write.table(mito_gene_df, file = "mitochondrial_genes.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)



path_General <- "/media/robb-e532/mr.twelve/Cellranger Count/"
path_KMY  <- paste0(path_General,"Sample_SRR11143096/outs/filtered_feature_bc_matrix/") %>% print()
path_KMO  <- paste0(path_General, "Sample_SRR11143097/outs/filtered_feature_bc_matrix/")   %>% print()
path_KMCR <- paste0(path_General, "Sample_SRR11143098/outs/filtered_feature_bc_matrix/")  %>% print()
path_KFY  <- paste0(path_General, "Aggr_99-102/outs/count/filtered_feature_bc_matrix/")   %>% print()
path_KFO  <- paste0(path_General, "Aggr_103-106/outs/count/filtered_feature_bc_matrix/")   %>% print()
path_KFCR <- paste0(path_General, "Aggr_107-110/outs/count/filtered_feature_bc_matrix/")  %>% print()

KMY <- CreateSeuratObject(Read10X(data.dir = path_KMY),
                          min.cells = 3,        #Include features detected in at least this many cells
                          min.features = 200)   #Include cells where at least this many features are detected
KMY$Sample <- "Young.M.AL"
KMY$Sex <- "Male"
KMY$Age <- "Young"
KMY$Diet <- "Ad.Libitum"

KMO <- CreateSeuratObject(Read10X(data.dir = path_KMO),
                          min.cells = 3,        #Include features detected in at least this many cells
                          min.features = 200)   #Include cells where at least this many features are detected
KMO$Sample <- "Old.M.AL"
KMO$Sex <- "Male"
KMO$Age <- "Old"
KMO$Diet <- "Ad.Libitum"

KMCR <- CreateSeuratObject(Read10X(data.dir = path_KMCR),
                           min.cells = 3,        #Include features detected in at least this many cells
                           min.features = 200)   #Include cells where at least this many features are detected
KMCR$Sample <- "Old.M.CR"
KMCR$Sex <- "Male"
KMCR$Age <- "Old"
KMCR$Diet <- "Restricted"

KFY <- CreateSeuratObject(Read10X(data.dir = path_KFY),
                          min.cells = 3,        #Include features detected in at least this many cells
                          min.features = 200)   #Include cells where at least this many features are detected
KFY$Sample <- "Young.F.AL"
KFY$Sex <- "Female"
KFY$Age <- "Young"
KFY$Diet <- "Ad.Libitum"

KFO <- CreateSeuratObject(Read10X(data.dir = path_KFO),
                          min.cells = 3,        #Include features detected in at least this many cells
                          min.features = 200)   #Include cells where at least this many features are detected
KFO$Sample <- "Old.F.AL"
KFO$Sex <- "Female"
KFO$Age <- "Old"
KFO$Diet <- "Ad.Libitum"

KFCR <- CreateSeuratObject(Read10X(data.dir = path_KFCR),
                           min.cells = 3,        #Include features detected in at least this many cells
                           min.features = 200)   #Include cells where at least this many features are detected
KFCR$Sample <- "Old.F.CR"
KFCR$Sex <- "Female"
KFCR$Age <- "Old"
KFCR$Diet <- "Restricted"




all.merged <- merge(KMY, c(KMO,KMCR,KFY,KFO,KFCR))

dim(all.merged)
table(all.merged$Sample, all.merged$orig.ident)


rm(path_KFCR, path_KFO, path_KFY, path_KMCR, path_KMO, path_KMY, gff_path, mitochondrial_genes_df, mitochondrial_genes)
rm(KMY,KMO,KMCR,KFY,KFO,KFCR, gff_data)


all.merged.copy <- all.merged

# Extract the vector of mitochondrial gene names from the data frame
mito_genes <- mito_gene_df$Name

#Getting orignal gene names as well as the index stored.
original_gene_names <- rownames(all.merged.copy@assays$RNA)

# Find mitochondrial genes in the object
mito_genes_in_all.merged <- mito_genes[mito_genes %in% original_gene_names]

# Create a new gene names vector
new_gene_names <- original_gene_names

# Renaming the mitochondrial genes by adding "MT-" prefix
indexes_of_mito_genes_in_all.merged <- match(mito_genes_in_all.merged, original_gene_names)
new_gene_names[indexes_of_mito_genes_in_all.merged] <- paste0("MT-", mito_genes_in_all.merged)

# Update the original gene names with the new names for mitochondrial genes
original_gene_names[indexes_of_mito_genes_in_all.merged] <- new_gene_names[indexes_of_mito_genes_in_all.merged]

# Now that original_gene_names has everything same, except the mitochondrial genes renamed, push it back to allmerged.copy
rownames(all.merged.copy@assays$RNA) <- original_gene_names
all.merged <- all.merged.copy

write.table(rownames(all.merged.copy@assays$RNA), file = "$Final.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

rm(indexes_of_mito_genes_in_all.merged, mito_genes, mito_genes_in_all.merged, new_gene_names, original_gene_names)


#Percentage feature set for Mitochondrial Genes

grep("^MT-",
     rownames(all.merged@assays$RNA$counts),
     value = T)

all.merged <- PercentageFeatureSet(all.merged, "^MT-",
                                   col.name = "percent.MT")
head(all.merged$percent.MT)


#Percentage feature set for and Ribosomal genes

grep("^Rp[sl]",
     rownames(all.merged@assays$RNA$counts.4),
     value = T)

all.merged <- PercentageFeatureSet(all.merged, "^Rp[sl]",
                                   col.name = "percent.Ribosomal")
head(all.merged$percent.Ribosomal)

# Percentage largest gene
#Run apply() over the columns (cells) and calculate what percentage of the data comes from the single most observed gene. Exclude Malat1 as it is usually highly expressed and affects results
# We will get:
# The count for the largest gene per cell
# The index position of the gene with the largest count
# The name of the most highly expressed gene per cell


all.merged[rownames(all.merged) != "Malat1",] -> SObject.F

apply(
  SObject.F@assays$RNA$counts.1,
  2,
  max) -> SObject.F$largest_count
apply(
  SObject.F@assays$RNA$counts.2,
  2,
  max) -> SObject.F$largest_count
apply(
  SObject.F@assays$RNA$counts.3,
  2,
  max) -> SObject.F$largest_count
apply(
  SObject.F@assays$RNA$counts.4,
  2,
  max) -> SObject.F$largest_count
apply(
  SObject.F@assays$RNA$counts.5,
  2,
  max) -> SObject.F$largest_count

head(SObject.F$largest_count)
apply(
  SObject.F@assays$RNA$counts.6,
  2,
  max) -> SObject.F$largest_count

apply(
  SObject.F@assays$RNA$counts.1,
  2,
  which.max) -> SObject.F$largest_index
head(SObject.F$largest_index)
tail(SObject.F$largest_index)
apply(
  SObject.F@assays$RNA$counts.2,
  2,
  which.max) -> SObject.F$largest_index
apply(
  SObject.F@assays$RNA$counts.3,
  2,
  which.max) -> SObject.F$largest_index
apply(
  SObject.F@assays$RNA$counts.4,
  2,
  which.max) -> SObject.F$largest_index
apply(
  SObject.F@assays$RNA$counts.5,
  2,
  which.max) -> SObject.F$largest_index
apply(
  SObject.F@assays$RNA$counts.6,
  2,
  which.max) -> SObject.F$largest_index

rownames(SObject.F)[SObject.F$largest_index] -> SObject.F$largest_gene
100 * SObject.F$largest_count / SObject.F$nCount_RNA -> SObject.F$percent.Largest.Gene

SObject.F$largest_gene -> all.merged$largest_gene

SObject.F$percent.Largest.Gene -> all.merged$percent.Largest.Gene
rm(SObject.F)



VlnPlot(all.merged,
        features=c("nFeature_RNA",
                   "nCount_RNA",
                   "percent.MT",
                   "percent.Ribosomal",
                   "percent.Largest.Gene"),
        ncol = 3)
VlnPlot(all.merged,
        features=c("percent.Largest.Gene"),
        ncol = 1) +
  scale_y_log10()

VlnPlot(all.merged,
        features=c("nFeature_RNA",
                   "nCount_RNA",
                   "percent.MT",
                   "percent.Ribosomal",
                   "percent.Largest.Gene"),
        group.by = "Sample",
        ncol = 3)
VlnPlot(all.merged,
        features=c("percent.Largest.Gene"),
        group.by = "Sample",
        ncol = 1) +
  scale_y_log10()





grid.arrange(
  FeatureScatter(all.merged,
                 feature1 = "nCount_RNA",
                 feature2 = "percent.Largest.Gene",
                 group.by = "orig.ident",
                 pt.size = 0.5),
  
  FeatureScatter(all.merged,
                 "nCount_RNA",
                 "nFeature_RNA",
                 group.by = "orig.ident",
                 pt.size = 0.5),
  ncol = 2)


## QC plots using ggplot2

qc.metrics <- as_tibble( all.merged[[]], rownames="Cell.Barcode")
head(qc.metrics)

## nFeatures vs. nCounts

qc.metrics %>%
  arrange(percent.MT) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.MT)) +
  geom_point() +
  scale_color_gradientn(colors=c("black" , "blue" , "green2" , "red" , "yellow")) +
  ggtitle("QC metrics") +
  geom_hline(yintercept = 500) +
  geom_hline(yintercept = 5500)


#Yintercept 5500 changed cause more cells are present over 5000 than the old matrices



qc.metrics %>%
  arrange(percent.MT) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.MT)) +
  geom_point(size=0.5) +
  scale_color_gradientn(colors=c("black" , "blue" , "green2" , "red" , "yellow")) +
  ggtitle("QC metrics") +
  geom_hline(yintercept = 500) +
  geom_hline(yintercept = 5500) +
  scale_y_continuous() +
  scale_x_log10()


## Complexity: log10(genes)/log10(counts)
...or a regression line through the cloud and then calculate the difference from the observed value to the expected.



qc.metrics %>%
  mutate(complexity=log10(nFeature_RNA) / log10(nCount_RNA)) -> qc.metrics

lm(log10(qc.metrics$nFeature_RNA)~log10(qc.metrics$nCount_RNA)) -> complexity.lm

qc.metrics %>%
  mutate(
    complexity_diff = log10(nFeature_RNA) - ((log10(qc.metrics$nCount_RNA)*complexity.lm$coefficients[2])+complexity.lm$coefficients[1])
  ) -> qc.metrics



qc.metrics %>%
  ggplot(aes(x=complexity_diff)) +
  geom_density(fill="yellow")

min(c(max(qc.metrics$complexity_diff),0-min(qc.metrics$complexity_diff))) -> complexity_scale


qc.metrics %>%
  mutate(complexity_diff=replace(complexity_diff,complexity_diff< -0.1,-0.1)) %>%
  ggplot(aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=complexity_diff)) +
  geom_point(size=0.5) +
  geom_abline(slope=complexity.lm$coefficients[2], intercept = complexity.lm$coefficients[1]) +
  scale_colour_gradient2(low="blue2",mid="grey",high="red2")

### Note: 10^2.75 = 562.... this may be a value that make sense for subset(nFeatures_RNA > 560)... Also, nCounts_RNA < 30000 (10^4.5)?
# For me, I believe 2.75 as lower bound bound for nFeatureRNA > 560 but upper bound 3.75(10^3.75 = 5600), nFeature_RNA < 5600 ~ 5500 at the above graph.
#Also for nCount_RNA upper bound maybe 4.5, 10^4.5 = 31600 ~ 30000.

# Largest Gene

qc.metrics %>%
  ggplot(aes(x=complexity_diff,
             y=percent.Largest.Gene)) +
  geom_point(size=0.5)



qc.metrics %>%
  group_by(largest_gene) %>%
  count() %>%
  arrange(desc(n)) -> largest_gene_list
largest_gene_list





largest_gene_list %>%
  filter(n>1000)%>%
  pull(largest_gene) -> largest_genes_to_plot
#1000 cause of histogram


qc.metrics %>%
  filter(largest_gene %in% largest_genes_to_plot) %>%
  mutate(largest_gene=factor(largest_gene, levels=largest_genes_to_plot)) %>%
  arrange(largest_gene) %>%
  ggplot(aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=largest_gene)) +
  geom_point(size=1) +
  scale_colour_manual(values=c("grey", RColorBrewer::brewer.pal(9,"Set1")))




qc.metrics %>%
  filter(largest_gene %in% largest_genes_to_plot) %>%
  mutate(largest_gene=factor(largest_gene, levels=largest_genes_to_plot)) %>%
  arrange(largest_gene) %>%
  ggplot(aes(x=complexity_diff,
             y=percent.Largest.Gene,
             colour=largest_gene)) +
  geom_point(size=1) +
  scale_colour_manual(values=c("grey", RColorBrewer::brewer.pal(9,"Set1")))



qc.metrics %>%
  arrange(percent.MT) %>%
  ggplot(aes(x=complexity_diff, y=percent.Largest.Gene, colour=percent.MT)) +
  geom_point() +
  scale_colour_gradient(low="grey", high="red2")



qc.metrics %>%
  arrange(percent.Ribosomal) %>%
  ggplot(aes(x=complexity_diff,
             y=percent.Largest.Gene,
             colour = percent.Ribosomal)) +
  geom_point() +
  scale_colour_gradient(low="grey", high="red2")






qc.metrics %>%
  ggplot(aes(percent.MT)) +
  geom_histogram(binwidth =0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage Mitochondrion") +
  geom_vline(xintercept =13)



qc.metrics %>%
  ggplot(aes(percent.Largest.Gene)) +
  geom_histogram(binwidth =0.7, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage Largest Gene") +
  geom_vline(xintercept =10)






###############################################
###############################################
# Filter Seurat object

all.merged

all.filtered <- subset(all.merged,
                       nFeature_RNA > 560 &
                         nFeature_RNA < 5500 &
                         nCount_RNA < 30000 &
                         percent.MT < 40 &
                         percent.Ribosomal < 25 &
                         percent.Largest.Gene < 25)
all.filtered


# Plot QC metrics in filtered object

VlnPlot(all.filtered,
        features=c("nFeature_RNA",
                   "nCount_RNA",
                   "percent.MT",
                   "percent.Ribosomal",
                   "percent.Largest.Gene"),
        ncol = 3)
VlnPlot(all.filtered,
        features=c("percent.Largest.Gene"),
        ncol = 1) +
  scale_y_log10()



VlnPlot(all.filtered,
        features=c("nFeature_RNA",
                   "nCount_RNA",
                   "percent.MT",
                   "percent.Ribosomal",
                   "percent.Largest.Gene"),
        group.by = "Sample",
        ncol = 3)
VlnPlot(all.filtered,
        features=c("percent.Largest.Gene"),
        group.by = "Sample",
        ncol = 1) +
  scale_y_log10()



grid.arrange(
  FeatureScatter(all.filtered,
                 feature1 = "nCount_RNA",
                 feature2 = "percent.Largest.Gene",
                 group.by = "Sample",
                 pt.size = 1),
  
  FeatureScatter(all.filtered,
                 "nCount_RNA",
                 "nFeature_RNA",
                 group.by = "Sample",
                 pt.size = 1),
  
  ncol=2)



VlnPlot(all.merged,
        features=c("nFeature_RNA",
                   "nCount_RNA",
                   "percent.MT",
                   "percent.Ribosomal",
                   "percent.Largest.Gene"),
        ncol = 5)
VlnPlot(all.filtered,
        features=c("nFeature_RNA",
                   "nCount_RNA",
                   "percent.MT",
                   "percent.Ribosomal",
                   "percent.Largest.Gene"),
        ncol = 5)
```

```{r, fig.width=12, fig.height=8}
grid.arrange(
  FeatureScatter(all.merged,
                 feature1 = "nCount_RNA",
                 feature2 = "percent.Largest.Gene",
                 group.by = "orig.ident",
                 pt.size = 0.5),
  FeatureScatter(all.filtered,
                 feature1 = "nCount_RNA",
                 feature2 = "percent.Largest.Gene",
                 group.by = "orig.ident",
                 pt.size = 0.5),
  FeatureScatter(all.merged,
                 "nCount_RNA",
                 "nFeature_RNA",
                 group.by = "orig.ident",
                 pt.size = 0.5),
  FeatureScatter(all.filtered,
                 "nCount_RNA",
                 "nFeature_RNA",
                 group.by = "orig.ident",
                 pt.size = 0.5),
  ncol=2)




all.filtered <- NormalizeData(all.filtered,
                              normalization.method ="LogNormalize",
                              verbose = T)



gene.expression <- apply(all.filtered@assays$RNA$data,
                         1, mean)


gene.expression <- sort(gene.expression,
                        decreasing = T)

layers <- list(
  all.filtered@assays$RNA$data.1,
  all.filtered@assays$RNA$data.2,
  all.filtered@assays$RNA$data.3,
  all.filtered@assays$RNA$data.4,
  all.filtered@assays$RNA$data.5,
  all.filtered@assays$RNA$data.6
)

combined_data <- matrix(0, nrow = nrow(layers[[1]]), ncol = ncol(layers[[1]]))

for (layer in layers) {
  combined_data <- combined_data + matrix(as.numeric(layer), nrow = nrow(layer))
}

head(gene.expression, n=25)



ggplot(mapping = aes(all.filtered@assays$RNA$data["Gapdh",])) +
  geom_histogram(binwidth =0.05, fill="yellow", colour="black") +
  ggtitle("Gapdh expression")
```
#####################
#####################

# Cell cycle scoring genes

# Seurat cell cycle stage marker genes
cc.genes.updated.2019
# Non-human format:
str_to_title(cc.genes.updated.2019)

## Predict cell cycle stage
all.filtered <- CellCycleScoring(all.filtered,
                                 s.features = str_to_title(cc.genes.updated.2019$s.genes),
                                 g2m.features = str_to_title(cc.genes.updated.2019$g2m.genes),
                                 set.ident = T)

all.filtered[[]]
```
```{r}
as_tibble(all.filtered[[]]) %>%
  ggplot(aes(Phase,
             fill = Sample)) +
  geom_bar( )
as_tibble(all.filtered[[]]) %>%
  ggplot(aes(Sample,
             fill = Phase)) +
  geom_bar( )
```
## Plot
```{r, fig.width=5,fig.height=4}
as_tibble(all.filtered[[]]) %>%
  ggplot(aes(x = S.Score,
             y = G2M.Score,
             color=Phase)) +
  geom_point(size = 1)
#coord_cartesian(xlim=c(-0.3,0.3),
#               ylim=c(-0.3,0.3))
```
## Filter high scoring cells
```{r}
all.filtered

all.filtered <- subset(all.filtered,
                       S.Score < 0.2 &
                         G2M.Score < 0.2)
all.filtered
```

#####################
#####################
# Find variable features
```{r}
all.filtered <- FindVariableFeatures(all.filtered,
                                     selection.method ="vst",
                                     nfeatures = 500)
head(VariableFeatures(all.filtered), 15)
```
## Top 15 most highly variable genes
```{r}
top15 <- head(VariableFeatures(all.filtered), 15)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(all.filtered)
plot2 <- LabelPoints(plot = plot1,
                     points = top15,
                     repel = T)
plot2
rm(top15,plot1, plot2)
```

```{r}
variance.data <- as_tibble(HVFInfo(all.filtered),
                           rownames ="Gene")

variance.data <- variance.data %>%
  mutate(hypervariable = Gene %in% VariableFeatures(all.filtered))

head(variance.data, n=10)
```

```{r, fig.width=6,fig.height=4}
variance.data %>%
  ggplot(aes(log(mean),log(variance),color=hypervariable)) +
  geom_point(size = 0.5) +
  scale_color_manual(values=c("black","red"))
```

##########################
##########################

# Scaling
Before dimensionality reduction genes will be scaled to have a mean of 0 and a variance of 1.
Here is the time to regress unwanted variables e.g.  vars.to.regress = c("S.Score", "G2M.Score")

```{r}
all.filtered <- ScaleData(all.filtered)
#vars.to.regress = c("S.Score", "G2M.Score","Age", "Diet", "Sex"))
```

###############################################################################
###############################################################################
# Dimensionality Reduction
# PCA
```{r}
all.filtered <- RunPCA(all.filtered,
                       #features = VariableFeatures(all.filtered),
                       nfeatures.print = 10)
```
```{r}
ElbowPlot(all.filtered,
          ndims = 50,
          reduction = "pca")
```


```{r, fig.width=15, fig.height=15}
grid.arrange(DimPlot(all.filtered,
                     reduction = "pca"),
             DimPlot(all.filtered,
                     reduction = "pca",
                     group.by = "Sex"),
             DimPlot(all.filtered,
                     reduction = "pca",
                     group.by = "Diet"),
             DimPlot(all.filtered,
                     reduction = "pca",
                     group.by = "Age"),
             
             ncol(2), nrow(2))
```


```{r, fig.height=3, fig.width=8}
DimHeatmap(all.filtered,
           dims=1:3, cells=10000)
```

# tSNE
The parameter 'perplexity value' is the expected number of nearest neighbors. By default this is set to 30.
Lower values will help resolve small clusters, at the expense of larger clusters becoming more diffuse. On the contrary, higher values will make the larger clusters more distinct, but may lose smaller clusters.

```{r}
all.filtered <- RunTSNE(all.filtered,
                        dims = 1:40,
                        seed.use = saved.seed,
                        perplexity = 30)
```

```{r, fig.width=15, fig.height=15}
grid.arrange(DimPlot(all.filtered,
                     reduction = "tsne"),
             DimPlot(all.filtered,
                     reduction = "tsne",
                     group.by = "Sex"),
             DimPlot(all.filtered,
                     reduction = "tsne",
                     group.by = "Diet"),
             DimPlot(all.filtered,
                     reduction = "tsne",
                     group.by = "Age"),
             ncol(2))
```


# UMAP
```{r}
all.filtered = RunUMAP(all.filtered,
                       dims = 1:40,
                       verbose = T)
```

```{r, fig.width=15, fig.height=15}
grid.arrange(DimPlot(all.filtered,
                     reduction = "umap"),
             DimPlot(all.filtered,
                     reduction = "umap",
                     group.by = "Sex"),
             DimPlot(all.filtered,
                     reduction = "umap",
                     group.by = "Diet"),
             DimPlot(all.filtered,
                     reduction = "umap",
                     group.by = "Age"),
             ncol(2))
```

#####################
#####################

# Find doblets

## Opytimize parameters with paramSweep
```{r}
sweep.list <- paramSweep_v3(all.filtered,
                            PCs = 1:40,
                            sct = F,
                            num.cores = 1)
```
```{r}
sweep.stats <- summarizeSweep(sweep.list,
                              GT = F)
Find.pK.stats <- find.pK(sweep.stats)
```
```{r}
barplot(Find.pK.stats$BCmetric, names.arg = Find.pK.stats$pK, las=2)
```
# Set up papameters
```{r}
nPCs <- 1:30                                 # Relevant principal components
nExp <- round(ncol(all.filtered) * 0.05)     # expect 5% doublets
pK <- 0.27                                   # x for max BCmetric in previous plot

nExp
pK
```
# Run Doublet finder
```{r}
all.filtered <- doubletFinder_v3(all.filtered,
                                 pN = 0.25,
                                 PCs = nPCs,    # Number of artificial doublets. Set by default
                                 pK = pK,       # Neighborhood size to compute artificial nearest neighbors
                                 nExp = nExp)   # Number of expected doublets
```

```{r, fig.width=20, fig.height=10}
# name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(all.filtered@meta.data)[grepl("DF.classification", colnames(all.filtered@meta.data))]
```


```{r, fig.width=10, fig.height=10}
grid.arrange(DimPlot(all.filtered,
                     reduction = "umap",
                     group.by = "orig.ident"),
             DimPlot(all.filtered,
                     reduction = "umap",
                     group.by = DF.name),
             ncol(2))
```

```{r}
VlnPlot(all.filtered,
        features = "nFeature_RNA",
        group.by = DF.name,
        pt.size = 0.1)
```

```{r}
dim(all.filtered)
all.filtered = all.filtered[, all.filtered@meta.data[, DF.name] == "Singlet"]
dim(all.filtered)





