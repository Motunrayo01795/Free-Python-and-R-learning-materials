#Package installation and library set up

#PACKAGES INSTALLATION
install.packages(c(
  "tidyverse",       # For data manipulation
  "matrixStats",     # For efficient computation of matrix statistics
  "cowplot",         # For combining multiple plots
  "tibble",          # Enhanced data frame structure
  "ggplot2",         # For creating visualizations
  "plotly",          # For creating plots
  "RColorBrewer",    # For generating color palettes for plots
  "gplots",          # For making plots
  "gameofthrones",   # For Game of Thrones-inspired color palettes
  "d3heatmap",       # For creating interactive heatmaps
  "gprofiler2",      # For functional enrichment analysis with g:Profiler
  "readxl"           # For reading Excel files
))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")  # For managing Bioconductor packages

install.packages("gplots")
BiocManager::install(c(
  "edgeR",           # For differential gene expression (DGE) analysis
  "limma",           # For linear modeling and DGE analysis
  "biomaRt",         # For accessing biological databases like Ensembl
  "org.Mm.eg.db",    # Annotation database for Mus musculus (mouse)
  "AnnotationDbi",   # For querying annotation databases
  "clusterProfiler", # For functional enrichment and pathway analysis
  "enrichplot"       # For visualizing enrichment analysis results
))


#import the packages
library(edgeR)        
library(limma)        
library(biomaRt)         
library(org.Mm.eg.db)  
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)  
library(tidyverse)
library(matrixStats)
library(cowplot)
library(tibble)
library(ggplot2)  
library(tidyverse)
library(plotly)
library(RColorBrewer)
library(gplots)
library(d3heatmap)
library(gprofiler2)
library(readxl)
#Set up your working directory
setwd("C:\\Users\\ADMIN\\Documents\\MyThirdRProject")
#Load your data
my_data <- read.csv("GSE261050_gene_counts_all.csv")

#Data preprocessing 
#Checking for missing values 
summary(my_data)
sum(is.na(my_data))
colSums(is.na(my_data))
str(my_data)     # Data types and structure
dim(my_data)     # Dimensions (rows, columns)
names(my_data)   # Column names

#Subsetting dataset
data1 <- my_data[ ,c("X","AN13999.BA32","AN06706.BA32","AN17402.BA32","AN05797.BA32",
                    "AN07303.BA32","AN14576.BA32","AN10747.BA32","AN11050.BA32",
                    "AN18901.BA32","AN18056.BA32","AN01733.BA32","AN04446.BA32",
                    "AN06400.BA32","AN02737.BA32","AN08855.BA32","AN08675.BA32",
                    "AN12979.BA32","AN02671.BA32", "AN07243.BA32","AN00331.BA32",
                    "AN05086.BA32")]


rownames(data1) <- data1$X
#Setting up the experimental design 

groups <- factor(c("AD","AD","AD","AD","AD","AD","AD","AD","AD","AD",
                   "Control","Control", "Control", "Control","Control","Control",
                   "Control","Control","Control","Control","Control"))
target <- data.frame(row.names = colnames(data1), groups)
sample_labels <- rownames(target)

#Gene Ontolology

# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")  # For human genes

# Map Ensembl Gene IDs to Gene Names
genes <- data1$X # Example Ensembl Gene IDs

# Query biomaRt to retrieve gene names
gene_names <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id",
                    values = genes,
                    mart = ensembl)
colnames(annotated_data)

annotated_data <- merge(data1, gene_names, by.x ="X", 
                   by.y = "ensembl_gene_id"  )
new_order <- c("X", "external_gene_name", "AN13999.BA32", "AN06706.BA32", "AN17402.BA32",
"AN05797.BA32","AN07303.BA32","AN14576.BA32", "AN10747.BA32", "AN11050.BA32", "AN18901.BA32",
"AN18056.BA32", "AN01733.BA32", "AN04446.BA32", "AN06400.BA32", "AN02737.BA32", 
"AN08855.BA32", "AN08675.BA32", "AN12979.BA32", "AN02671.BA32", "AN07243.BA32", 
"AN00331.BA32", "AN05086.BA32")
annotated_data <- annotated_data[new_order]
annotated_data <- annotated_data[, -1]
colnames(annotated_data)[1] <- "genes"

#Saved selected and annotated data into your PC 
write.csv(annotated_data, "Annotated_data.csv")

#Normalization and transformation 

#Getting the DEG list 
my_deglist <- DGEList(annotated_data)
#filtering to remove the lowly expressed genes
cpm <- cpm(my_deglist, log = TRUE)
cpm.df <- as.data.frame(cpm)

#Seting threshold for lowly expressed genes
fil_threshold <- rowSums(cpm>1)>=2 #Eplain 
my_deglist.filtered <- cpm.df[fil_threshold,]


#Normalize filtered data 

my_deglist.filtered.Norm <- calcNormFactors(my_deglist.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(my_deglist.filtered.Norm, log = TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
log

#Hierarchical clustering

distance <- dist(t(log2.cpm.filtered.norm), method = "euclidean")
head(distance)
png("Hierachical Clustering of Samples.png",width = 400, height = 500)
distance <- dist(t(log2.cpm.filtered.norm), method = "euclidean")
clusters <- hclust(distance,method = "complete")
plot(clusters, labels = sample_labels)
dev.off()

#Principal Component analysis 
pca.res <- prcomp(t(log2.cpm.filtered.norm),scale.= F, retx = T)

#eigenValues
pc.var <- pca.res$sdev^2
pc.percentage <- round(pc.var/sum(pc.var )*100, 1)
pca.res.df <- as_tibble(pca.res$x)

#Plotting Principal Component Analysis
png("PCA.png", width= 400, height = 500)
ggplot(pca.res.df, aes(x=PC1,y=PC2, label = sample_labels,colour = groups)) +
  geom_point(size= 4) + geom_label() +
  xlab(paste0("PC1(",pc.percentage[1],"%",")"))+
  ylab(paste0("PC2(",pc.percentage[2],"%",")"))+
  labs(title = "Principal Component Analysis") +
  theme_bw()
dev.off()

png("PCA2.png", width= 400, height = 500)
ggplot(pca.res.df, aes(x=PC3,y=PC4, label = sample_labels,colour = groups)) +
  geom_point(size= 4) + geom_label() +
  xlab(paste0("PC3(",pc.percentage[3],"%",")"))+
  ylab(paste0("PC4(",pc.percentage[4],"%",")"))+
  labs(title = "Principal Component Analysis") +
  theme_bw()
dev.off()

png("PCA3.png", width= 400, height = 500)
ggplot(pca.res.df, aes(x=PC5,y=PC6, label = sample_labels,colour = groups)) +
  geom_point(size= 4) + geom_label() +
  xlab(paste0("PC5(",pc.percentage[5],"%",")"))+
  ylab(paste0("PC6(",pc.percentage[6],"%",")"))+
  labs(title = "Principal Component Analysis") +
  theme_bw()
dev.off()
group <- factor(target$groups)
design <-model.matrix(~0 + group)
colnames(design) <- levels(group)


