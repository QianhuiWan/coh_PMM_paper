
# Note: use test2 output as filnal output, i.e. lightskyblue

library(tidyverse)
library(biomaRt)
library(WGCNA)
library(DESeq2)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 8)

#===============================================================================
#
#  prepare input files
#
#==============================================================================

# count matrix and sample metadata needed


#===============================================================================
#
#  Read the gene counts table and plot the sample tree
#
#===============================================================================
setwd("~/githubRepo/coh_PMM/RNAseq_scripts/gene_counts_stat/WGCNA")

# Read the gene counts table 
data0 <- read.table("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/gene_count_matrix_STARcounts_PMM.tsv", header=T,row.names=1,sep="\t")

PMM_metadata <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/PMM_sampleSheet.tsv")
L1_tetx_num <- read_tsv(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/activationScore_R2/filtered_L1/num_active_all_L1s_L1PA2_L1PA3_pmm.tsv")
L1_TEtx_num <- L1_tetx_num %>% 
  dplyr::select(sample = sampleID, L1, L1PA2, L1PA3) %>% unique() %>% 
  `colnames<-`(c("sample", "L1_as", "L1PA2_as", "L1PA3_as"))

sample_metadata <- PMM_metadata %>% 
  `colnames<-`(c("sample_ID", "group", "L1PA2_st", "L1PA3_st")) %>% 
  left_join(L1_TEtx_num, by = c("sample_ID" = "sample")) 

# Normalization with lcpm
lcpm_PMM <- readr::read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/lcpm_TMM_PMM_STARcounts.tsv")

dataExpr_lcpm <- lcpm_PMM %>% 
  pivot_wider(names_from = sample, 
              values_from = lcpm) %>% 
  tibble::column_to_rownames(var = "gene_name")

lcpm_matrix = as.matrix(dataExpr_lcpm)

### Remove rows with all NAs
# lcpm_matrix <- lcpm_matrix[rowSums(is.na(lcpm_matrix)) != ncol(lcpm_matrix), ]

### Remove genes with very low average FPKM
# keep_expr <- rowMeans(lcpm_matrix, na.rm = TRUE) > 0
keep_expr <- rowSums(lcpm_matrix>0, na.rm = TRUE) > 3
lcpm_matrix_filt <- lcpm_matrix[keep_expr, ] #26925

### Calculate variance per gene
# gene_var <- apply(fpkm_matrix_filt, 1, var)

### Keep top N most variable genes
# top_n <- 32000  # or 5000 for faster performance
# top_genes <- names(sort(gene_var, decreasing = TRUE)[1:top_n])

# fpkm_matrix_filt_final <- fpkm_matrix_filt[top_genes, ]
lcpm_matrix_filt_final <- lcpm_matrix_filt

datExpr = t(log2(lcpm_matrix_filt_final+1))

head(datExpr[1:5,1:5]) # samples in row, genes in column
match(sample_metadata$sample_ID, colnames(data0))
# datExpr <- datExpr[,1:5000]



# Calculate sample distance and cluster the samples
sampleTree = hclust(dist(datExpr), method = "average");
# plot sample tree
pdf(file = "1-n-sampleClustering_test2.pdf", width = 40, height = 9);
par(cex = 1.3);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()


#===============================================================================
#
#  Choose soft threshold parameter
#
#===============================================================================

# Choose a set of soft threshold parameters
powers = c(c(1:20), seq(from = 22, to=30, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, 
                        networkType = "signed",
                        # networkType = "signed hybrid",
                        # networkType = "unsigned",
                        verbose = 5)
# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file = "2-n-sft_signned_test2.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#===============================================================================
#
#  Turn data expression into topological overlap matrix
#
#===============================================================================

# Turn data expression into topological overlap matrix
power=sft$powerEstimate #12
# power=12

# Option 1: automatic
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = power, 
                       TOMType = "signed",
                       networkType = "signed",
                       # networkType = "signed hybrid",
                       # TOMType = "unsigned",  # still works well here
                       minModuleSize = 30,
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       verbose = 3)
cor<- stats::cor
# unsigned -> nodes with positive & negative correlation are treated equally 
# signed -> nodes with negative correlation are considered *unconnected*, treated as zero

sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
pdf(file = "12-module_tree_blockwise_test2.pdf", width = 8, height = 6);
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

saveRDS(net, file = "WGCNA_test2_net.rds")
################################################################################
################################################################################
# Option 2a: step-by-step
power = power
adjacency = adjacency(datExpr, power = power)
TOM = TOMsimilarity(adjacency); # Turn adjacency into topological overlap
dissTOM = 1-TOM


# Option 2b: 
TOM = TOMsimilarityFromExpr(datExpr, power = power)
dissTOM = 1-TOM 
dim(dissTOM)


#===============================================================================
#
#  Construct modules (proceed with the gene tree from option 2b)
#
#===============================================================================
# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file = "gene_cluster_signned_test2.pdf", width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()

# Module identification using dynamic tree cut
# We like large modules, so we set the minimum module size relatively high:
# minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
                            pamRespectsDendro = FALSE, minClusterSize = 30);
table(dynamicMods)
length(table(dynamicMods)) 
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf(file = "module_tree_signned_test2.pdf", width = 8, height = 6);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, 
                    guideHang = 0.05,main = "Gene dendrogram and module colors")
dev.off()

#===============================================================================
#
#  Merge modules
#
#===============================================================================
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Merge close modules
MEDissThres=0.40
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
mergedColors = merge$colors  
mergedMEs = merge$newMEs  
# Plot merged module tree
pdf(file = "merged_Module_Tree_signned_test2.pdf", width = 12, height = 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
dev.off()
#write.table(merge$oldMEs,file="oldMEs.txt");
#write.table(merge$newMEs,file="newMEs.txt");

#===============================================================================
#
#  Export of networks to external software
#
#===============================================================================

# Export the gene list of old modules 
for (i in 1:length(merge$oldMEs)){
  modules = c(substring(names(merge$oldMEs)[i], 3));
  genes = colnames(datExpr)
  inModule = is.finite(match(dynamicColors,modules))
  modGenes = genes[inModule]
  modTOM=TOM[inModule,inModule]
  dimnames(modTOM)=list(modGenes,modGenes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("output_for_cytoscape/orign_CytoscapeInput-edges-", paste(modules, collapse="-"), "_test2.txt", sep=""),
                                 nodeFile = paste("output_for_cytoscape/orign_CytoscapeInput-nodes-", paste(modules, collapse="-"), "_test2.txt", sep=""),
                                 weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = dynamicColors[inModule]);
}


# Export the gene lists of merged modules to Cytoscape format
modules <- unique(mergedColors)

for (mod in modules) {
  inModule <- mergedColors %in% mod
  modGenes <- colnames(datExpr)[inModule]
  
  # Subset TOM
  modTOM <- TOM[inModule, inModule]
  dimnames(modTOM) <- list(modGenes, modGenes)
  
  # Export to Cytoscape
  exportNetworkToCytoscape(
    modTOM,
    edgeFile = paste0("output_for_cytoscape/merge_CytoscapeInput-edges-", mod, "_test2.txt"),
    nodeFile = paste0("output_for_cytoscape/merge_CytoscapeInput-nodes-", mod, "_test2.txt"),
    weighted = TRUE,
    threshold = -1,  # Export all edges
    nodeNames = modGenes,
    nodeAttr = mergedColors[inModule]
  )
}



# find ZNF141 and ZNF382 module
# 1. Extract gene symbols from rownames (assuming format like ZNF141_ENSG000001...)
gene_symbols <- sub("_.*", "", colnames(datExpr))
# 2. List of genes you care about
genes_of_interest <- c("ZNF141", "ZNF382")
# 3. Find matching indices
match_idx <- match(genes_of_interest, gene_symbols)
# 4. Check their assigned modules
dynamicColors[match_idx]

gene2module <- data.frame(
  gene_symbol = gene_symbols,
  full_gene_id = colnames(datExpr),
  module = dynamicColors,
  stringsAsFactors = FALSE
)
subset(gene2module, gene_symbol %in% genes_of_interest)


gene2mergedModule <- data.frame(
  gene_symbol = gene_symbols,
  full_gene_id = colnames(datExpr),
  merged_module = mergedColors,
  stringsAsFactors = FALSE
)
subset(gene2mergedModule, gene_symbol %in% genes_of_interest)
#lightskyblue


#===============================================================================
#
#  PART 1: Correlate module eigen-genes and samples (or other discrete data)
#
#===============================================================================
# Heatmap of old module eigen-genes and samples
pdf(file="oldMEs_signned.pdf")
library("pheatmap")
rownames(merge$oldMEs)=names(data0)
pheatmap(merge$oldMEs,cluster_col=T,cluster_row=T,show_rownames=F,show_colnames=T,fontsize=6)
dev.off()


# Heatmap of new module eigen-genes and sample trait (e.g. Zone)
col_ann <- sample_metadata[,c(1,2)] %>% as.data.frame()
rownames(col_ann) <- col_ann[,1]
col_ann <- data.frame(col_ann)
col_ann$group <- as.factor(col_ann$group)
col_ann <- col_ann[order(col_ann$group),]
col_ann$sample_ID <- NULL
head(col_ann)
ann_color <- list("col_ann" = c("BMPC" = "green",
                                "PMM" = "red"))

data <- data.frame(merge$newMEs)
data <- data[order(match(rownames(data), rownames(col_ann))),]
dim(merge$newMEs)

pdf(file="newMEs_signned.pdf")
rownames(merge$newMEs)=names(data0)
pheatmap(data,cluster_col=T,cluster_row=F,show_rownames=F,
         show_colnames=T,fontsize=12,
         annotation_row = col_ann, annotation_colors = ann_color)
dev.off()

#=====================================================================================
#
#  PART 2: Correlation between gene modules and microbial traits (continuous data)
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs = orderMEs(MEs0)

# Read microbial data as traits
# bac_traits = read.table("b_order_234.txt", header = T, sep = "\t")
bac_traits <- sample_metadata[,c(1,6)] %>% as.data.frame()
rownames(bac_traits) = bac_traits[, 1]
bac_traits$sample_ID <- NULL
# sample names should be consistent in eigen genes and traits !!!!
# bac_traits = bac_traits[match(rownames(MEs), rownames(bac_traits)), ]
table(rownames(MEs) == rownames(bac_traits))

# Calculate pearson correlation coefficients between module eigen-genes and traits
moduleTraitCor = cor(MEs, bac_traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#write.table(moduleTraitCor,file="moduleTrait_correlation.txt");
#write.table(moduleTraitPvalue,file="moduleTrait_pValue.txt");


#  Plot heatmap of module-traits relationship
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf("5-module-traits-order_signned.pdf", width = 80, height = 15)
par(mar = c(15, 12, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(bac_traits),
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#  Plot heatmap of module-traits relationship
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf("5-module-traits-order1_signned.pdf", width = 20, height = 10)
par(mar = c(15, 12, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(bac_traits),
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#=====================================================================================
#
#   Intramodular analysis: identifying genes with high geneModuleMembership & geneTraitSignficance
#
#=====================================================================================


# Define variable Verru containing the Verrucomicrobiales column of bac_traits
Verru = as.data.frame(bac_traits$L1PA2_as);
names(Verru) = "L1PA2"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

MET = orderMEs(cbind(MEs, Verru))

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, Verru, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Verru), sep="");
names(GSPvalue) = paste("p.GS.", names(Verru), sep="");

# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(10,10,1,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


module = "blue4"
# Rename to moduleColors
moduleColors = mergedColors
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for L1PA2",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "blue")

## Draw bubble plot for particular module
colsum_bac_traits <- colSums(bac_traits)
colsum_bac_traits <- data.frame(colsum_bac_traits)
colsum_bac_traits$b_order <- rownames(colsum_bac_traits)
library(tidyr)
moduleTraitCor_long <- data.frame(moduleTraitCor)
moduleTraitCor_long$module <- rownames(moduleTraitCor)
moduleTraitCor_long <- moduleTraitCor_long
# moduleTraitCor_long <- gather(moduleTraitCor_long, L1PA2_as, factor_key = TRUE)
moduleTraitCor_long <- gather(moduleTraitCor_long, "b_order", "PCC", L1PA2_as, factor_key = TRUE)


moduleTraitPvalue_long <- data.frame(moduleTraitPvalue)
moduleTraitPvalue_long$module <- rownames(moduleTraitPvalue)
moduleTraitPvalue_long <- moduleTraitPvalue_long
moduleTraitPvalue_long <- gather(moduleTraitPvalue_long, "b_order", "pval", L1PA2_as, factor_key = TRUE)

moduleTrait_long <- merge(moduleTraitCor_long, moduleTraitPvalue_long, by = c("module","b_order"))

bubble_Data <- merge(moduleTrait_long, colsum_bac_traits, by = "b_order")
#just want module = "lightgreen"
bubble_Data_lightgreen <- bubble_Data[which(bubble_Data$module == "MEblue4"),]

library(ggplot2)
ggplot(bubble_Data_lightgreen, aes(x= colsum_bac_traits, y= PCC, size = colsum_bac_traits,
                                   color = PCC, label = b_order)) +
  geom_text(hjust = 1, size=3) +
  geom_point(alpha=1) + ylab("Module-taxon correlation") + xlab("Relative abundance (sum)") +
  theme_bw()

############# Summary ###################################

head(datExpr)[1:5,1:5] # transcriptome data

head(sample_metadata)[1:5,] # metadata (sample info)
# head(bac_traits)[1:5,1:5] # external trait


#=====================================================================================
#
#   Cytoscape
#
#=====================================================================================


#if(!"RCy3" %in% installed.packages()){
#  install.packages("BiocManager")
#  BiocManager::install("RCy3")
#}

# https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/
library(RCy3)

cytoscapePing () # make sure cytoscape is open
cytoscapeVersionInfo ()

###### for yellow module of the merged data (newMEs) #################################
edge <- read.delim("output_for_cytoscape/merge_CytoscapeInput-edges-lightgreen.txt")
colnames(edge)
colnames(edge) <- c("source", "target","weight","direction","fromAltName","toAltName")

node <- read.delim("output_for_cytoscape/merge_CytoscapeInput-nodes-lightgreen.txt")
colnames(node)  
colnames(node) <- c("id","altName","node_attributes") 

createNetworkFromDataFrames(node,edge[1:50,], title="my first network", collection="DataFrame Example")

################ customise the network visualization ##################################
# use other pre-set visual style
setVisualStyle('Marquee')

# set up my own style
style.name = "myStyle"
defaults <- list(NODE_SHAPE="diamond",
                 NODE_SIZE=30,
                 EDGE_TRANSPARENCY=120,
                 NODE_LABEL_POSITION="W,E,c,0.00,0.00")
nodeLabels <- mapVisualProperty('node label','id','p')
nodeFills <- mapVisualProperty('node fill color','node_attributes','d',c("A","B"), c("#FF9900","#66AAAA"))
arrowShapes <- mapVisualProperty('Edge Target Arrow Shape','interaction','d',c("activates","inhibits","interacts"),c("Arrow","T","None"))
edgeWidth <- mapVisualProperty('edge width','weight','p')

createVisualStyle(style.name, defaults, list(nodeLabels,nodeFills,arrowShapes,edgeWidth))
setVisualStyle(style.name)




