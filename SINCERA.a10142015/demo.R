# demo.R
# Author: Minzhe Guo (minzhe.guo@cchmc.org)
# Description: 
#				A demonstration of using SINCERA to analyse E16.5 mouse lung single cell data.
# 				The data are RNA-seq of 148 single cells from two independent sample preparations 
#               (86 cells from sample 1 and 62 cells from sample 2) from E16.5 mouse lung and RNA
#               expression values were calculated using the FPKM.
#				Through the pipeline analysis, we distinguished major cell types of fetal mouse lung, 
#               including epithelial, endothelial, smooth muscle, pericyte, and fibroblast-like cell types, 
#               and identified cell type specific gene signatures, bioprocesses, and key regulators.
# Dependencies: 
#             SINCERA a10142015
#             E16.5.Rda
# Date: October 14, 2015

# ============================= 

# TODO Add expression matrix to hierarchical clustering
# TODO Remove useless genes (low variance, not expressed)

#######################################################
#             Load Sincera and SC3 Code      #
#######################################################

# library(SingleCellExperiment)
# library(SC3)
# library(scater)
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("ALL")

#NOTE: please make sure that demo.R, sincera.R, and E16.5.Rda are placed in the same folder and the folder is the current working directory of R

source("sincera.R")
load("E16.5.Rda")

head(cells)

head(genes)

expressions[1:5,1:5]

#expressions <- read.table(file = 'JD81.tsv', sep = '\t', header = TRUE)
#expressions <- as.matrix(read.table('JD81.tsv', header=TRUE, sep="\t", row.names=1, as.is=TRUE))
#cells <- read.table(file = 'cells.tsv', sep = '\t', header = TRUE)
#genes <- read.table(file = 'genes.tsv', sep = '\t', header = TRUE)

#expressions <- read.table(file = 'JD81.tsv', header=TRUE, sep="\t", row.names=1, as.is=TRUE)
#cells <- read.table(file = 'cells.tsv', header=TRUE, sep="\t", row.names=1, as.is=TRUE)
#genes <- read.table(file = 'genes.tsv', header=TRUE, sep="\t", row.names=1, as.is=TRUE)

# expressions_n <- read.table(file = 'magic_output.tsv', header=TRUE, sep="\t", row.names=1, as.is=TRUE)
# cells_n <- read.table(file = 'magic_cells.tsv', header=TRUE, sep="\t", row.names=1, as.is=TRUE)
# genes_n <- read.table(file = 'magic_genes.tsv', header=TRUE, sep="\t", row.names=1, as.is=TRUE)

expressions_n <- read.table(file = 'magic_output_k_1.tsv', header=TRUE, sep="\t", row.names=1, as.is=TRUE)
cells_n <- read.table(file = 'magic_k_cells.tsv', header=TRUE, sep="\t", row.names=1, as.is=TRUE)
genes_n <- read.table(file = 'magic_k_genes.tsv', header=TRUE, sep="\t", row.names=1, as.is=TRUE)

#######################################################
#      Organized the Data using ExpressionSet         #
#######################################################

# cells is a data frame containing the cell information; rownames will be used as the identifiers of cells
pd <- new("AnnotatedDataFrame", data = cells_n)

# genes is a data frame containing the gene information; rownames will be used as the identifiers of genes; a "SYMBOL" column is required for encoding the official symbols of genes
fd <- new("AnnotatedDataFrame", data = genes_n)

# expressions is a 2D data matrix encoding the epxression values
# rownames of expressions should be exactly the same as the rownames of genes
# colnames of expressions should be exactly the same as the rownames of cells
ES <- new("ExpressionSet", exprs = as.matrix(expressions_n), phenoData = pd, featureData = fd)

ES

cat("\nE16.5 Demonstration\n\n")

#######################################################
#                Gene Pre-filtering                   #
#######################################################

# truncate expression values out of bounds
# ES <- exprs.truncate(ES, upper.bound=NULL, lower.bound=0.0001)

# determine the threshold for specificity filter based on the expression of a set of reference genes, e.g., ribosomal genes
# the function first identifies ubiquitously expressed ref genes: reference genes with expression >= exp.t in at least exp.p percentage of cells per sample
# measure the specificity of ubiquitously expressed ref genes
# finding a threshold between 0 and 1 that filters out at least (1-tolerance) percentage of ubiquitously expressed ref genes
#mouse.ribosomal.genes
#ref.idx <- which(fData(ES)$SYMBOL %in% mouse.ribosomal.genes)
#specificity.threshold=specificity.criterion.selection(ES, group.by="SAMPLE", groups=NULL, ref.idx=ref.idx, exp.t=5, exp.p = 0.95, tolerance=0.05, step=0.01)
#specificity.threshold

# calculate per-sample specificity and abundancy scores for each profile
#ES <- prefiltering(ES, group.by="SAMPLE", groups=NULL, specificity.threshold=specificity.threshold, abundancy.threshold=5, abundancy.count=2, mode="intersect", fig.res=150, do.batch.analysis=TRUE, export=T)

# select profiles that passed prefiltering filters for further analysis
#ES <- ES[which(fData(ES)[,PREFILTERING.PREFIX]==1),]


#######################################################
#                Cluster Assignment                   #
#######################################################

# # perform gene-by-gene per-sample z-score transformation
# ESz <- normalization.zscore(ES, group.by="SAMPLE", groups=NULL, verbose=TRUE) 
# 
# # hierarchical clustering
# ret <- cluster.assignment(ESz, cluster.label="CLUSTER", 
#                        clustering.method="hc", 
#                        distance.method="pearson", #"spearman","euclidean" 
#                        linkage.method="average", 
#                        num.singleton=0, h=0.5, k=NULL, 
#                        do.shift=TRUE, do.plot=TRUE, 
#                        verbose=TRUE, export=TRUE, export.components="pd")
# 
# 
# # consensus clustering
# #notrun  ret <- cluster.assignment(ESz, cluster.label="consensus_cluster", clustering.method="consensus", distance.method="pearson", linkage.method="average", min.area.increase=0.2, maxK=10, reps=10, pItem=0.8, pFeature=1, clusterAlg="hc",verbose=TRUE, export=TRUE, export.components="pd")
# 
# # tight clustering
# #notrun  ret <- cluster.assignment(ESz, cluster.label="tight_cluster", clustering.method="tight", target=1, k.min=25, verbose=TRUE, export=TRUE, export.components="pd")
#    
# 
# # the cluster membership shall be encoded in ret$cell.cluster
# cell_cluster <- ret$cell.cluster
# 
# # if hc is used, cell order in the dendrogram shall be encoded in ret$cell.order
# cell_order <- ret$cell.order
# 
# # load cell type markers from extenal files
# # tow columns are required: SYMBOL-official symbol of markers, TYPE-cell type
# # example:   SYMBOL TYPE
# #            Epcam  Epithelial
# #            Cdh1   Epithelial
# #            Pecam1 Endothelial
# #            Emcn   Endothelial
# markers <- read.table(file="markers.txt", sep="\t", head=T)
# 
# head(markers)
# 
# # log normalization: increase expression by 1 and log2 normalized
# ES.log2 <- normalization.log(ES, log.base=2, increase=1)
# 
# # order the cells in the order as they appeared in the dendrogram if hc was used for clustering
# ES.log2 <- cluster.ordering(ES.log2, col.order=cell_order)
# ESz <- cluster.ordering(ESz, col.order=cell_order)
# 
# # export clustering results for visualization in GENE-E
# # export log2 expression of all genes
# expression.export4viz(ES.log2, file.prefix="allgenes.log2.")
# # export zscore expression of cell type markers
# expression.export4viz(ESz, signature=markers, cell.attrs=c("SAMPLE", "CLUSTER"), file.prefix="markers.zscore.")
# 
# # print("1")
# 
# ES.sqrt <- normalization.sqrt(ES)
# ES.sqrt <- cluster.ordering(ES, col.order=cell_order)
# # plot the squqre root normalized gene expression across cells ordered by clusters
# #plotProfiles(ES.sqrt, genes=c("Foxm1","Tbx3", "Acta2","Actg2", "Pdgfrb", "Fn1","Tcf21", "Emcn","Kdr","Fcer1g", "Fcgr2b", "Cdh1", "Epcam"), fig.filename="marker_profiles.tif", ylabel="fpkm, sqrt")
# 
# # print("2")
# 
# # permutation analysis to determine significance of the derived clustering scheme
# # n - the number of permutations
# # for quick demonstration
# cluster.permutation.analysis(ES, group.by=CLUSTER.LABEL, n=20, distance.method="euclidean", log.base=2, verbose=T) 
# 
# # for reproducing the permutation results in the paper, may take more than one hour 
# if (FALSE) {
#     set.seed(10)
#     cluster.permutation.analysis(ES, group.by="CLUSTER", n=5000, distance.method="euclidean", log.base=2, verbose=T) 
# }

#######################################################
#                Differential Expression              #
#######################################################

welch.groups <- c()
wilcoxon.groups <- c()

clusters <- sort(unique(pData(ES)$CLUSTER))
for (i in clusters) {
    # if both the numbers of cells in the cluster and outside the cluster are greater than 5,
    # use welch's t-test, otherwise, use wilcoxon rank sum test
    if (length(which(pData(ES)$CLUSTER %in% i))>5 & length(which(!(pData(ES)$CLUSTER %in% i)))>5) {
        welch.groups <- c(welch.groups, i)
    } else {
        wilcoxon.groups <- c(wilcoxon.groups, i)
    }
}

welch.groups
wilcoxon.groups

# per cluster differential expression test
if (length(welch.groups)>0) {
    ES <- diff.test(ES, group.by="CLUSTER", groups=welch.groups, method="welch", diffexpr.prefix="test_", do.fdr=FALSE, verbose=T, export=T, export.components="fd")
}
if (length(wilcoxon.groups)>0) {
    ES <- diff.test(ES, group.by="CLUSTER", groups=wilcoxon.groups, method="wilcoxon", diffexpr.prefix="test_", do.fdr=FALSE, verbose=T, export=T, export.components="fd")
}

# differential expression using SAMseq
# ES <- diff.test(ES, group.by="CLUSTER", groups=NULL, method="samseq", diffexpr.prefix="samseq_", do.fdr=FALSE, verbose=T, export=T, export.components="fd")

# retrieve differentially expressed genes
diffgenes <- get.diff.genes(ES, group.by="CLUSTER", groups=NULL, diffexpr.prefix="test_", threshold=0.01) 

# export expression of differentially expressed genes for visualization
dg <- diffgenes[,c("SYMBOL", "CLUSTER")]
colnames(dg) <- c("SYMBOL", "TYPE")
expression.export4viz(ESz, signature=dg, cell.attrs=c("SAMPLE", "CLUSTER"), file.prefix="diffgenes.zscore.")

print("differential enrichment done")
#######################################################
#                Cell Type Enrichment                 #
#######################################################

# selecting cluster specific differentially expressed genes for cell type enrichment analysis
groups <- sort(unique(pData(ES)$CLUSTER))
prefix="use_for_celltype_"
for (i in groups) {
	i.de.name <- paste(prefix, i, sep="")
	i.test.name <- paste(DIFF.EXPR.PREFIX, i, sep="")
	fData(ES)[,i.de.name] <- 0
	if (i %in% c("3")) {
		fData(ES)[which(fData(ES)[,i.test.name] < 0.03), i.de.name] <- 1 
	} else {
		fData(ES)[which(fData(ES)[,i.test.name] < 0.01), i.de.name] <- 1 
	}
}

# initialize knowledge base for cell type enrichment analysis

genome <- NULL
associations <- NULL
KB <- NULL

genome <- rownames(genes)
associations <- associations.01112014


#' Contruct the knowledge base for cell type enrichment analysis
#' 
#' associations (data.frame) a data frame containing the gene and cell type associations
#' genome (character) the full set of genes in a scRNA-seq data
#' species (character) the species of the genome: MUSMU - mouse, HOMSA - human
#' id.type (character) the type of ids of the genes in the genome: ENSEMBL - ENSEMBL gene id, SYMBOL -  Entrez Gene Symbol, EG - Entrez Gene Id
#' return a list of 3 items: celltype.gene.association - genome related associations, celltype.genome.count - genome-wide cell type associations, genome - symbol mapping for the genome
KB <- celltype.enrichment.initKB(associations, genome, species="MUSMU", id.type="ENSEMBL")


#' Contruct the knowledge base for cell type enrichment analysis
#' 
#' ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' KB (list) the knowledge base prepared by celltype.enrichment.initKB()
#' group.by (character) the name of the column that contains the cluster information
#' groups (character) the clusters for cell type enrichment
#' species (character) the species of the genome: MUSMU - mouse, HOMSA - human
#' id.type (character) the type of ids of the genes in the genome: ENSEMBL - ENSEMBL gene id, SYMBOL -  Entrez Gene Symbol, EG - Entrez Gene Id
#' celltype.enrichment.prefix (character) the prefix of columns encoding the cluster-specific gene list for cell type enrichment analysis
celltype.enrichment(ES, KB, group.by="CLUSTER", groups=NULL, species="MUSMU", id.type="ENSEMBL", celltype.enrichment.prefix=CELLTYPE.ENRICHMENT.PREFIX, verbose=TRUE)



#######################################################
#                Marker-based Validation              #
#######################################################

ESz <- normalization.zscore(ES, group.by=SAMPLE.LABEL, groups=NULL, verbose=TRUE)
 
#' Cell type validation using the expression patterns of known biomarkers
#' 
#' ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' group.by (character) the name of the column that contains the cluster information
#' groups (character) the clusters for cell type validation
#' threshold (numeric) the threshold of expression; 0 means that markers will not be used to rank cells where their expression < 0
#' marker.prefix (character) the prefix of columns encoding the cell type markers
celltype.validation(ESz, group.by=CLUSTER.LABEL, groups=c(1,2,3,5,7,8,9), marker.prefix=MARKER.PREFIX, threshold=0)


#######################################################
#                 Signature Prediction                #
#######################################################

#' Cell type specific signature prediction
#' 
#' ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' group.by (character) the name of the column that contains the cluster information
#' groups (character) the clusters for signature prediction
#' trainset.prefix (character) prefix for labeling the columns encoding the training instances for each group
#' train.class.prefix (character) prefix for labeling the columns encoding the class of group specific training instances
#' marker.prefix (character) prefix for labeling the columns encoding the biomarkers for each group
#' testset.prefix (character) prefix for labeling the columns encoding the testing instances for each group
#' common.prefix (character) prefix for labeling the columns encoding the results of common metric for each group
#' common.threshold (numeric) in common gene metric, expression > COMMON.TRESHOLD will considered as expressed
#' common.percentage (numeric) in common gene metric, genes that express in >COMMON.PRECENTAGE cluster cells will be considered as a common
#' unique.prefix (character) prefix for labeling the columns encoding the results of unique metric for each group
#' unique.ratio (numeric) parameter for unique gene metric
#' unique.quantile (numeric) parameters for unique gene metric
#' test.statistic.metric.prefix (character) prefix for labeling the columns encoding the results of test statistic metric for each group 
#' diff.expr.prefix (character) prefix for labeling the columns encoding the results of differential test
#' diff.expr.threshold (numeric) genes with differential expression p-value<0.05 will be selected as candidates for signature 
#' syn.sim.prefix (character) prefix for labeling the columns encoding the results of synthetic profile similarity metric for each group
#' signature.prefix (character) prefix for labeling the columns encoding the results of signature prediction
#' gene.symbol.label (character) the name of the column encoding gene symbols
#' verbose (logical)
#' export (logical) wheterh to export the ExpressionSet
#' export.components (character) the components of an ExpressionSet object that will be exported: fd- gene information, pd-cell information, expr-expression
#' return an ExpressionSet object containing the results of signature prediction in the attributes of fData
ES <- signature.analysis(ES, group.by = CLUSTER.LABEL, groups=c(1,2,3,5,7,8,9),
                                # training set
                                trainset.prefix = TRAINSET.PREFIX,
                                train.class.prefix = TRAIN.CLASS.PREFIX,
                                marker.prefix=MARKER.PREFIX,
                                # testing set
                                testset.prefix=TESTSET.PREFIX,
                                # common gene metric
                                common.prefix = COMMON.PREFIX,
                                common.threshold=COMMON.TRESHOLD, common.percentage=COMMON.PERCENTAGE,
                                # unique gene metric
                                unique.prefix = UNIQUE.PREFIX,
                                unique.ratio=UNIQUE.RATIO, unique.quantile=UNIQUE.QUANTILE,
                                # test statistic metric
                                test.statistic.metric.prefix = TEST.STATS.METRIC.PREFIX, log.base=2,
                                diff.expr.prefix=DIFF.EXPR.PREFIX,
                                diff.expr.threshold = 0.05,
                                # synthetic profile similarity metric 
                                syn.sim.prefix = SYN.SIM.PREFIX,
                                signature.prefix = SIGNATURE.PREFIX,
                                gene.symbol.label = GENE.SYMBOL.LABEL,
                                verbose=TRUE,
                                export=TRUE, export.components="fd" )




# cross validation of signature prediction for cluster 1,2,3,5,7,8,9

# configuring the columns for validation of signature prediciton
sig.groups <- c(1,2,3,5,7,8,9)

# specifying the number of signature genes for each cluster
sig.n <- rep(0, length(sig.groups))
j <- 1
for (i in sig.groups) {
    if (i==3) { # for cluster 3, selecting the top k as signature genes, where k equals to the number of cluster 3 specific differentially expressed genes with p-value<0.03
        test.i <- which(fData(ES)[, paste("test_", i, sep="")] < 0.03) 
    } else { # for cluster 1,2,5,7,8,9, selecting the top k as signature genes, where k equals to the number of cluster specific differentially expressed genes with p-value<0.01
        test.i <- which(fData(ES)[, paste("test_", i, sep="")] < 0.01)
    }
    sig.n[j] <- length(test.i)
    j <- j + 1
}

# for quick demonstration
# diff.threshold=1 means that for each cluster, we consider the genes with pvalue(diff.test)<1 as the candidates for signature prediction
# percentages=c(0.2) means that for each cluster, we sample 80% as training data, the remaining 20% as testing data
# repeats = 2 means that for each cluster, we repeat subsampling twice
signature.validation(ES, group.by="CLUSTER", groups=sig.groups, signature.prefix = SIGNATURE.PREFIX, diff.threshold=1, sig.n=sig.n, percentages=c(0.2), repeats=2, verbose=FALSE, export=FALSE) 

# for reproducing the cross validation results in the paper; this may take about 2 hours
if (FALSE) {
    set.seed(10)
    signature.validation(ES, group.by="CLUSTER", groups=sig.groups, signature.prefix = SIGNATURE.PREFIX, diff.threshold=1, sig.n=sig.n, percentages=c(0.2), repeats=100, verbose=TRUE, export=FALSE) 
}

#######################################################
#                Driving Force Analysis               #
#######################################################


groups <- c(9) # performing the driving force analysis for the epithelial cells (cluster 9)
prefix.tg=TG.PREFIX # the prefix for labeling columns encoding the potential cluster-specific transcription factors
prefix.tf=TF.PREFIX # the prefix for labeling columns encoding the potential cluster-specific regulatory targets
tf.name = TF.LABEL # the name of the column encoding the transcritpion factor annotation

# selecting cluster specific differentially expressed genes as candidate regulatory targets
# selecting cluster specific differentially expressed TFs or commonly expressed TFs as candidate TFs

option = 1 # configurations for quick demonstration of driving force analysis
# option = 2 # configurations for reproducing the driving force analysis results in the manuscript (platform: Windows 7, R 3.2.0)

if (option == 1) {     # for quick demonstration

	for (i in groups) {
		i.test.name <- paste("test_", i, sep="") # the name of the column encoding the p-values of differential expression 
		i.common.name <- paste("common_", i, sep="") # the name of the column encoding the common gene metric
		i.tg.name <- paste(prefix.tg, i, sep="") 
		i.tf.name <- paste(prefix.tf, i, sep="")
		fData(ES)[,i.tg.name] <- 0
		fData(ES)[,i.tf.name] <- 0
		
		# selecting genes with pvalue(diff.test)<0.0001 as the potential regulatory targets
		fData(ES)[which(fData(ES)[,i.test.name] < 0.0001), i.tg.name] <- 1 
        
		# selecting genes with pvalue(diff.test)<0.005 and with transcription factor annotation as the potential transcription factors
		fData(ES)[which(fData(ES)[,i.test.name] < 0.005), i.tf.name] <- 1 
		fData(ES)[which(fData(ES)[,tf.name] == 0), i.tf.name] <- 0	
	}

} else if (option == 2) { 

    # configurations for reproducing the driving force analysis results in the manuscript
	# this would take about 1 hour for running
    
	for (i in groups) {
		i.test.name <- paste("test_", i, sep="") # the name of the column encoding the p-values of differential expression 
		i.common.name <- paste("common_", i, sep="") # the name of the column encoding the common gene metric
		i.tg.name <- paste(prefix.tg, i, sep="")
		i.tf.name <- paste(prefix.tf, i, sep="")
		fData(ES)[,i.tg.name] <- 0
		fData(ES)[,i.tf.name] <- 0
		
        
		# selecting genes with pvalue(diff.test)<0.01 as the potential regulatory targets
		fData(ES)[which(fData(ES)[,i.test.name] < 0.01), i.tg.name] <- 1
        
		# selecting cluster-specific common genes with pvalue(diff.test)<0.005 and with transcription factor annotation as the potential transcription factors
		fData(ES)[which(fData(ES)[,i.test.name] < 0.05), i.tf.name] <- 1
		fData(ES)[which(fData(ES)[,i.common.name] == 1), i.tf.name] <- 1
		fData(ES)[which(fData(ES)[,tf.name] == 0), i.tf.name] <- 0	
	}
}

# after selecting the potential transcription factors and regulatory targets 
driving.force.analysis(ES, group.by="CLUSTER", groups=c(9), unique.by="SYMBOL")


################################################################################
#      Consensus Maximization for Nkx2-1 Regulatory Targets Prediction         #
################################################################################

consensus_maximization(Nkx2.1_data, prior=c(1,0, 1,0.5, 1,0.5, 1,0.5), epslon=0, max.iter=100)


cat("\n\nDemonstration completed\n\n")
