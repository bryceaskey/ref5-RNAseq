# Differential gene expression analysis using edgeR package
# Goal - perform pairwise comparisons of expression data
library(edgeR)
library(tibble)

# Specify paths to data files
justCounts <- "C:/Users/Bryce/Research/ref5-RNAseq/data/expression-counts/"
geneLengths <- "C:/Users/Bryce/Research/ref5-RNAseq/data/gene-lengths.txt"
araport11 <- read.csv("C:/Users/Bryce/Research/ref5-RNAseq/data/Araport11.csv", header=FALSE)

# Load gene expression and gene length data
countFiles <- paste(justCounts, dir(justCounts), sep="")
countLabels <- c("H10_rep1", "H10_rep2", "H10_rep3", "H55_rep1", "H55_rep2", "H55_rep3", "R5_rep1",
                 "R5_rep2", "R5_rep3", "WT_rep1", "WT_rep2", "WT_rep3", "Y6_rep1", "Y6_rep2", "Y6_rep3")
countGroups <- c("A", "A", "A", "B", "B", "B", "C", "C", "C", "D", "D", "D", "E", "E", "E")
expressionData <- readDGE(countFiles, labels=countLabels, group=countGroups)
expressionData$genes <- read.delim(geneLengths, row.names=1)

# Only keep genes with significant level of expression
keep <- filterByExpr(expressionData, group=countGroups)
expressionData <- expressionData[keep, , keep.lib.sizes=FALSE]

# Normalize for library size
expressionData <- calcNormFactors(expressionData, method="TMM")

# Create MDS plot to verify integrity of RNAseq data
plotMDS(expressionData)

# Perform DEG analysis to compare gene expression of NAA treated plants to control
expDesign <- model.matrix(~0+group, data=expressionData$samples)
expressionData <- estimateDisp(expressionData, design=expDesign)
colnames(expDesign) <- levels(expressionData$samples$group)
modelFit <- glmQLFit(expressionData, expDesign)

contrast <- makeContrasts(A-D, levels=expDesign)
H10_WT <- glmQLFTest(modelFit, contrast=contrast)
H10_DEGs <- topTags(H10_WT, n=nrow(expressionData), p.value=1)$table
H10_DEGs <- rownames_to_column(H10_DEGs, var="locus")

contrast <- makeContrasts(B-D, levels=expDesign)
H55_WT <- glmQLFTest(modelFit, contrast=contrast)
H55_DEGs <- topTags(H55_WT, n=nrow(expressionData), p.value=1)$table
H55_DEGs <- rownames_to_column(H55_DEGs, var="locus")

contrast <- makeContrasts(C-D, levels=expDesign)
R5_WT <- glmQLFTest(modelFit, contrast=contrast)
R5_DEGs <- topTags(R5_WT, n=nrow(expressionData), p.value=1)$table
R5_DEGs <- rownames_to_column(R5_DEGs, var="locus")

contrast <- makeContrasts(E-D, levels=expDesign)
Y6_WT <- glmQLFTest(modelFit, contrast=contrast)
Y6_DEGs <- topTags(Y6_WT, n=nrow(expressionData), p.value=1)$table
Y6_DEGs <- rownames_to_column(Y6_DEGs, var="locus")

# Load gene annotations
araport11 <- araport11[, c(1, 3, 4, 13)]
colnames(araport11) <- c("locus", "short_name", "name", "aliases")

# Add gene annotations to DEG dataframes
H10_DEGs <- merge(araport11, H10_DEGs, by="locus", all.y=TRUE)
H55_DEGs <- merge(araport11, H55_DEGs, by="locus", all.y=TRUE)
R5_DEGs <- merge(araport11, R5_DEGs, by="locus", all.y=TRUE)
Y6_DEGs <- merge(araport11, Y6_DEGs, by="locus", all.y=TRUE)

write.csv(H10_DEGs, file="../Research/ref5-RNAseq/data/DEGs/H10_all_genes.csv", row.names=FALSE)
write.csv(H55_DEGs, file="../Research/ref5-RNAseq/data/DEGs/H55_all_genes.csv", row.names=FALSE)
write.csv(R5_DEGs, file="../Research/ref5-RNAseq/data/DEGs/R5_all_genes.csv", row.names=FALSE)
write.csv(Y6_DEGs, file="../Research/ref5-RNAseq/data/DEGs/Y6_all_genes.csv", row.names=FALSE)


# # Calculate fpkm for all genes
# fpkmAllGenes <- as.data.frame(rpkm(expressionData$counts, gene.length=expressionData$genes$Length, normalized.lib.sizes=TRUE, log=FALSE))
# fpkmAllGenes <- rownames_to_column(fpkmAllGenes, var="locus")
# 
# # Calculate tmm for all genes
# tmmAllGenes <- as.data.frame(cpm(expressionData, log=FALSE))
# tmmAllGenes <- rownames_to_column(tmmAllGenes, var="locus")
# 
# # Add tmm and fpkm to Cyp79A2 DEG dataframe ----
# fpkm_WT <- vector(mode="numeric", length=nrow(Cyp79A2_DEGs))
# fpkmSE_WT <- vector(mode="numeric", length=nrow(Cyp79A2_DEGs))
# tmm_WT <- vector(mode="numeric", length=nrow(Cyp79A2_DEGs))
# tmmSE_WT <- vector(mode="numeric", length=nrow(Cyp79A2_DEGs))
# fpkm_Cyp79A2 <- vector(mode="numeric", length=nrow(Cyp79A2_DEGs))
# fpkmSE_Cyp79A2 <- vector(mode="numeric", length=nrow(Cyp79A2_DEGs))
# tmm_Cyp79A2 <- vector(mode="numeric", length=nrow(Cyp79A2_DEGs))
# tmmSE_Cyp79A2 <- vector(mode="numeric", length=nrow(Cyp79A2_DEGs))
# 
# for(i in 1:nrow(Cyp79A2_DEGs)){
#   geneLocus <- Cyp79A2_DEGs$locus[i]
#   
#   fpkm_WT[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 8:10]))
#   fpkmSE_WT[i] <- sd(c(unlist(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 8:10])))/sqrt(3)
#   
#   tmm_WT[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 8:10]))
#   tmmSE_WT[i] <- sd(c(unlist(tmmAllGenes[tmmAllGenes$locus==geneLocus, 8:10])))/sqrt(3)
#   
#   fpkm_Cyp79A2[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 2:4]))
#   fpkmSE_Cyp79A2[i] <- sd(c(unlist(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 2:4])))/sqrt(3)
#   
#   tmm_Cyp79A2[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 2:4]))
#   tmmSE_Cyp79A2[i] <- sd(c(unlist(tmmAllGenes[tmmAllGenes$locus==geneLocus, 2:4])))/sqrt(3)
# }
# 
# Cyp79A2_DEGs$fpkm_WT <- fpkm_WT
# Cyp79A2_DEGs$fpkmSE_WT <- fpkmSE_WT
# 
# Cyp79A2_DEGs$tmm_WT <- tmm_WT
# Cyp79A2_DEGs$tmmSE_WT <- tmmSE_WT
# 
# Cyp79A2_DEGs$fpkm_Cyp79A2 <- fpkm_Cyp79A2
# Cyp79A2_DEGs$fpkmSE_Cyp79A2 <- fpkmSE_Cyp79A2
# 
# Cyp79A2_DEGs$tmm_Cyp79A2 <- tmm_Cyp79A2
# Cyp79A2_DEGs$tmmSE_Cyp79A2 <- tmmSE_Cyp79A2
# 
# Cyp79A2_DEGs <- Cyp79A2_DEGs[c(1:5,11:18,6,9)]
# colnames(Cyp79A2_DEGs) <- c("locus", "short_name", "name", "aliases", "length",
#                             "fpkm_WT", "fpkmSE_WT", "tmm_WT", "tmmSE_WT",
#                             "fpkm_Cyp79A2", "fpkmSE_Cyp79A2", "tmm_Cyp79A2", "tmmSE_Cyp79A2", "log2FC_Cyp79A2", "pValue_Cyp79A2")
# 
# # Add tmm and fpkm to PAOx DEG dataframe ----
# fpkm_WT <- vector(mode="numeric", length=nrow(PAOx_DEGs))
# fpkmSE_WT <- vector(mode="numeric", length=nrow(PAOx_DEGs))
# tmm_WT <- vector(mode="numeric", length=nrow(PAOx_DEGs))
# tmmSE_WT <- vector(mode="numeric", length=nrow(PAOx_DEGs))
# fpkm_PAOx <- vector(mode="numeric", length=nrow(PAOx_DEGs))
# fpkmSE_PAOx <- vector(mode="numeric", length=nrow(PAOx_DEGs))
# tmm_PAOx <- vector(mode="numeric", length=nrow(PAOx_DEGs))
# tmmSE_PAOx <- vector(mode="numeric", length=nrow(PAOx_DEGs))
# 
# for(i in 1:nrow(PAOx_DEGs)){
#   geneLocus <- PAOx_DEGs$locus[i]
#   
#   fpkm_WT[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 8:10]))
#   fpkmSE_WT[i] <- sd(c(unlist(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 8:10])))/sqrt(3)
#   
#   tmm_WT[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 8:10]))
#   tmmSE_WT[i] <- sd(c(unlist(tmmAllGenes[tmmAllGenes$locus==geneLocus, 8:10])))/sqrt(3)
#   
#   fpkm_PAOx[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 5:7]))
#   fpkmSE_PAOx[i] <- sd(c(unlist(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 5:7])))/sqrt(3)
#   
#   tmm_PAOx[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 5:7]))
#   tmmSE_PAOx[i] <- sd(c(unlist(tmmAllGenes[tmmAllGenes$locus==geneLocus, 5:7])))/sqrt(3)
# }
# 
# PAOx_DEGs$fpkm_WT <- fpkm_WT
# PAOx_DEGs$fpkmSE_WT <- fpkmSE_WT
# 
# PAOx_DEGs$tmm_WT <- tmm_WT
# PAOx_DEGs$tmmSE_WT <- tmmSE_WT
# 
# PAOx_DEGs$fpkm_PAOx <- fpkm_PAOx
# PAOx_DEGs$fpkmSE_PAOx <- fpkmSE_PAOx
# 
# PAOx_DEGs$tmm_PAOx <- tmm_PAOx
# PAOx_DEGs$tmmSE_PAOx <- tmmSE_PAOx
# 
# PAOx_DEGs <- PAOx_DEGs[c(1,15:18,6,9)]
# colnames(PAOx_DEGs) <- c("locus",
#                          "fpkm_PAOx", "fpkmSE_PAOx", "tmm_PAOx", "tmmSE_PAOx", "log2FC_PAOx", "pValue_PAOx")
# 
# 
# # Merge DEG dataframes and save as an RDS ----
# DEGs_Cyp79A2 <- merge(Cyp79A2_DEGs, PAOx_DEGs, by="locus")
# for(col in 6:ncol(DEGs_Cyp79A2)){
#   DEGs_Cyp79A2[,col] <- as.numeric(formatC(DEGs_Cyp79A2[,col], width=5, format="G"))
# }
# saveRDS(DEGs_Cyp79A2, file="C:/Users/bca08_000/Documents/interactive-data/RNAseq/data/DEGs_Cyp79A2.rds")