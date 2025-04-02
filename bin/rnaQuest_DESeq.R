#This script performs differential expression analysis on the Pai et al. dataset (DOI: 10.1186/s40478-022-01453-1) using DESeq2.

library(ggplot2)
library(limma)
library(Glimma)
library(viridis)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DESeq2)
library(tidyverse)
library(treemap)
library(vsn)
library(plotly)
library(pheatmap)
library(EnhancedVolcano)
library(htmlwidgets)

#####################
###Data Management###
#####################

###Import read count matrix for TLE study
setwd(rstudioapi::selectDirectory(caption = "Select Project Directory"))
dat <- read.delim('data/merged.tsv', header = TRUE, sep="\t", check.names = FALSE)
dat[1:6,1:6]

#Remove special tagged rows, the first four, from the head of dat.
length(dat[,1])
head(dat)[,1:3]
dat <- dat[5:length(dat[,1]),]
head(dat)[,1:3]
length(dat[,1])
#Extract Ensembl ID and gene symbols
geneID <- as.data.frame(cbind(rownames(dat), dat$gene_symbols))
colnames(geneID) <- c('EnsID', 'gene_symbols')
head(geneID)
dat <- dat[,2:7]
head(dat)

#Import sample metadata.
metadata <- read.delim('data/SraRunMetadata.csv', header=TRUE, sep=",")
head(metadata)

#Metadata condition as factor
metadata$condition <- factor(metadata$condition)
levels(metadata$condition) <- c('TL','TLE')
metadata$ID <- c('Ctrl1', 'Ctrl2', 'Ctrl3', 'TLE1', 'TLE2', 'TLE3')

#Create sample information table
coldata <- metadata[,c('ID','condition', 'Run')]
rownames(coldata) <- metadata$ID
coldata
#Assign row names of coldata to column names of dat
colnames(dat) <- rownames(coldata)
head(dat)


#######################
###DESeq DE analysis###
#######################

#Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = dat,
                              colData = coldata,
                              design = ~ condition)
dds

#Replace Ensembl ID with gene symbol were available
IDs <- rownames(dds)
head(IDs)
head(geneID)
colnames(geneID)[1:2] <- c('ensemblID', 'geneID')
#check that IDs and geneID match each other
all(IDs == geneID$ensemblID)
length(which(IDs!=geneID$ensemblID))
length(IDs)
length(geneID$ensemblID)
tail(IDs)
tail(geneID$ensemblID)
#Create a new gene feature column that uses gene symbol, if available, and Ensembl ID otherwise.
geneID$geneNew <- ifelse(is.na(geneID$geneID), 
                         geneID$ensemblID, 
                         geneID$geneID)

#Check that geneNew shows desired output for missing gene symbols
head(geneID)
#Spot check random IDs and ensemblIDs for sanity check
rndm <- round(runif(10, min = 0, max = length(IDs)), digits = 0)
IDs[rndm]
geneID[rndm,]

#Add gene annotations
gene_info <- AnnotationDbi::select(org.Hs.eg.db, 
                    keys = geneID$ensemblID, 
                    columns = c("SYMBOL", "MAP", "GENETYPE", "GENENAME"), 
                    keytype = "ENSEMBL")
head(gene_info)
#return only unique mappings
gene_info <- gene_info[isUnique(gene_info$ENSEMBL),]
#geneID and gene_info do not match yet
all(geneID$ensemblID == gene_info$ENSEMBL)
gene_info <- merge(geneID, gene_info, 
                   by.x = "ensemblID",  # Column in colData to match
                   by.y = "ENSEMBL",    # Column in gene_info to match
                   all.x = TRUE)  # Keep all colData rows
# Restore original order of colData
gene_info <- gene_info[match(geneID$ensemblID, gene_info$ensemblID), ]
#Now they do match!
all(geneID$ensemblID == gene_info$ENSEMBL)
head(geneID)
head(gene_info)

#Add geneID to dds metadata
mcols(dds) <- cbind(mcols(dds), gene_info)
rownames(dds) <- geneID$geneNew

#Pre-filter low count genes
#Smallest group size is the number of samples in each group (3 in this ).
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

#Confirm factor levels
dds$ID
dds$condition
dds$Run
# mcols(dds)$MAP

#Differential expression analysis
dds <- DESeq(dds)
#Remove nonconverged genes
dds <- dds[which(mcols(dds)$betaConv),]
boxplot(log2(counts(dds, normalized = TRUE)+1), axes = FALSE, ylab = 'log2 counts')
axis(1, labels = FALSE, at = c(seq(1,40,1)))
axis(2, labels = TRUE)
text(x = seq_along(colnames(dds)), y = -2, labels = colnames(dds), 
     # srt = 45,    #rotate
     adj = 0.5,    #justify
     xpd = TRUE, #print outside of plot area
     cex = 0.75)  #smaller font size

#Plot dispersion estimates
plotDispEsts(dds)

#Save an rds file to start the script here
# saveRDS(dds, 'bin/dds_DE_TL-TLE.rds')
# dds <- readRDS('bin/dds_DE_TL-TLE.rds')


#########################
###Build results table###
#########################

levels(dds$condition)
#Set contrast groups, reference level should be listed last.
contrast <- c('condition', 'TLE', 'TL')
res <- results(dds,
               contrast = contrast,
               alpha = 0.1,
               pAdjustMethod = 'BH')
res
summary(res)

#How many adjusted p-values are less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)
table(res$padj < 0.1)
#Create significant results table ordered by fold change
res10 <- res[which(res$padj < 0.1),]
res10 <- res10[order(-res10$log2FoldChange),]
#How many adjusted p-values are less than 0.05?
sum(res$padj < 0.05, na.rm=TRUE)
table(res$padj < 0.05)

#How many up-regulated (FDR < 0.1)
length(which(res10$log2FoldChange > 0))
#How many down-regulated (FDR < 0.1)
length(which(res10$log2FoldChange < 0))
bar <- data.frame(direction = c('up','down'), genes = c(length(which(res10$log2FoldChange > 0)),length(which(res10$log2FoldChange < 0))))
bar

#Log fold change shrinkage for visualization and ranking, especially since we only have 3 samples per group
plotMA(res, ylim=c(-3,3))
resultsNames(dds)
resLFC <- lfcShrink(dds, contrast = contrast, type = 'ashr')
resLFC
plotMA(resLFC, ylim=c(-3,3))
res <- resLFC

#Add geneID variables
res$geneNew <- row.names(res)
res_annotated <- merge(as.data.frame(res), gene_info, by="geneNew", all.x=TRUE)
head(mcols(dds)[3])
head(res_annotated)

#Write out results table
#Reorder columns and order results table by smallest p-value
colnames(res_annotated)[1] <- 'gene'
colnames(res_annotated)[6] <- 'FDR'
res_annotated <- res_annotated[order(res_annotated$FDR),]
head(res_annotated)
# write.csv(as.data.frame(res_annotated), file=paste0('results/','TL','_','TLE','_allGenes_DEresults.csv'), row.names = TRUE)

#Glimma MD plot
res.df <- as.data.frame(res_annotated)
res.df$log10MeanNormCount <- log10(res.df$baseMean)
idx <- rowSums(counts(dds)) > 0
res.df <- res.df[idx,]
res.df$FDR[is.na(res.df$FDR)] <- 1
res.df <- res.df[isUnique(res.df$gene),]
exm <- res.df[,c(1:6,13)]
rownames(exm) <- exm$gene
anno <- res.df[,c(1,7:12)]
rownames(anno) <- anno$gene
# color_vector <- viridis(100)[cut(exm$FDR, breaks=100)]
# de_colors <- setNames(viridis(3), c('DE', 'Not Sig', 'Not DE'))
de_colors <- viridis(3)
status <- ifelse(exm$FDR < 0.05 & abs(exm$log2FoldChange) > 1, 'DE',
                 ifelse(exm$FDR > 0.05 & abs(exm$log2FoldChange) < 1, 'Not Sig', 'Not DE'))

glMDPlot(exm,
         xval = 'log10MeanNormCount',
         yval = 'log2FoldChange',
         anno = anno,
         groups = dds$condition,
         # cols = de_colors,
         cols = c('#BEBEBEFF', '#440154FF', '#21908CFF'),
         status = status,
         samples = colnames(dds)
)


#########
###EDA###
#########

#Create treemap of gene biotypes for all genes
biotype_all <- data.frame(table(mcols(dds)$GENETYPE))
biotype_all$Var1 <- gsub('_',' ',biotype_all$Var1)
biotype_all <- biotype_all[order(biotype_all$Freq, decreasing = TRUE),]
# png(filename='results/figs/tree.png',width=800, height=1295)
treemap(biotype_all,
        index = 'Var1',
        vSize = 'Freq',
        type = 'index',
        palette = viridis(length(biotype_all$Var1)),
        aspRatio = 1/1.618,
        title = '',
        inflate.labels = TRUE,
        lowerbound.cex.labels = 0
)
# dev.off()

#Create treemap of gene biotypes for significant results genes
res_annotated$GENETYPE <- as.factor(res_annotated$GENETYPE)
levels(res_annotated$GENETYPE)
biotype_sig <- data.frame(
  table(
    res_annotated$GENETYPE[which(res_annotated$FDR < 0.05)]
    )
  )
biotype_sig$Var1 <- gsub('_',' ',biotype_sig$Var1)
biotype_sig <- biotype_sig[order(biotype_sig$Freq, decreasing = TRUE),]
# png(filename='results/figs/tree.png',width=800, height=1295)
treemap(biotype_sig,
        index = 'Var1',
        vSize = 'Freq',
        type = 'index',
        palette = viridis(length(biotype_sig$Var1)),
        aspRatio = 1/1.618,
        title = '',
        inflate.labels = TRUE,
        lowerbound.cex.labels = 0
)
# dev.off()
biotype <- merge(biotype_all, biotype_sig, by = 'Var1')
colnames(biotype) <- c('biotype','AllGenes','SigGenes')
biotype <- biotype[order(biotype$AllGenes, decreasing = TRUE),]
# Combine rows with low counts into 'other' for the Chi-squared test.
biotype <- rbind(biotype,c('other', sum(biotype$AllGenes[4:8]), sum(biotype$SigGenes[4:8])))
biotype <- biotype[c(1:3,9),]
biotype$AllGenes <- as.integer(biotype$AllGenes)
biotype$SigGenes <- as.integer(biotype$SigGenes)
biotype$AllProp <- round(biotype$AllGenes/sum(biotype$AllGenes), digits = 3)
biotype$SigProp <- round(biotype$SigGenes/sum(biotype$SigGenes), digits = 3)
biotype <- biotype[,c(1,2,4,3,5)]
biotype
# Perform Chi-Squared Goodness-of-Fit Test
# Compute total counts
total_all <- sum(biotype$AllGenes)  # Total genes
total_sig <- sum(biotype$SigGenes)  # Total significant genes
# Compute expected counts under the null hypothesis
biotype$Expected <- total_sig * (biotype$AllGenes / total_all)
chi_sq_test <- chisq.test(biotype$SigGenes, p = biotype$AllGenes / total_all, rescale.p = TRUE)
# Print results
print(chi_sq_test)

#VST for visualization and ranking, generally scales well for large datasets
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
meanSdPlot(assay(vsd))
#NTD for visualization and ranking 
ntd <- normTransform(dds)
head(assay(ntd), 3)
meanSdPlot(assay(ntd))
#rlog for visualization and ranking, generally works well for small datasets
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
meanSdPlot(assay(rld))
#I'll use rld for visualization moving forward, because of the small sample size.

#Dimensional reduction EDA
#PCA plot
plotPCA(rld, intgroup='condition')
plotPCA(rld, intgroup='condition', pcsToUse=3:4)

#MDS plot
# rld <- calcNormFactors(rld)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=condition,shape=condition)) + geom_point(size=3)

#######################
###Visualize results###
#######################

#Create a function that will output a standalone html file of a plotly figure for integration into the subsequent markdown report.
# setwd('results/figs/')
# outdir <- getwd()
saveInteractivePlot <- function(plot, outdir, filename) {
  # Ensure the output directory exists
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  # Construct full file path
  file_path <- file.path(outdir, filename)
  # Save the widget
  saveWidget(plot, file_path, selfcontained = TRUE)
  # Read the saved HTML file
  html_content <- readLines(file_path)
  # Write the modified HTML file back
  writeLines(html_content, file_path)
  message("Interactive plot saved successfully: ", file_path)
}

#Data Quality assessment by sample clustering and visualization
ann_colors = viridis(3, begin = 0, end = 1, direction = -1)
ann_colors = list(Condition = c(TL=ann_colors[2], TLE=ann_colors[3]))
# df <- as.data.frame(colData(dds)[,'condition'])
df <- data.frame(Condition = colData(dds)[, 'condition'])
rownames(df) <- colnames(dds)  # Ensure row names match colnames in heatmap

select <- order(rowMeans(counts(dds, normalized = TRUE)),
                decreasing = TRUE)[1:10]
#Declare expression matrix
# rld <- rld[,c(which(rld$condition=='TL'),which(rld$condition=='TLE'))]
rld <- rld[,c(which(rld$condition==contrast[2]),which(rld$condition==contrast[3]))]
# Extract sample names in the correct order (TL first, then TLE)
ordered_samples <- colnames(rld)[order(rld$condition)]
pheatmap(assay(rld)[select,ordered_samples],
         cluster_rows = TRUE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         annotation_col = df,
         annotation_colors = ann_colors,
         color = viridis(n=256, begin = 0, end = 1, direction = -1),
         labels_row = mcols(rld)$geneNew
         )

# select <- assay(rld)[head(order(res$FDR),100),]
#Up-regulated genes
geneNum <- 40
selectUp <- assay(rld)[rownames(res10),ordered_samples][1:geneNum,]
selectUpNames <- rownames(res10)[1:geneNum]
selectUpNames <- which(rownames(dds)%in%selectUpNames)
zUp <- (selectUp - rowMeans(selectUp))/sd(selectUp) #z-score
pheatmap(zUp,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         annotation_col = df,
         annotation_colors = ann_colors,
         labels_row = mcols(rld)$geneNew[selectUpNames],
         color = viridis(n=256, begin = 0, end = 1, direction = -1))

#Down-regulated genes
geneNum <- 40
selectDown <- assay(rld)[rownames(res10),ordered_samples][(length(res10[,1])-geneNum):length(res10[,1]),]
selectDownNames <- rownames(res10)[(length(res10[,1])-geneNum):length(res10[,1])]
selectDownNames <- which(rownames(dds)%in%selectDownNames)
zDown <- (selectDown - rowMeans(selectDown))/sd(selectDown) #z-score
pheatmap(zDown,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         annotation_col = df,
         annotation_colors = ann_colors,
         labels_row = mcols(rld)$geneNew[selectDownNames],
         color = viridis(n=256, begin = 0, end = 1, direction = -1))

#Up and Down-regulated genes
selectUpDownNames <- c(selectUpNames,selectDownNames)
zUpDown <- rbind(zUp, zDown) #z-score
# setwd('results/figs/')
# tiff(filename= paste0('TLE','_','TL','_upDownHeatmap.tiff'),width=1200, height=1000)
gap_index <- sum(df$Condition == "TL")
pheatmap(zUpDown,
         fontsize = 8,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         annotation_col = df,
         annotation_colors = ann_colors,
         # gaps_col = gap_index,
         labels_row = mcols(rld)$geneNew[selectUpDownNames],
         labels_col = "",
         show_colnames = TRUE,
         color = viridis(n=256, begin = 0, end = 1, direction = -1))
# dev.off()

#Volcano plot
EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue',
                # title = paste0(contrast[2], ' vs. ', contrast[3]),
                title = bquote('TL vs. TLE'),
                subtitle = NULL,
                legendPosition = 'none',
                pCutoff = 10e-6,
                FCcutoff = 1,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                colConnectors = 'black',
                widthConnectors = 1.0,
                arrowheads = FALSE,
                col = c('grey30', viridis(3, direction = -1)), colAlpha = 0.7,
                gridlines.major = FALSE,
                gridlines.minor = FALSE) #+ coord_flip()


#Interactive volcano plot
# Set thresholds
fdr_cutoff <- 0.001
logfc_cutoff <- 1
# Set colors
sig_colors <- viridis(3, direction = -1)
# Set tooltip text
res_annotated$tooltip <- paste0(
  res_annotated$gene,
  "<br>log₂FC: ", round(res_annotated$log2FoldChange, 3),
  "<br>FDR: ", signif(res_annotated$FDR, 3)
)

res_annotated$Significant <- ifelse(
  res_annotated$FDR < fdr_cutoff & abs(res_annotated$log2FoldChange) > logfc_cutoff, "Sig",
  ifelse(
    res_annotated$FDR > fdr_cutoff & abs(res_annotated$log2FoldChange) >= logfc_cutoff, "nonDE",
    "nonSig"
  )
)

volcano <- ggplot(res_annotated, aes(x = log2FoldChange, y = -log10(pvalue), text = tooltip)) +
  geom_point(aes(color = Significant), alpha = 0.7) +
  scale_color_manual(values = c("nonSig" = "gray40", "Sig" = sig_colors[3], "nonDE" = sig_colors[1])) +
  labs(
    x = "log₂ Fold Change",
    y = "-log₁₀(p-value)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  ) +
  theme(legend.position = "none") +
  # Horizontal line (significance)
  geom_hline(
    yintercept = 5,
    linetype = "dashed",
    color = "black",
    linewidth = 0.3
  ) +
  # Vertical lines (logfc)
  geom_vline(
    xintercept = -logfc_cutoff,
    linetype = "dashed",
    color = "black",
    linewidth = 0.3
  ) +
  geom_vline(
    xintercept = logfc_cutoff,
    linetype = "dashed",
    color = "black",
    linewidth = 0.3
  )

ggplotly(volcano, tooltip = "text")


#Plot counts highest foldChange gene
plotCounts(dds, gene = rownames(res)[which.max(res$log2FoldChange)], intgroup = 'condition', normalized = TRUE, transform = TRUE)

plt <- plotCounts(dds, gene = rownames(res)[which.max(res$log2FoldChange)], intgroup = 'condition', normalized = TRUE, transform = TRUE, returnData=TRUE)

ggplot(plt, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) +
  labs(title = rownames(res)[which.max(res$log2FoldChange)]) +
  # labs(title = rownames(res)[which.min(res$padj)]) +
  theme(plot.title = element_text(hjust = 0.5)) +
  #scale_y_log10(breaks=c(25,100,400))
  scale_y_log10()

#Create annotation function for significance bracket
annotate_significance <- function(p) {
  # Extract gene of interest data
  y_max <- max(plt$count, na.rm = TRUE)
  padj_value <- res[ens, "padj"]
  # Compute positions for bracket and text
  y_bracket <- y_max * 1.30  # 5% above max expression
  y_text <- y_max * 1.5     # 10% above max
  # Add annotation to the plot
  p <- p +
    # Horizontal line (bracket)
    annotate("segment", x = 1, xend = 2, y = y_bracket, yend = y_bracket) +
    # Vertical ticks at each end
    annotate("segment", x = 1, xend = 1, y = y_bracket, yend = y_bracket * 0.90) +
    annotate("segment", x = 2, xend = 2, y = y_bracket, yend = y_bracket * 0.90) +
    # P-value text annotation
    annotate("text", x = 1.5, y = y_text, 
             label = paste("FDR =", format(padj_value, digits=3)))
  p  # Return the modified plot
}

#Plot counts gene of interest
geneOfInt <- 'LGR6'
ens <- which(res$geneNew == geneOfInt)
# plotCounts(dds, gene = rownames(res)[ens], intgroup = 'group')

plt <- plotCounts(dds, gene = rownames(res)[ens], intgroup = 'condition', normalized = TRUE, transform = TRUE, returnData=TRUE)

p <- ggplot(plt, aes(x=condition, y=count, color = condition)) + 
  # geom_boxplot(color = 'black', outlier.shape = NA) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_color_viridis_d(begin = 0.75, end = 0.25, direction = 1) +
  labs(title = geneOfInt) +
  # labs(title = rownames(res)[which.min(res$padj)]) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line('black')) +
  #scale_y_log10(breaks=c(25,100,400))
  scale_y_log10()
p

max(plt$count)
# setwd('results/figs/')
# tiff(filename= paste0(geneOfInt,'_plotCounts.tiff'),width=800, height=1295)
p <- ggplot(plt, aes(x=condition, y=count, color = condition)) + 
  geom_boxplot(color = 'black', outlier.shape = NA, width = 0.25) +
  geom_point(cex = 5, position=position_jitter(w=0.1,h=0)) +
  scale_color_viridis_d(begin = 0.75, end = 0.25, direction = 1) +
  labs(title = bquote(italic(.(geneOfInt))),
       y = 'Counts',
       x = NULL) +
  theme_bw() +
  theme(text = element_text(size=30),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line('black'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_log10()
p

#Add significance bracket
p + 
  annotate("text",
         x=1.5,
         y=max(plt$count)*1.5,
         label=paste("FDR =",format(res[ens, "padj"],digits=3))
         ) +
  annotate("segment", x = 1, xend = 2, y = max(plt$count) * 1.30, yend = max(plt$count) * 1.30, size = 0.5)

# dev.off()

annotate_significance(p)

p <- ggplot(plt, aes(x=condition, y=count, color = condition)) + 
  geom_boxplot(color = 'black', outlier.shape = NA, width = 0.25) +
  geom_point(cex = 5, position=position_jitter(w=0.1,h=0)) +
  scale_color_viridis_d(begin = 0.75, end = 0.25, direction = 1) +
  labs(title = geneOfInt,
       y = 'Counts',
       x = NULL) +
  theme_bw() +
  theme(text = element_text(size=30),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line('black'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_log10()
  # annotate("text",
  #          x=1.5,
  #          y=max(plt$count)*1.05,
  #          label=paste("FDR =",format(res[ens, "padj"],digits=3))
  #          )

ggplotly(p)
p_pltly <- ggplotly(annotate_significance(p))
p_pltly

#Save the interactive plot out as a standalone html file
#Uncomment the next line to run the function from source
# source("bin/rnaQuest_html.R")
saveInteractivePlot(p_pltly, 'results/figs/', paste0(geneOfInt, "_plotlyCounts.html"))

#Export an rld expression file for use in the gene expression dashboard app (dash app)

#Get expression data from rld
rld_matrix <- assay(rld)  # genes x samples
rld_long <- as.data.frame(rld_matrix) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression")

#Add sample metadata
sample_meta <- colData(dds) %>%
  as.data.frame() %>%
  rownames_to_column("sample")

rld_joined <- rld_long %>%
  left_join(sample_meta, by = "sample")

rld_joined <- rld_joined %>%
  mutate(gene = toupper(gene))

#Add DE results (log2FC, FDR)
de_results <- as.data.frame(res_annotated) %>%
  select(gene, log2FoldChange, FDR)

#Remove duplicates and keep the row with the lowest FDR
de_results <- de_results %>%
  filter(!is.na(FDR)) %>%
  group_by(gene) %>%
  slice_min(order_by = FDR, n = 1, with_ties = FALSE) %>%
  ungroup()

de_results <- de_results %>%
  mutate(gene = toupper(gene))

#Merge everything
final_df <- rld_joined %>%
  left_join(de_results, by = "gene")

final_df %>%
  filter(!is.na(FDR)) %>%
  select(gene, log2FoldChange, FDR) %>%
  distinct() %>%
  head()
#rows missing log2foldchange of FDR are not significantly DE.

#Save for Dash
write_csv(final_df, "results/dash_expression_data.csv")


#Gene set enrichment analysis
#Coming soon!
