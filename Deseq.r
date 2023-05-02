library(DESeq2)
library(apeglm)

# loading count matrix from D:\Research Work\Cancer\Temp_data\all_counts.csv

dat <- read.csv("D:/STAR/FEATURE_COUNTS.csv", header = TRUE, row.names = 1)

# loading colData from D:\Research Work\Cancer\Temp_data\colData.txt

info <- read.table("D:/STAR/Classification.txt", header = TRUE, sep = '\t')

dds <- DESeqDataSetFromMatrix(dat, info, ~ condition)

# removing lowly expressed genes
keep <- rowSums(counts(dds)) >= 0
dds <- dds[keep,]

# main analysis
ddsDE <- DESeq(dds)

#exporting normalized read counts
normCounts <- counts(ddsDE, normalized= T)
write.csv(normCounts, file = "D:/STAR/DESEQ2_normCounts.csv")

res <- results(ddsDE)
resOrdered <- res[order(res$padj),]
write.csv(resOrdered, file = "D:/STAR/FEATURE_COUNTS.csv")

summary(res)

# Plotting MA plot for p-value adjusted < 0.05 and padj > 0.05

library (ggplot2)
library (grid)

deSeqRes <- read.csv("D:/STAR/FEATURE_COUNTS.csv", row.names = 1)

# MA plot for p-value adjusted < 0.05 and padj > 0.05
deSeqRes$sig <- ifelse(deSeqRes$padj <= 0.05, "yes", "no")


# Volcano plot 
ggplot(deSeqRes, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = sig), size = 0.5) +
  scale_color_manual(values = c("red", "black")) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("log2 Fold Change") +
  ylab("-log10(padj)")

# Plot MA with range on y axis as -15 to 15 for log2FoldChange, also make a very thick 0 axis line, and remove the legend, legend as p-value adjusted <= 0.05 as red and padj > 0.05 as black
ggplot(deSeqRes, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = sig), size = 0.5) +
  scale_color_manual(values = c("red", "black")) +
  theme_bw() +
  xlab("log2 Fold Change") +
  ylab("-log10(padj)") +
  ylim(-15, 15) +
  geom_hline(yintercept = 0, color = "black", size = 2) 

library(org.Hs.eg.db)
res.df <- as.data.frame(res)
res.df$symbol <- mapIds(org.Hs.eg.db, keys = rownames(res.df), keytype = "ENSEMBL", column = "SYMBOL")
res.df

#Enhanced Volcano 
library('EnhancedVolcano')
EnhancedVolcano(res.df, x ="log2FoldChange", y = "padj",xlim = c(-1,1), ylim = c(0, 15), lab = res.df$symbol)


EnhancedVolcano(res.df,
                lab = res.df$symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                main = "Volcano Plot",
                selectLab = c('ACSM2A', 'OR4D9', 'OR2T7', 'VENTXP1', 'H3C2', 'SCAMP3', 'LDB2', 'SYNPO2'),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-14,
                FCcutoff = 2.0,
                pointSize = 2.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 3/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')


