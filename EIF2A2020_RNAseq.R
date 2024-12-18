library(tidyverse)
#BiocManager::install("DESeq2")
library("DESeq2")
#library(EnhancedVolcano)
library(org.Mm.eg.db)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) ##Set working directory to where this file is.
getwd()

con0_2 <- read.table(file="input/CON0-2.counts.genes") %>% rename(control0_2=V2)
con0_3 <- read.table(file="input/CON0-3.counts.genes") %>% rename(control0_3=V2)
con0_4 <- read.table(file="input/CON0-4.counts.genes") %>% rename(control0_4=V2)
con0_5 <- read.table(file="input/CON0-5.counts.genes") %>% rename(control0_5=V2)

con1_2 <- read.table(file="input/CON1-2.counts.genes") %>% rename(control1_2=V2)
con1_3 <- read.table(file="input/CON1-3.counts.genes") %>% rename(control1_3=V2)
con1_4 <- read.table(file="input/CON1-4.counts.genes") %>% rename(control1_4=V2)
con1_5 <- read.table(file="input/CON1-5.counts.genes") %>% rename(control1_5=V2)

con3_2 <- read.table(file="input/CON3-2.counts.genes") %>% rename(control3_2=V2)
con3_3 <- read.table(file="input/CON3-3.counts.genes") %>% rename(control3_3=V2)
con3_4 <- read.table(file="input/CON3-4.counts.genes") %>% rename(control3_4=V2)
con3_5 <- read.table(file="input/CON3-5.counts.genes") %>% rename(control3_5=V2)

eIF2A0_2 <- read.table(file="input/eIF2A0-2.counts.genes") %>% rename(eIF2A0_2=V2)
eIF2A0_3 <- read.table(file="input/eIF2A0-3.counts.genes") %>% rename(eIF2A0_3=V2)
eIF2A0_4 <- read.table(file="input/eIF2A0-4.counts.genes") %>% rename(eIF2A0_4=V2)
eIF2A0_5 <- read.table(file="input/eIF2A0-5.counts.genes") %>% rename(eIF2A0_5=V2)

eIF2A1_2 <- read.table(file="input/eIF2A1-2.counts.genes") %>% rename(eIF2A1_2=V2)
eIF2A1_3 <- read.table(file="input/eIF2A1-3.counts.genes") %>% rename(eIF2A1_3=V2)
eIF2A1_4 <- read.table(file="input/eIF2A1-4.counts.genes") %>% rename(eIF2A1_4=V2)
eIF2A1_5 <- read.table(file="input/eIF2A1-5.counts.genes") %>% rename(eIF2A1_5=V2)

eIF2A3_2 <- read.table(file="input/eIF2A3-2.counts.genes") %>% rename(eIF2A3_2=V2)
eIF2A3_3 <- read.table(file="input/eIF2A3-3.counts.genes") %>% rename(eIF2A3_3=V2)
eIF2A3_4 <- read.table(file="input/eIF2A3-4.counts.genes") %>% rename(eIF2A3_4=V2)
eIF2A3_5 <- read.table(file="input/eIF2A3-5.counts.genes") %>% rename(eIF2A3_5=V2)

raw.counts <- bind_cols(con0_2,con0_3,con0_4,con0_5,con1_2,con1_3,con1_4,con1_5,con3_2,con3_3,con3_4,con3_5,
                        eIF2A0_2,eIF2A0_3,eIF2A0_4,eIF2A0_5,eIF2A1_2,eIF2A1_3,eIF2A1_4,eIF2A1_5,eIF2A3_2,eIF2A3_3,eIF2A3_4,eIF2A3_5)
dim(raw.counts)
head(raw.counts)
rownames(raw.counts) <- raw.counts$V1...1
raw.counts <- raw.counts[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48)]
write.table(raw.counts, sep="\t",file="input/rawcounts.txt", 
          row.names=TRUE,col.names=NA,quote=FALSE)

sample= c("control0_2","control0_3","control0_4","control0_5","control1_2","control1_3","control1_4","control1_5","control3_2","control3_3","control3_4","control3_5",
         "eIF2A0_2","eIF2A0_3","eIF2A0_4","eIF2A0_5","eIF2A1_2","eIF2A1_3","eIF2A1_4","eIF2A1_5","eIF2A3_2","eIF2A3_3","eIF2A3_4","eIF2A3_5")
treatment=c("control","control","control","control","control","control","control","control","control","control","control","control",
       "eIF2A","eIF2A","eIF2A","eIF2A","eIF2A","eIF2A","eIF2A","eIF2A","eIF2A","eIF2A","eIF2A","eIF2A")
time=c(0,0,0,0,1,1,1,1,3,3,3,3,0,0,0,0,1,1,1,1,3,3,3,3)
group=paste(treatment,time,sep="_")
meta.data <- data.frame(sample,treatment,time,group)
rownames(meta.data) <- meta.data$sample
head(meta.data)
write.table(meta.data, sep="\t",file="input/metadata.txt", 
            row.names=TRUE,col.names=NA,quote=FALSE)

meta.data <- read.table("input/metadata.txt", row.names = 1)
raw.counts <- read.table("input/GSE166829_rawcounts.txt")
View(raw.counts)
View(meta.data)

### create DESeq matrix
count.data.set = DESeqDataSetFromMatrix(countData=raw.counts, 
                                        colData=meta.data, design= ~ group) 

keep <- rowSums(counts(count.data.set)>=5) >= ncol(raw.counts)*0.25 # genes counts more than 5 in at least 10 samples
count.filter <- count.data.set[keep,]
nrow(count.filter) # 14486 genes

# create DESeq object
#count.data.set.object <- DESeq(count.data.set)
count.data.set.object <- DESeq(count.filter)
count.data.set.object

# normalized counts (without VST)
dds <- estimateSizeFactors(count.data.set.object)
normalized.counts <- counts(dds, normalized=TRUE)
write.csv(normalized.counts, file="output/norm_counts.csv")

# 'vst' normalization (varianceStabilizingTransformation)
vsd <- vst(count.data.set.object)

### extract normalized counts
norm.data = assay(vsd)
head(norm.data)
write.table(norm.data, sep="\t",file="output/Norm_data_all_genes_noFilter.txt", 
            row.names=TRUE,col.names=NA,quote=FALSE)

write.table(norm.data, sep="\t",file="output/Norm_data_all_genes_Filter.txt", 
            row.names=TRUE,col.names=NA,quote=FALSE)

### hierarchical clustering analyses and to plot a dendrogram. Evaluate dissimilarities (calculate Euclidean distance) between all eight replicates based on their normalized gene counts.
sampleDists <- dist(t(norm.data),  method = "euclidean")
### Having the distance (dissimilarity) we can finally perform hierarchical cluster analysis using hclust function
clusters=hclust(sampleDists)
plot(clusters)
### plots PCA for the first two principal components
getMethod("plotPCA","DESeqTransform")
plotPCA.format <- function (object, ...) 
{
  .local <- function (object, intgroup = "condition", 
                      ntop = 500, returnData = FALSE) 
  {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                       length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                                 drop = FALSE])
    group <- if (length(intgroup) > 1) {
      factor(apply(intgroup.df, 1, paste, collapse = ":"))
    }
    else {
      colData(object)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                    intgroup.df, name = colnames(object))
    if (returnData) {
      attr(d, "percentVar") <- percentVar[1:2]
      return(d)
    }
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", 
                                fill = "group")) + geom_point(size = 4, shape=21,alpha=0.6) + 
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
      coord_fixed() #+  geom_label_repel((aes(label=sample)))
  }
  .local(object, ...)
}


#------------------------------------------------------------------------

# use color blind friendly palette
cbPalette <- c("#E69F00", #lightorange
               "#56B4E9", #blue
               "#D55E00", #darkorange
               "#009E73", #green
               "#CC79A7", #magenta
               "#0072B2", #darkblue
               "#F0E442", #yellow
               "#999999" #grey
)

#install.packages("ggrepel")
library(ggrepel) # https://ggrepel.slowkow.com/articles/examples.html

plotPCA.format(vsd, intgroup=c("group"))+ 
  geom_text_repel(aes(label=colData(vsd)$sample),size=3,
                  color="grey50",
                  box.padding   = 0.4,
                  point.padding = 0,
                  #force=1,
                  #force_pull=10,
                  max.overlaps = Inf, # always show all label, regardless of overlap
                  #min.segment.length = 0, # always draw line
                  segment.color = 'grey50')+
  scale_fill_manual(values=cbPalette) +
  theme_bw()+
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x  = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank()
  )+
  theme(aspect.ratio=1/1)
# FYI, this is the default function plotPCA()
#plotPCA(vsd, intgroup=c("group"))

ggsave(filename="figures/PCA2.png",width=9,height=7,units="cm",dpi=400)



library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix( sampleDists )
#rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
#colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pDist<-pheatmap(sampleDistMatrix,
                clustering_distance_rows = sampleDists,
                clustering_distance_cols = sampleDists,
                col = colors)
save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(pDist, "pDist.png")
dev.off()

### DEG

res <- results(count.data.set.object, 
               contrast=c("group","eIF2A_3","control_3"),
               alpha = 0.05)

res <- results(count.data.set.object, 
               contrast=c("group","eIF2A_1","control_1"),
               alpha = 0.05)

res <- results(count.data.set.object, 
               contrast=c("group","eIF2A_0","control_0"),
               alpha = 0.05)



res <- results(count.data.set.object, 
               contrast=c("group","eIF2A_1","eIF2A_0"),
               alpha = 0.05)

res <- results(count.data.set.object, 
               contrast=c("group","control_1","control_0"),
               alpha = 0.05)

res <- results(count.data.set.object, 
               contrast=c("group","eIF2A_3","eIF2A_0"),
               alpha = 0.05)

res <- results(count.data.set.object, 
               contrast=c("group","control_3","control_0"),
               alpha = 0.05)

summary(res)
out <- capture.output(summary(res))

cat("eIF2A_0 vs control_0", out, file="output/results_summary.txt", sep="\n", append=TRUE)
cat("eIF2A_1 vs control_1", out, file="output/results_summary.txt", sep="\n", append=TRUE)
cat("eIF2A_3 vs control_3", out, file="output/results_summary.txt", sep="\n", append=TRUE)


cat("eIF2A_0 vs control_0 filtered", out, file="output/results_summary.txt", sep="\n", append=TRUE)
cat("eIF2A_1 vs control_1 filtered", out, file="output/results_summary.txt", sep="\n", append=TRUE)
cat("eIF2A_3 vs control_3 filtered", out, file="output/results_summary.txt", sep="\n", append=TRUE)

cat("eIF2A_1 vs eIF2A_0 filtered", out, file="output/results_summary.txt", sep="\n", append=TRUE)
cat("control_1 vs control_0 filtered", out, file="output/results_summary.txt", sep="\n", append=TRUE)

cat("eIF2A_3 vs eIF2A_0 filtered", out, file="output/results_summary.txt", sep="\n", append=TRUE)
cat("control_1 vs control_0 filtered", out, file="output/results_summary.txt", sep="\n", append=TRUE)

##
##

res = as.data.frame(na.omit(res))
res.ord = res[order(res$padj),]
##
res.ord$symbol=rownames(res.ord)
res.ord$entrez = mapIds(org.Mm.eg.db, keys=rownames(res.ord), column="ENTREZID", keytype="SYMBOL", multiVals="first") %>%
  as.character()
res.ord$description = mapIds(org.Mm.eg.db, keys=res.ord$entrez, column="GENENAME", keytype="ENTREZID", multiVals="first") %>%
  as.character()

#res.ord$description <- gsub(x = res.ord$description, pattern = "\\ ",replacement = "_")
View(res.ord)
write.csv(res.ord, file="output/Results_Eif2a3_ctrl3.csv",row.names=FALSE)
write.csv(res.ord, file="output/Results_Eif2a1_ctrl1.csv",row.names=FALSE)
write.csv(res.ord, file="output/Results_Eif2a0_ctrl0.csv",row.names=FALSE)

write.csv(res.ord, file="output/Results_Eif2a1_Eif2a0.csv",row.names=FALSE)
write.csv(res.ord, file="output/Results_ctrl1_ctrl0.csv",row.names=FALSE)
write.csv(res.ord, file="output/Results_Eif2a3_Eif2a0.csv",row.names=FALSE)
write.csv(res.ord, file="output/Results_ctrl3_ctrl0.csv",row.names=FALSE)

##
write.table(res.ord, sep="\t",file="output/Results_Eif2a3_ctrl3.txt",row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.ord, sep="\t",file="output/Results_Eif2a1_ctrl1.txt",row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.ord, sep="\t",file="output/Results_Eif2a0_ctrl0.txt",row.names=TRUE,col.names=NA,quote=FALSE)
##
##


res.sig = res.ord[res.ord$padj <= 0.05,]
head(res.sig)

res.sig$entrez = mapIds(org.Mm.eg.db, keys=rownames(res.sig), column="ENTREZID", keytype="SYMBOL", multiVals="first")
res.sig$symbol=rownames(res.sig)
res.sig$description = mapIds(org.Mm.eg.db, keys=res.sig$entrez, column="GENENAME", keytype="ENTREZID", multiVals="first")
res.sig$description <- gsub(x = res.sig$description, pattern = "\\ ",replacement = "_")

res.sig <- as.data.frame(res.sig)
View(res.sig)

dup=res.sig %>% group_by(entrez) %>% dplyr::filter(n() > 1)
View(dup)

write.table(res.sig, sep="\t",file="data/Results_sig_3.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.sig, sep="\t",file="data/Results_sig_1.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.sig, sep="\t",file="data/Results_sig_0.txt", row.names=TRUE,col.names=NA,quote=FALSE)


Sig.up = subset(res.sig, log2FoldChange >0)
Sig.down = subset(res.sig, log2FoldChange <0)
##
##
write.table(Sig.up, sep="\t",file="data/Sig_up_3.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.table(Sig.down, sep="\t",file="data/Sig_down_3.txt", row.names=TRUE,col.names=NA,quote=FALSE)

write.table(Sig.up, sep="\t",file="data/Sig_up_1.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.table(Sig.down, sep="\t",file="data/Sig_down_1.txt", row.names=TRUE,col.names=NA,quote=FALSE)

write.table(Sig.up, sep="\t",file="data/Sig_up_0.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.table(Sig.down, sep="\t",file="data/Sig_down_0.txt", row.names=TRUE,col.names=NA,quote=FALSE)

##
##
# make a gene logFC list for networkanalyst.ca
#res.sig = read.table(file="data/Results_sig_3.txt", row.names = 1)
head(res.sig)
resSig.sub = dplyr::select(res.sig, symbol,log2FoldChange)
head(resSig.sub)
##
##
write.table(resSig.sub, sep="\t",file="data/Results_sig_3_networkanalyst.txt", row.names=FALSE,quote = FALSE)
write.table(resSig.sub, sep="\t",file="data/Results_sig_1_networkanalyst.txt", row.names=FALSE,quote = FALSE)
write.table(resSig.sub, sep="\t",file="data/Results_sig_0_networkanalyst.txt", row.names=FALSE,quote = FALSE)
##
##
library("clusterProfiler")
res.sig0 <- read.table(file="data/Results_sig_0.txt",row.names = 1)
res.sig1 <- read.table(file="data/Results_sig_1.txt",row.names = 1)
res.sig3 <- read.table(file="data/Results_sig_3.txt",row.names = 1)
head(res.sig0)

norm.data <- read.table(file="data/Norm_data_all_genes.txt",row.names = 1)
norm.data$entrez <- mapIds(org.Mm.eg.db, keys=rownames(norm.data), column="ENTREZID", keytype="SYMBOL", multiVals="first")
View(norm.data)
dup=norm.data %>% group_by(entrez) %>% dplyr::filter(n() > 1)
View(dup)
nrow(dup)
all.genes <- norm.data$entrez

KEGG<- enrichKEGG(gene         = res.sig1$entrez,
                      organism     = 'mmu',
                      universe      = all.genes,
                      pvalueCutoff=0.05, pAdjustMethod="BH", 
                      qvalueCutoff=0.1)
KEGG <- setReadable(KEGG, OrgDb = org.Mm.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
head(KEGG)
KEGG.df=as.data.frame(KEGG)
KEGG.df$Description <- gsub(x = KEGG.df$Description, pattern = "\\ ",replacement = "_")
View(KEGG.df)


KEGG<- enrichKEGG(gene         = res.sig3$entrez,
                  organism     = 'mmu',
                  universe      = all.genes,
                  pvalueCutoff=0.05, pAdjustMethod="BH", 
                  qvalueCutoff=0.1)
KEGG <- setReadable(KEGG, OrgDb = org.Mm.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
head(KEGG)
KEGG.df=as.data.frame(KEGG)
KEGG.df$Description <- gsub(x = KEGG.df$Description, pattern = "\\ ",replacement = "_")
View(KEGG.df)



KEGG<- enrichKEGG(gene         = res.sig0$entrez,
                  organism     = 'mmu',
                  universe      = all.genes,
                  pvalueCutoff=0.05, pAdjustMethod="BH", 
                  qvalueCutoff=0.1)
KEGG <- setReadable(KEGG, OrgDb = org.Mm.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
head(KEGG)
KEGG.df=as.data.frame(KEGG)
KEGG.df$Description <- gsub(x = KEGG.df$Description, pattern = "\\ ",replacement = "_")
View(KEGG.df)

##
##
write.table(KEGG.df, sep="\t",file="data/KEGG_3.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.csv(KEGG.df, sep="\t",file="data/KEGG_3.csv", row.names=TRUE,col.names=NA,quote=FALSE)

write.table(KEGG.df, sep="\t",file="data/KEGG_1.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.csv(KEGG.df, sep="\t",file="data/KEGG_1.csv", row.names=TRUE,col.names=NA,quote=FALSE)

write.table(KEGG.df, sep="\t",file="data/KEGG_0.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.csv(KEGG.df, sep="\t",file="data/KEGG_0.csv", row.names=TRUE,col.names=NA,quote=FALSE)
##
##

#=================================== plot KEGG network (not in paper)


geneList = res.sig0[,2]
names(geneList) = as.character(res.sig0$symbol)

geneList = res.sig1[,2]
names(geneList) = as.character(res.sig1$symbol)

geneList = res.sig3[,2]
names(geneList) = as.character(res.sig3$symbol)

geneList = sort(geneList, decreasing = TRUE)
head(geneList)

View(KEGG.df)
#
p <- cnetplot(KEGG, showCategory=14,
              node_label="gene",
              categorySize="pvalue", colorEdge = TRUE,
              order=TRUE, by="pvalue", foldChange=geneList)+
  scale_color_gradient2(high='red', low='blue', midpoint=0, mid="white")+
  labs(color = "log2 fold change",size="Gene count")

cowplot::plot_grid(p,ncol=1, rel_widths=1)

ggsave(
  filename="cnet_KEGG_BS_BG.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8,
  height = 7,
  units = "in",
  dpi = 600,
  limitsize = F
)


## plot KEGG terms
KEGG0 <- read.table(file="data/KEGG_0.txt",row.names = 1)
KEGG1 <- read.table(file="data/KEGG_1.txt",row.names = 1)
KEGG3 <- read.table(file="data/KEGG_3.txt",row.names = 1)

KEGG0$Description <- gsub(x = KEGG0$Description, pattern = "\\_",replacement = " ")
KEGG1$Description <- gsub(x = KEGG1$Description, pattern = "\\_",replacement = " ")
KEGG3$Description <- gsub(x = KEGG3$Description, pattern = "\\_",replacement = " ")

colnames(KEGG0) <- paste(colnames(KEGG0),"0", sep = "_")
colnames(KEGG1) <- paste(colnames(KEGG1),"1", sep = "_")
colnames(KEGG3) <- paste(colnames(KEGG3),"3", sep = "_")
View(KEGG0)
KEGG.merge<- right_join(KEGG0, KEGG1, by =c("Description_0"="Description_1"))
KEGG.merge<- left_join(KEGG.merge, KEGG3, by =c("Description_0"="Description_3"))
View(KEGG.merge)
write.table(KEGG.merge, sep="\t",file="data/KEGG_merge1.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.csv(KEGG.merge, sep="\t",file="data/KEGG_merge1.csv", row.names=TRUE,col.names=NA,quote=FALSE)

KEGG.merge <- read.csv(file="data/KEGG_merge1.csv", row.names=1)

##
KEGG.m <- KEGG.merge %>% mutate(x1="0h") %>% mutate(x2="1h") %>% mutate(x3="3h") %>%
  mutate(Log10adj.P_0=-log10(p.adjust_0)) %>% 
  mutate(Log10adj.P_1=-log10(p.adjust_1)) %>%
  mutate(Log10adj.P_3=-log10(p.adjust_3))
View(KEGG.m)

p <- ggplot(KEGG.m, 
            aes(y = fct_reorder(Description_0, Log10adj.P_1))) + 
  geom_point(aes(size = Count_0, color = Log10adj.P_0, x=x1)) +
  geom_point(aes(size = Count_1, color = Log10adj.P_1, x=x2))+
  geom_point(aes(size = Count_3, color = Log10adj.P_3, x=x3))+
  theme_bw(base_size = 10) +
  theme(panel.grid.major = element_blank())+
  scale_color_gradient2(#limits=c(-8,5) ,
                        midpoint = 0, low = "blue", mid = "lightsalmon",
                        high = "red4", space = "Lab" )+
  ylab(NULL)+
  xlab(NULL)+
  xlim("0h","1h","3h")+
  theme(axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black",size=10),
        axis.text.x = element_text(colour = "black",size=10),
        legend.title = element_text(color = "black", size = 10))+
  labs(color = "-log10 adj. p",size="Gene count")
# legend.text = element_blank())+
# 
p
p <- ggplot(KEGG.m, 
            aes(y = fct_reorder(Description_0, Log10adj.P_1))) + 
  geom_point(aes(size = GeneRatio_0, color = Log10adj.P_0, x=x1)) +
  geom_point(aes(size = GeneRatio_1, color = Log10adj.P_1, x=x2))+
  geom_point(aes(size = GeneRatio_3, color = Log10adj.P_3, x=x3))+
  theme_bw(base_size = 10) +
  scale_color_gradient2(#limits=c(-8,5) ,
    midpoint = 0, low = "blue", mid = "white",
    high = "red", space = "Lab" )+
  ylab(NULL)+
  xlab(NULL)+
  xlim("0h","1h","3h")+
  theme(axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black",size=10),
        axis.text.x = element_text(colour = "black",size=10),
        legend.title = element_text(color = "black", size = 10))+
  labs(color = "-log10 adj. p",size="Gene ratio")


p
ggsave(
  filename="KEGG_BS-ST-AS_BG.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8,
  height = 7,
  units = "in",
  dpi = 600,
  limitsize = F
)


#================================== Compare 3 time points
library(VennDiagram)
t0 = read.table(file="data/Results_sig_0.txt", row.names = 1)
t1 = read.table(file="data/Results_sig_1.txt", row.names = 1)
t3 = read.table(file="data/Results_sig_3.txt", row.names = 1)

grid.newpage()
grid::grid.draw(VennDiagram::venn.diagram(list("0h"=t0$symbol,
                                               "1h"=t1$symbol,
                                               "3h"=t3$symbol), NULL)) #6.5x5.5in


venn <- get.venn.partitions(list("0h"=t0$symbol,
                                 "1h"=t1$symbol,
                                 "3h"=t3$symbol))
str(venn)
head(venn)

list <- select(venn,..set..,..values..)
View(list)
str(list)

#### separate nested lists, tranfer wide to long format
library(splitstackshape)
df_venn <- listCol_w(venn, "..values..", fill = NA)
View(df_venn)
venn_long <- df_venn %>% 
  tidyr::pivot_longer(cols=..values.._fl_0001:..values.._fl_1243,values_to = "gene")%>% 
  na.omit()%>% dplyr::select(!name)
View(venn_long)
write.csv(venn_long,"Venn_long.csv", row.names = FALSE)
write.table(venn_long,file="Venn_long.txt",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)

write.table(venn_long$gene,file="Venn_long_gene.txt",sep="\t",row.names=F,quote=FALSE)

#========== Filter out TFs that are not expressed in our RNAseq
TFall <- read.csv(file="data/node_table_TFall.csv")
TFall.p <- read.csv(file="data/node_table_TFall_protein.csv")
TFall.g <- read.csv(file="data/node_table_TFall_GTEx.csv")
head(TFall)
norm.data <- read.table(file="data/Norm_data_all_genes.txt",row.names = 1)
head(norm.data)
norm.data$entrez <- mapIds(org.Mm.eg.db, keys=rownames(norm.data), column="ENTREZID", keytype="SYMBOL", multiVals="first")
norm.data.noNA <- na.omit(norm.data)

grid.newpage()
grid::grid.draw(VennDiagram::venn.diagram(list("TF"=TFall$Id,
                                               "norm_counts"=norm.data.noNA$entrez), NULL)) #6.5x5.5in
venn <- get.venn.partitions(list("TF"=TFall$Id,
                                 "counts"=norm.data.noNA$entrez))
View(venn)
mapIds(org.Mm.eg.db, keys=c("23871", "14460"), column="SYMBOL", keytype="ENTREZID", multiVals="first")



#====================================
res.sig0 <- read.table(file="data/Results_sig_0.txt",row.names = 1)
res.sig1 <- read.table(file="data/Results_sig_1.txt",row.names = 1)
res.sig3 <- read.table(file="data/Results_sig_3.txt",row.names = 1)

head(res.sig0)

colnames(res.sig0) <- paste(colnames(res.sig0),"0h", sep = "_")
colnames(res.sig1) <- paste(colnames(res.sig1),"1h", sep = "_")
colnames(res.sig3) <- paste(colnames(res.sig3),"3h", sep = "_")

res.sig.merge <- merge(res.sig0, res.sig1, by=0, all=TRUE)
rownames(res.sig.merge) <- res.sig.merge$Row.names

res.sig.merge <- merge(res.sig.merge, res.sig3, by=0, all=TRUE)

res.sig.merge <- res.sig.merge[,-2] %>% dplyr::rename(gene=Row.names)
View(res.sig.merge)

res.sig.allFC <- res.sig.merge %>% dplyr::select(gene,log2FoldChange_0h,log2FoldChange_1h,log2FoldChange_3h,padj_0h,padj_1h,padj_3h)

View(res.sig.allFC)
write.table(res.sig.allFC, sep="\t",file="data/sig-0-1-3h_allFC-p.txt", 
            row.names=TRUE,col.names=NA,quote=FALSE)
write.csv(res.sig.allFC, sep="\t",file="data/sig-0-1-3h_allFC-p.csv", 
            row.names=TRUE,col.names=NA,quote=FALSE)
#
res.sig.allFC.long <- res.sig.allFC %>% 
  tidyr::pivot_longer(cols=2:4,values_to = "log2FC")%>% 
  na.omit()
res.sig.allFC.long$absFC <- abs(res.sig.allFC.long$log2FC)
View(res.sig.allFC.long)

res.sig.allFC.long1 <- res.sig.allFC.long[order(res.sig.allFC.long[,'gene'],-res.sig.allFC.long[,'absFC']),] 
View(res.sig.allFC.long1)

res.sig.allFC.long1 <- res.sig.allFC.long1[!duplicated(res.sig.allFC.long1$gene),] 
dup1=res.sig.allFC.long1 %>% group_by(gene) %>% dplyr::filter(n() > 1) 
View(dup1)

res.sig.allFC.long1 <- res.sig.allFC.long1[order(-res.sig.allFC.long1[,'absFC']),]
res.sig.allFC.long.top <- res.sig.allFC.long1[c(1:50),]
View(res.sig.allFC.long.top)

res.sig.allFC.top <- res.sig.allFC %>% dplyr::filter(gene %in% c(res.sig.allFC.long.top$gene))

colnames(res.sig.allFC.top) <- c("gene","0h","1h","3h")
res.sig.allFC.top[is.na(res.sig.allFC.top)] <- 0

View(res.sig.allFC.top)

library("gplots")
m <- as.matrix(as.data.frame(lapply(res.sig.allFC.top[,-1], as.numeric),check.names=F))
View(m)
hclust.ave <- function(x) hclust(x, method="average")
par(oma=c(2,2,2,2))
heatmap.2(m,
          labRow = res.sig.allFC.top$gene,
          scale = "none", 
          col = bluered(100), 
          trace = "none", 
          density.info = "none",
          cexRow = 0.7,
          cexCol = 0.8,
          offsetRow = -0.2,
          offsetCol = 0,
          distfun = dist,
          #hclustfun = hclust,
          hclustfun=hclust.ave,
          Colv=FALSE,
          dendrogram='row',
          key=TRUE, keysize=0.75, key.title = F,key.xlab="log2 fold change",cex.lab=5.0, cex.axis=5.0,
          #key.par=list(mar=c(1,1,1,1)),
          lhei=c(2,12), lwid=c(2,2),  
          key.par = list(cex=0.6)
          )
dev.off()
?heatmap.2()

####################### 
# Select top 50 DEGs with lowest adj P value among all 3 timepoints
res.sig.allFC.long <- res.sig.allFC %>% 
  tidyr::pivot_longer(cols=5:7,values_to = "adjP") 
View(res.sig.allFC.long)
res.sig.allFC.long1 <- res.sig.allFC.long[-c(2:4)] %>% na.omit()
res.sig.allFC.long1 <- res.sig.allFC.long1[order(res.sig.allFC.long1[,'gene'],res.sig.allFC.long1[,'adjP']),] 
View(res.sig.allFC.long1)

res.sig.allFC.long1 <- res.sig.allFC.long1[!duplicated(res.sig.allFC.long1$gene),] 
dup1=res.sig.allFC.long1 %>% group_by(gene) %>% dplyr::filter(n() > 1) 
View(dup1)

res.sig.allFC.long1 <- res.sig.allFC.long1[order(res.sig.allFC.long1[,'adjP']),]
res.sig.allFC.long.top <- res.sig.allFC.long1[c(1:50),]
View(res.sig.allFC.long.top)

res.sig.allFC.top <- res.sig.allFC %>% dplyr::filter(gene %in% c(res.sig.allFC.long.top$gene))

View(res.sig.allFC.top)
##
sig.m <- res.sig.allFC.top %>% mutate(x1="0h") %>% mutate(x2="1h") %>% mutate(x3="3h") %>%
  mutate(Log10adj.P_0=-log10(padj_0h)) %>% 
  mutate(Log10adj.P_1=-log10(padj_1h)) %>%
  mutate(Log10adj.P_3=-log10(padj_3h))
View(sig.m)

p <- ggplot(sig.m, 
            aes(y = fct_reorder(gene, Log10adj.P_0))) + 
  geom_point(aes(color = log2FoldChange_0h, size = Log10adj.P_0, x=x1)) +
  geom_point(aes(color = log2FoldChange_1h, size = Log10adj.P_1, x=x2))+
  geom_point(aes(color = log2FoldChange_3h, size = Log10adj.P_3, x=x3))+
  theme_bw(base_size = 8) +
  theme(panel.grid.major = element_blank())+
  scale_color_gradient2(#limits=c(-8,5) ,
    midpoint = 0, low = "blue4", mid = "white",
    high = "red4", space = "Lab" )+
  ylab(NULL)+
  xlab(NULL)+
  xlim("0h","1h","3h")+
  theme(axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        legend.title = element_text(color = "black", size = 8))+
  labs(size = "-log10 adj. p",color="log2 fold change") # 4.5x7.5in
# legend.text = element_blank())+
# 
p
##
## rank by fold change then p
head(sig.m)
sig.m1 <- sig.m

sig.m1 <- sig.m1%>% mutate(rankFC=coalesce(log2FoldChange_0h,log2FoldChange_1h))
sig.m1 <- sig.m1%>% mutate(rankP=coalesce(Log10adj.P_0,Log10adj.P_1))
sig.m1.up <- sig.m1 %>% filter(rankFC>0)
View(sig.m1.up)

sig.m1.down <- sig.m1 %>% filter(rankFC<0)
View(sig.m1.down)


sig.m1.up <- sig.m1.up[order(-sig.m1.up[,'Log10adj.P_0'],-sig.m1.up[,'Log10adj.P_1'],-sig.m1.up[,'Log10adj.P_3']),] 
sig.m1.down <- sig.m1.down[order(-sig.m1.down[,'Log10adj.P_0'],-sig.m1.down[,'Log10adj.P_1'],-sig.m1.down[,'Log10adj.P_3']),] 
sig.m1.down=sig.m1.down[order(nrow(sig.m1.down):1),]

sig.m2 <- rbind(sig.m1.up,sig.m1.down)
View(sig.m2)
sig.m2$y <- c(1:50)
p <- ggplot(sig.m2, 
            aes(y = fct_reorder(gene,-y))) + 
  geom_point(aes(color = log2FoldChange_0h, size = Log10adj.P_0, x=x1)) +
  geom_point(aes(color = log2FoldChange_1h, size = Log10adj.P_1, x=x2))+
  geom_point(aes(color = log2FoldChange_3h, size = Log10adj.P_3, x=x3))+
  theme_bw(base_size = 8) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_color_gradient2(#limits=c(-8,5) ,
    midpoint = 0, low = "blue4", mid = "white",
    high = "red4", space = "Lab" )+
  ylab(NULL)+
  xlab(NULL)+
  xlim("0h","1h","3h")+
  theme(axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        legend.title = element_text(color = "black", size = 8))+
  labs(size = "-log10 adj. p",color="log2 fold change") # 2.5x7.5in
p
p <- ggplot(sig.m2, aes(y = gene)) + 
  geom_point(aes(color = log2FoldChange_0h, size = Log10adj.P_0, x=x1)) +
  geom_point(aes(color = log2FoldChange_1h, size = Log10adj.P_1, x=x2))+
  geom_point(aes(color = log2FoldChange_3h, size = Log10adj.P_3, x=x3))+
  theme_bw(base_size = 8) +
  theme(panel.grid.major = element_blank())+
  scale_color_gradient2(#limits=c(-8,5) ,
    midpoint = 0, low = "blue4", mid = "white",
    high = "red4", space = "Lab" )+
  ylab(NULL)+
  xlab(NULL)+
  xlim("0h","1h","3h")+
  theme(axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        legend.title = element_text(color = "black", size = 8))+
  labs(size = "-log10 adj. p",color="log2 fold change") # 4.5x7.5in

# legend.text = element_blank())+
# 
p
?factor()
##
##
library("gplots")
m <- as.matrix(as.data.frame(lapply(res.sig.allFC.top[,-1], as.numeric),check.names=F))
View(m)
hclust.ave <- function(x) hclust(x, method="average")
par(oma=c(2,2,2,2))
heatmap.2(m,
          labRow = res.sig.allFC.top$gene,
          scale = "none", 
          col = bluered(100), 
          trace = "none", 
          density.info = "none",
          cexRow = 0.7,
          cexCol = 0.8,
          offsetRow = -0.2,
          offsetCol = 0,
          distfun = dist,
          #hclustfun = hclust,
          hclustfun=hclust.ave,
          Colv=FALSE,
          dendrogram='row',
          key=TRUE, keysize=0.75, key.title = F,key.xlab="log2 fold change",cex.lab=5.0, cex.axis=5.0,
          #key.par=list(mar=c(1,1,1,1)),
          lhei=c(2,12), lwid=c(2,2),  
          key.par = list(cex=0.6)
)
dev.off()
?heatmap.2()
