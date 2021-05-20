## step3--数据标准化，用于后续聚类分析
rld <- rlog(dds,blind = FALSE)
vsd <- vst(dds, blind=FALSE)
##a = rld
##b = vsd
##assay(a) <- limma::removeBatchEffect(assay(a),c(RNAseq_SampleAnno$batch))
##assay(b) <- limma::removeBatchEffect(assay(b),c(RNAseq_SampleAnno$batch))
save(rld,a,file = "rlog.RData")
save(vsd,b,file = "vsd.RData")
DESeq2::plotPCA(rld)
##DESeq2::plotPCA(a)
DESeq2::plotPCA(vsd)
##DESeq2::plotPCA(b)

### step4--用rld和vst数据做PCA和MDS 聚类
# PCA1,最终用的是rld算法
RNAseq_SampleAnno$condition = factor(RNAseq_SampleAnno$condition,levels = c("pES","6dSSEA4pos","6dSSEA4neg",
                                                                            "8dAPRnegOCT4neg","10dAPRnegOCT4neg","12dAPRnegOCT4neg","14dAPRnegOCT4neg",
                                                                            "8dAPRposOCT4neg","10dAPRposOCT4neg","10dAPRposOCT4pos",
                                                                            "12dAPRposOCT4neg","12dAPRposOCT4pos",
                                                                           "14dAPRposOCT4neg","14dAPRposOCT4pos","nES"))

pca_result_x <- plotPCA(rld,returnData=T)
ggplot(pca_result_x,aes(x=PC1,y=PC2,color = RNAseq_SampleAnno$condition)) + 
  geom_point(cex=7,alpha = 1) + 
  theme_bw() + theme(panel.grid=element_blank()) + 
  xlab("PC1: 71% variance") + ylab("PC2: 13% variance") +
  scale_color_manual(values = rev(colorRampPalette(brewer.pal(11,"Spectral")[2:11])(15)))
  
# PCA2,未采用
pca_result <- as.data.frame(assay(vsd)) %>% t() %>% prcomp(center = TRUE, scale. = TRUE)
pca_result_x <- as.data.frame(pca_result$x)
ggplot(pca_result_x,aes(x=PC1,y=PC2,color = RNAseq_SampleAnno$condition)) + 
  geom_point(cex=6) + 
  geom_text_repel(aes(label=RNAseq_SampleAnno$condition),size = 3.5) + 
  theme_bw()
eigenvals <- (pca_result$sdev)^2
PC1var_explained <- eigenvals[1] / sum(eigenvals)
PC2var_explained <- eigenvals[2] / sum(eigenvals)

###MDS聚类
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rownames(sample_anno)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#距离聚类
?pheatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
plot(hclust(sampleDists))
#MDS聚类
library(stats)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = condition)) +
  geom_point(size = 6) + coord_fixed() + ggtitle("MDS with vsd data") + 
  geom_text_repel(aes(label=RNAseq_SampleAnno$condition),size = 3.5) +
  xlab("M1") + ylab("M2") + theme_bw()

### step5 对标准化的数据做heatmap聚类分析
library(pheatmap)
a = assay(rld)
pheatmap(cor(a),
         border_color = NA,
         cellwidth = 12,
         cellheight = 12)

pheatmap(cor(a),
         border_color = NA,
         cellwidth = 12,
         cellheight = 12,
         color = rev(colorRampPalette(c(brewer.pal(11,"Spectral")[1:11]))(100)))

