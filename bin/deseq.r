#!/usr/bin/env Rscript

################################################
################################################
## REQUIREMENTS                               ##
################################################
################################################

## THIS IS PROGRAM IS PROVIDED BY ZDQ,YOU CAN###
## COMMUNICATE ANY PROBLEM WITH 
## zhangdongqin2@126.com

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(DESeq2)
library(BiocParallel)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationHub)
library(ggrepel)
library(amap)
library(clusterProfiler)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(
    make_option(c("-i", "--count_file"    ), type="character", default=NULL    , metavar="path"   , help="Count file matrix where rows are genes and columns are samples."                        ),
    make_option(c("-f", "--count_col"     ), type="integer"  , default=2       , metavar="integer", help="First column containing sample count data."                                             ),
    make_option(c("-d", "--id_col"        ), type="integer"  , default=1       , metavar="integer", help="Column containing identifiers to be used."                                              ),
    make_option(c("-r", "--sample_suffix" ), type="character", default=''      , metavar="string" , help="Suffix to remove after sample name in columns e.g. '.rmDup.bam' if 'DRUG_R1.rmDup.bam'."),
    make_option(c("-o", "--outdir"        ), type="character", default='./'    , metavar="path"   , help="Output directory."                                                                      ),
    make_option(c("-p", "--outprefix"     ), type="character", default='deseq2', metavar="string" , help="Output prefix."                                                                         ),
    make_option(c("-v", "--vst"           ), type="logical"  , default=FALSE   , metavar="boolean", help="Run vst transform instead of rlog."                                                     ),
    make_option(c("-c", "--cores"         ), type="integer"  , default=1       , metavar="integer", help="Number of cores."                                                                       ),
    make_option(c("-e", "--coldata_file"  ), type="character", default=NULL    , metavar="path"   , help="Coldata file which contain columns are samples and condition."                          ),
    make_option(c("-b", "--sample_con"    ), type="character", default=''      , metavar="string" , help="the comparison names of DEG analysis. eg.'tumor-normal'"								  )
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$count_file)){
    print_help(opt_parser)
    stop("Please provide a counts file.", call.=FALSE)
}
if (is.null(opt$coldata_file)){
    print_help(opt_parser)
    stop("Please provide a coldata file.", call.=FALSE)
}

################################################
################################################
## READ IN COUNTS FILE                        ##
################################################
################################################

count.table           <- read.table(file=opt$count_file,header=TRUE,sep="\t",row.names=1)
col.table		      <- read.table(file=opt$coldata_file,header=TRUE,sep="\t")
#rownames(count.table) <- count.table[,opt$id_col]
#count.table           <- count.table[,opt$count_col:ncol(count.table),drop=FALSE]
#colnames(count.table) <- gsub(opt$sample_suffix,"",colnames(count.table))
#colnames(count.table) <- gsub(pattern='\\.$', replacement='', colnames(count.table))

################################################
################################################
## RUN DESEQ2                                 ##
################################################
################################################

if (file.exists(opt$outdir) == FALSE) {
    dir.create(opt$outdir,recursive=TRUE)
}
setwd(opt$outdir)

#samples.vec <- sort(colnames(count.table))
#groups      <- sub("_[^_]+$", "", samples.vec)
#groups      <- col.table$condition
#if (length(unique(groups)) == 1 || length(unique(groups)) == length(samples.vec)) {
#    quit(save = "no", status = 0, runLast = FALSE)
#}

DDSFile <- paste(opt$outprefix,".dds.RData",sep="")
if (file.exists(DDSFile) == FALSE) {
    counts  <- as.matrix(count.table)
    #counts  <- count.table[,samples.vec,drop=FALSE]
    coldata <- col.table
    #coldata <- data.frame(row.names=colnames(counts), condition=groups)
    dds     <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~ condition)
    dds     <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(opt$cores))
    if (!opt$vst) {
        vst_name <- "rlog"
        rld      <- rlog(dds)
    } else {
        vst_name <- "vst"
        rld      <- varianceStabilizingTransformation(dds)
    }
    assay(dds, vst_name) <- assay(rld)
    save(dds,file=DDSFile)
} else {
    load(DDSFile)
    vst_name <- intersect(assayNames(dds), c("vst", "rlog"))
    if (length(vst_name)==0) { # legacy might mean vst was saved as a separate object called rld
        vst_name <- "loaded_rld"
        assay(dds, vst_name) <- assay(rld)
    } else {
        vst_name==vst_name[1]
    }
}

################################################
################################################
## PLOT QC                                    ##
################################################
################################################

##' PCA pre-processeor
##'
##' Generate all the necessary information to plot PCA from a DESeq2 object
##' in which an assay containing a variance-stabilised matrix of counts is
##' stored. Copied from DESeq2::plotPCA, but with additional ability to
##' say which assay to run the PCA on, and adds an assessment of how well
##' each PC explains the experimental grouping of the data.
##' 
##' @param object The DESeq2DataSet object.
##' @param intgroup interesting groups: a character vector of names in 'colData(x)' to use for grouping.
##' @param ntop number of top genes to use for principla components, selected by highest row variance.
##' @param assay the name or index of the assay that stores the variance-stabilised data.
##' @return A data.frame containing the projected data alongside the grouping columns.
##' A 'percentVar' attribute is set which includes the percentage of variation each PC explains,
##' and additionally how much the variation within that PC is explained by the grouping variable.
##' @author Gavin Kelly
plotPCA_vst <- function (object, intgroup = "condition", ntop = 500, assay=length(assays(object))) {
    rv         <- rowVars(assay(object, assay))
    select     <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca        <- prcomp(t(assay(object, assay)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = ":"))
    }  else {
        colData(object)[[intgroup]]
    }
    d <- cbind(pca$x, group = group, intgroup.df, name = colnames(object))
    percentFrame <- data.frame(PC=seq(along=percentVar), percentVar=100*percentVar, groupR=0.0)
    for (ipc in seq(along=percentVar)) {
        fit1 <- lm(pca$x[,ipc]  ~ group)
        percentFrame$groupR[ipc] <- 100*summary(fit1)$r.squared
    }
    attr(d, "percentVar") <- percentFrame
    return(d)
}

PlotFile <- paste(opt$outprefix,".plots.pdf",sep="")
if (file.exists(PlotFile) == FALSE) {
    pdf(file=PlotFile,onefile=TRUE,width=7,height=7)

    ## PCA
    ntop <- c(500, Inf)
    for (n_top_var in ntop) {
        pca.data      <- plotPCA_vst(dds, assay=vst_name,intgroup=c("condition"),ntop=n_top_var)
        percentVar    <- round(attr(pca.data, "percentVar")$percentVar)
        plot_subtitle <- ifelse(n_top_var==Inf, "All genes", paste("Top", n_top_var, "genes"))
        pl <- ggplot(pca.data, aes(PC1, PC2, color=condition)) +
              geom_point(size=3) +
              xlab(paste0("PC1: ",percentVar[1],"% variance")) +
              ylab(paste0("PC2: ",percentVar[2],"% variance")) +
              labs(title = paste0("First PCs on ", vst_name, "-transformed data"), subtitle = plot_subtitle) + 
              theme(legend.position="top",
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=1))
        print(pl)

        pl <- ggplot(attr(pca.data, "percentVar"), aes(x=PC, y=percentVar)) +
              geom_line(aes(colour="explained by PC")) +
              geom_line(aes(y=groupR, colour="of PC explained by condition")) +
              scale_x_continuous(breaks=seq(along=percentVar), minor_breaks=NULL)  +
              labs(title="Diagnostics of PCs", subtitle=plot_subtitle, x="Component", y="Percentage explaned", colour="Percentage variation") +
              theme_bw() +
              theme(legend.position="top")
        print(pl)

        pc_r <- order(attr(pca.data, "percentVar")$groupR, decreasing=TRUE)
        pl <- ggplot(pca.data, aes_string(paste0("PC", pc_r[1]), paste0("PC", pc_r[2]), color="condition")) +
              geom_point(size=3) +
              xlab(paste0("PC", pc_r[1], ": ",percentVar[pc_r[1]],"% variance")) +
              ylab(paste0("PC", pc_r[2], ": ",percentVar[pc_r[2]],"% variance")) +
              labs(title = paste0("Group-Explanatory PCs of ", vst_name, "-tranformed data"), subtitle = plot_subtitle) + 
              theme(legend.position="top",
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=1))
        print(pl)
    } # at end of loop, we'll be using the user-defined ntop if any, else all genes
    
    ## WRITE PC1 vs PC2 VALUES TO FILE
    pca.vals           <- pca.data[,1:2]
    colnames(pca.vals) <- paste0(colnames(pca.vals), ": ", percentVar[1:2], '% variance')
    pca.vals           <- cbind(sample = rownames(pca.vals), pca.vals)
    write.table(pca.vals,file=paste(opt$outprefix,".pca.vals.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=TRUE)

    ## SAMPLE CORRELATION HEATMAP
    sampleDists      <- dist(t(assay(dds, vst_name)))
    sampleDistMatrix <- as.matrix(sampleDists)
    colors           <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pheatmap(
        sampleDistMatrix,
        clustering_distance_rows=sampleDists,
        clustering_distance_cols=sampleDists,
        col=colors,
        main=paste("Euclidean distance between", vst_name, "of samples")
    )

    ## WRITE SAMPLE DISTANCES TO FILE
    write.table(cbind(sample = rownames(sampleDistMatrix), sampleDistMatrix),file=paste(opt$outprefix,".sample.dists.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

    dev.off()
}

################################################
################################################
## SAVE SIZE FACTORS                          ##
################################################
################################################

SizeFactorsDir <- "size_factors/"
if (file.exists(SizeFactorsDir) == FALSE) {
    dir.create(SizeFactorsDir,recursive=TRUE)
}

NormFactorsFile <- paste(SizeFactorsDir,opt$outprefix,".size_factors.RData",sep="")
if (file.exists(NormFactorsFile) == FALSE) {
    normFactors <- sizeFactors(dds)
    save(normFactors,file=NormFactorsFile)

    for (name in names(sizeFactors(dds))) {
        sizeFactorFile <- paste(SizeFactorsDir,name,".txt",sep="")
        if (file.exists(sizeFactorFile) == FALSE) {
            write(as.numeric(sizeFactors(dds)[name]),file=sizeFactorFile)
        }
    }
}

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

RLogFile <- "R_sessionInfo.log"
if (file.exists(RLogFile) == FALSE) {
    sink(RLogFile)
    a <- sessionInfo()
    print(a)
    sink()
}

################################################
################################################
################################################
################################################
################################################
################################################
## DESEQ2 RESULTS OUTPUT                      ##
################################################
################################################

dds <- estimateSizeFactors(dds) #对dds进行参数估计
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

normalized_counts <- counts(dds, normalized=TRUE) ##提取counts
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
write.table(normalized_counts, file="zdq_trans.Count_matrix.xls.DESeq2.normalized.xls",
quote=F, sep="\t", row.names=T, col.names=T)
rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
write.table(rlogMat, file="zdq_trans.Count_matrix.xls.DESeq2.normalized.rlog.xls",
quote=F, sep="\t", row.names=T, col.names=T)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
hc <- hcluster(t(rlogMat), method="pearson")
pdf("zdq_trans.Count_matrix.xls.DESeq2.normalized.rlog.pearson.pdf", pointsize=10)
heatmap(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none",
col=hmcol, margins=c(11,11), main="The pearson correlation of each
sample")
dev.off()

sampleA = unlist(strsplit(opt$sample_con, "-"))[1] 
sampleB = unlist(strsplit(opt$sample_con, "-"))[2]
print (opt$sample_con)
print (sampleA)
print (sampleB)



contrastV <- c("condition", sampleA, sampleB)
print (contrastV)


res <- results(dds)
baseA <- counts(dds, normalized=TRUE)[, colData(dds)$condition == sampleA ]

if (is.vector(baseA)){
    baseMeanA <- as.data.frame(baseA)
} else {
    baseMeanA <- as.data.frame(rowMeans(baseA))
}
colnames(baseMeanA) <- sampleA
#head(baseMeanA)
baseB <- counts(dds, normalized=TRUE)[, colData(dds)$condition == sampleB ]
if (is.vector(baseB)){
        baseMeanB <- as.data.frame(baseB)
} else {
        baseMeanB <- as.data.frame(rowMeans(baseB))
}
colnames(baseMeanB) <- sampleB
#head(baseMeanB)
res <- cbind(baseMeanA, baseMeanB, as.data.frame(res))
res <- cbind(ID=rownames(res), as.data.frame(res))
res$baseMean <- rowMeans(cbind(baseA, baseB))
res$padj[is.na(res$padj)] <- 1
res <- res[order(res$pvalue),]
comp314 <- paste(sampleA, "_vs_", sampleB, sep=".")
# 生成文件名
file_base <- paste("zdq_trans.Count_matrix.xls.DESeq2", comp314, sep=".")
file_base1 <- paste(file_base, "results.xls", sep=".")

write.table(as.data.frame(res), file=file_base1, sep="\t", quote=F, row.names=F)

res_de <- subset(res, res$padj<0.1, select=c('ID', sampleA,
        sampleB, 'log2FoldChange', 'padj'))
res_de_up <- subset(res_de, res_de$log2FoldChange>=1)

file <- paste("zdq_trans.Count_matrix.xls.DESeq2",sampleA,"_UP_",sampleB,
        'xls', sep=".") 
write.table(as.data.frame(res_de_up), file=file, sep="\t", quote=F, row.names=F)
res_de_dw <- subset(res_de, res_de$log2FoldChange<=(-1)*1)

file <- paste("zdq_trans.Count_matrix.xls.DESeq2",sampleA, "_DOWN_", sampleB, 
        'xls', sep=".") 
write.table(as.data.frame(res_de_dw), file=file, sep="\t", quote=F, row.names=F)
# 差异基因ID
res_de_up_id = data.frame(ID=res_de_up$ID, 
        type=paste(sampleA,"_higherThan_", sampleB, sep="."))
res_de_dw_id = data.frame(ID=res_de_dw$ID, 
        type=paste(sampleA,"_lowerThan_", sampleB, sep="."))
de_id = rbind(res_de_up_id, res_de_dw_id)
file <- "zdq_trans.Count_matrix.xls.DESeq2.all.DE"
write.table(as.data.frame(de_id), file=file, sep="\t", quote=F, row.names=F, col.names=F)

res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.01),'sig'] <- 'up'
res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.01),'sig'] <- 'down'
res1[which(abs(res1$log2FoldChange) <= 1 | res1$padj >= 0.01),'sig'] <- 'none'

pdf("zdq_trans.DESeq2.volcano.pdf")
ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
geom_point(size = 1) +  #绘制散点图
scale_color_manual(values = c('red', 'gray', 'green'), limits = c('up', 'none', 'down')) +  #自定义点的颜色
labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = 'control vs treat', color = '') +  #坐标轴标题
theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), #背景色、网格线、图例等主题修改
    panel.background = element_rect(color = 'black', fill = 'transparent'), 
    legend.key = element_rect(fill = 'transparent')) +
geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +  #添加阈值线
geom_hline(yintercept = 2, lty = 3, color = 'black') +
xlim(-12, 12) + ylim(0, 35)
dev.off()
deg_expr <- normalized_counts[rownames(normalized_counts) %in% de_id$ID,]

pdf("zdq_trans.DESeq2.DEG.heatmap.pdf")
pheatmap(deg_expr, cluster_row=T, scale="row")
dev.off()
################################################
################################################
################################################
################################################
################################################
################################################
## GO AND KEGG ENRICHMENT ANALYSIS            ##
################################################
################################################
gene_symbol=as.vector(de_id$ID)
up_symbol=as.vector(res_de_up_id$ID)
down_symbol=as.vector(res_de_dw_id$ID)


gene_id= mapIds(org.Hs.eg.db,keys=gene_symbol,column="ENTREZID",keytype="SYMBOL",multiVals="first")
GO.CC <- enrichGO(OrgDb="org.Hs.eg.db", gene = gene_id, ont = "CC", pvalueCutoff = 1, readable= TRUE)
pdf("zdq.trans.go.CC.circle.pdf")
dotplot(GO.CC,showCategory=10,title="Enrichment GO.CC Top10") #泡泡图
dev.off()
pdf("zdq.trans.go.CC.barplot.pdf")
barplot(GO.CC, showCategory=20,title="EnrichmentGO.CC")  #柱状图
dev.off()


GO.MF <- enrichGO(OrgDb="org.Hs.eg.db", gene = gene_id, ont = "MF", pvalueCutoff = 1, readable= TRUE)
pdf("zdq.trans.go.MF.circle.pdf")
dotplot(GO.MF,showCategory=10,title="Enrichment GO.MF Top10") #泡泡图
dev.off()
pdf("zdq.trans.go.MF.barplot.pdf")
barplot(GO.MF, showCategory=20,title="EnrichmentGO.MF")  #柱状图
dev.off()


GO.BP <- enrichGO(OrgDb="org.Hs.eg.db", gene = gene_id, ont = "BP", pvalueCutoff = 1, readable= TRUE)
pdf("zdq.trans.go.BP.circle.pdf")
dotplot(GO.BP,showCategory=10,title="Enrichment GO.BP Top10") #泡泡图
dev.off()
pdf("zdq.trans.go.BP.barplot.pdf")
barplot(GO.BP, showCategory=20,title="EnrichmentGO.BP")  #柱状图
dev.off()


write.csv(as.data.frame(GO.BP), 'DEGS_go.BP.csv', sep = '\t', row.names = FALSE, quote = FALSE)
write.csv(as.data.frame(GO.CC), 'DEGS_go.CC.csv', sep = '\t', row.names = FALSE, quote = FALSE)
write.csv(as.data.frame(GO.MF), 'DEGS_go.MF.csv', sep = '\t', row.names = FALSE, quote = FALSE)


KEGG <- enrichKEGG(gene= gene_id,organism  = 'hsa', qvalueCutoff = 0.05)
pdf("DEGS_KEGG.pdf")
dotplot(KEGG,font.size=8)
dev.off()
write.csv(as.data.frame(KEGG), 'DEGS_KEGG.csv', sep = '\t', row.names = FALSE, quote = FALSE)
################################################
################################################
################################################
################################################
################################################
################################################
## UPGENES GO AND KEGG ENRICHMENT ANALYSIS    ##
################################################
################################################

up_id= mapIds(org.Hs.eg.db,keys=up_symbol,column="ENTREZID",keytype="SYMBOL",multiVals="first")
upGO.CC <- enrichGO(OrgDb="org.Hs.eg.db", gene = up_id, ont = "CC", pvalueCutoff = 1, readable= TRUE)
pdf("up.go.CC.circle.pdf")
dotplot(upGO.CC,showCategory=10,title="Enrichment GO.CC Top10") #泡泡图
dev.off()
pdf("up.go.CC.barplot.pdf")
barplot(upGO.CC, showCategory=20,title="EnrichmentGO.CC")  #柱状图
dev.off()


upGO.MF <- enrichGO(OrgDb="org.Hs.eg.db", gene = up_id, ont = "MF", pvalueCutoff = 1, readable= TRUE)
pdf("up.go.MF.circle.pdf")
dotplot(upGO.MF,showCategory=10,title="Enrichment GO.MF Top10") #泡泡图
dev.off()
pdf("up.go.MF.barplot.pdf")
barplot(upGO.MF, showCategory=20,title="EnrichmentGO.MF")  #柱状图
dev.off()


upGO.BP <- enrichGO(OrgDb="org.Hs.eg.db", gene = up_id, ont = "BP", pvalueCutoff = 1, readable= TRUE)
pdf("up.go.BP.circle.pdf")
dotplot(upGO.BP,showCategory=10,title="Enrichment GO.BP Top10") #泡泡图
dev.off()
pdf("up.go.BP.barplot.pdf")
barplot(upGO.BP, showCategory=20,title="EnrichmentGO.BP")  #柱状图
dev.off()


write.csv(as.data.frame(upGO.BP), 'upDEGS_go.BP.csv', sep = '\t', row.names = FALSE, quote = FALSE)
write.csv(as.data.frame(upGO.CC), 'upDEGS_go.CC.csv', sep = '\t', row.names = FALSE, quote = FALSE)
write.csv(as.data.frame(upGO.MF), 'upDEGS_go.MF.csv', sep = '\t', row.names = FALSE, quote = FALSE)


upKEGG <- enrichKEGG(gene= up_id,organism  = 'hsa', qvalueCutoff = 0.05)
pdf("upDEGS_KEGG.pdf")
dotplot(upKEGG,font.size=8)
dev.off()
write.csv(as.data.frame(upKEGG), 'upDEGS_KEGG.csv', sep = '\t', row.names = FALSE, quote = FALSE)

################################################
################################################
################################################
################################################
################################################
################################################
## DOWNGENES GO AND KEGG ENRICHMENT ANALYSIS  ##
################################################
################################################

down_id = mapIds(org.Hs.eg.db,keys=down_symbol,column="ENTREZID",keytype="SYMBOL",multiVals="first")
downGO.CC <- enrichGO(OrgDb="org.Hs.eg.db", gene = down_id, ont = "CC", pvalueCutoff = 1, readable= TRUE)
pdf("down.go.CC.circle.pdf")
dotplot(downGO.CC,showCategory=10,title="Enrichment GO.CC Top10") #泡泡图
dev.off()
pdf("down.go.CC.barplot.pdf")
barplot(downGO.CC, showCategory=20,title="EnrichmentGO.CC")  #柱状图
dev.off()


downGO.MF <- enrichGO(OrgDb="org.Hs.eg.db", gene = down_id, ont = "MF", pvalueCutoff = 1, readable= TRUE)
pdf("down.go.MF.circle.pdf")
dotplot(downGO.MF,showCategory=10,title="Enrichment GO.MF Top10") #泡泡图
dev.off()
pdf("down.go.MF.barplot.pdf")
barplot(downGO.MF, showCategory=20,title="EnrichmentGO.MF")  #柱状图
dev.off()


downGO.BP <- enrichGO(OrgDb="org.Hs.eg.db", gene = down_id, ont = "BP", pvalueCutoff = 1, readable= TRUE)
pdf("down.go.BP.circle.pdf")
dotplot(downGO.BP,showCategory=10,title="Enrichment GO.BP Top10") #泡泡图
dev.off()
pdf("down.go.BP.barplot.pdf")
barplot(downGO.BP, showCategory=20,title="EnrichmentGO.BP")  #柱状图
dev.off()


write.csv(as.data.frame(downGO.BP), 'downDEGS_go.BP.csv', sep = '\t', row.names = FALSE, quote = FALSE)
write.csv(as.data.frame(downGO.CC), 'downDEGS_go.CC.csv', sep = '\t', row.names = FALSE, quote = FALSE)
write.csv(as.data.frame(downGO.MF), 'downDEGS_go.MF.csv', sep = '\t', row.names = FALSE, quote = FALSE)


downKEGG <- enrichKEGG(gene= down_id,organism  = 'hsa', qvalueCutoff = 0.05)
pdf("downDEGS_KEGG.pdf")
dotplot(downKEGG,font.size=8)
dev.off()
write.csv(as.data.frame(downKEGG), 'downDEGS_KEGG.csv', sep = '\t', row.names = FALSE, quote = FALSE)
