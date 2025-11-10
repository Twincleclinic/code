library(tidyverse)
library(VennDiagram)

aspartame <- read.csv("Aspartame_targets.csv", stringsAsFactors = FALSE)
mafld <- read.csv("MAFLD_targets.csv", stringsAsFactors = FALSE)

aspartame_genes <- unique(aspartame$Gene)
mafld_genes <- unique(mafld$Gene)

intersect_genes <- intersect(aspartame_genes, mafld_genes)
write.csv(data.frame(Gene = intersect_genes), "Intersection_genes.csv", row.names = FALSE)

venn.plot <- draw.pairwise.venn(
  area1 = length(aspartame_genes),
  area2 = length(mafld_genes),
  cross.area = length(intersect_genes),
  category = c("Aspartame", "MAFLD"),
  fill = c("#f4b266", "#b8643c"),
  cat.cex = 1.2,
  fontface = "bold",
  euler.d = FALSE,
  scaled = FALSE
)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(forcats)

A <- read.csv("A_intersection_genes.csv", stringsAsFactors = FALSE)
B <- read.csv("B_first_screen_genes.csv", stringsAsFactors = FALSE)
C <- read.csv("C_second_screen_genes.csv", stringsAsFactors = FALSE)
D <- read.csv("D_second_screen_genes.csv", stringsAsFactors = FALSE)

A_entrez <- bitr(A$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
B_entrez <- bitr(B$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
C_entrez <- bitr(C$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
D_entrez <- bitr(D$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ego_A <- enrichGO(gene=A_entrez$ENTREZID, OrgDb=org.Hs.eg.db, ont="ALL", pvalueCutoff=0.05, qvalueCutoff=0.05, readable=TRUE)
ekegg_A <- enrichKEGG(gene=A_entrez$ENTREZID, organism="hsa", pvalueCutoff=0.05)
df_GO_A <- as.data.frame(ego_A)
df_KEGG_A <- as.data.frame(ekegg_A)
df_GO_A$Category <- factor(df_GO_A$ONTOLOGY, levels=c("BP","CC","MF"))
df_GO_A <- df_GO_A %>% group_by(ONTOLOGY) %>% slice_min(order_by=p.adjust, n=10)
df_KEGG_A <- df_KEGG_A %>% slice_min(order_by=p.adjust, n=10)
ggplot(df_GO_A, aes(x=reorder(Description, -log10(p.adjust)), y=-log10(p.adjust), fill=Category)) +
  geom_bar(stat="identity") + coord_flip() + theme_bw() + 
  scale_fill_manual(values=c("#FF7F50","#FFD700","#87CEEB")) + 
  labs(x=NULL, y="-log10(p.adjust)")
ggsave("A_GO_barplot.pdf", width=8, height=6)

ggplot(df_KEGG_A, aes(x=reorder(Description, -log10(p.adjust)), y=-log10(p.adjust))) +
  geom_bar(stat="identity", fill="#8FBC8F") + coord_flip() + theme_bw() + 
  labs(x=NULL, y="-log10(p.adjust)")
ggsave("A_KEGG_barplot.pdf", width=8, height=6)

ego_B <- enrichGO(gene=B_entrez$ENTREZID, OrgDb=org.Hs.eg.db, ont="ALL", pvalueCutoff=0.05, qvalueCutoff=0.05, readable=TRUE)
ekegg_B <- enrichKEGG(gene=B_entrez$ENTREZID, organism="hsa", pvalueCutoff=0.05)
df_B_GO <- as.data.frame(ego_B) %>% separate_rows(geneID, sep="/")
df_B_KEGG <- as.data.frame(ekegg_B) %>% separate_rows(geneID, sep="/")
df_B_all <- bind_rows(df_B_GO %>% mutate(Type="GO"), df_B_KEGG %>% mutate(Type="KEGG"))
ggplot(df_B_all, aes(axis1=geneID, axis2=Description, y=Count)) +
  geom_alluvium(aes(fill=p.adjust)) +
  geom_stratum(width=1/8, fill="gray90", color="gray") +
  geom_text(stat="stratum", aes(label=after_stat(stratum)), size=2.3) +
  theme_void() + scale_fill_gradient(low="#91c9ff", high="#ff0000")
ggsave("B_GO_KEGG_sankey.pdf", width=9, height=5)

ego_C <- enrichGO(gene=C_entrez$ENTREZID, OrgDb=org.Hs.eg.db, ont="ALL", pvalueCutoff=0.05, qvalueCutoff=0.05, readable=TRUE)
df_C <- as.data.frame(ego_C) %>% group_by(ONTOLOGY) %>% slice_min(order_by=pvalue, n=8)
ggplot(df_C, aes(x=GeneRatio, y=fct_reorder(Description, GeneRatio))) +
  geom_point(aes(size=Count, color=pvalue)) +
  facet_wrap(~ONTOLOGY, scales="free_y") +
  scale_color_gradient(low="blue", high="red") + theme_bw() +
  labs(x="GeneRatio", y=NULL)
ggsave("C_GO_bubble.pdf", width=8, height=6)

ekegg_D <- enrichKEGG(gene=D_entrez$ENTREZID, organism="hsa", pvalueCutoff=0.05)
df_D <- as.data.frame(ekegg_D) %>% slice_min(order_by=pvalue, n=15)
ggplot(df_D, aes(x=GeneRatio, y=fct_reorder(Description, GeneRatio))) +
  geom_point(aes(size=Count, color=pvalue)) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw() + labs(x="GeneRatio", y=NULL)
ggsave("D_KEGG_bubble.pdf", width=8, height=6)


library(DESeq2)
library(ggplot2)
library(ggsci)

counts <- read.csv("counts_matrix.csv", row.names = 1)
meta <- read.csv("metadata.csv", row.names = 1)
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = meta,
                              design = ~ Group)
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$group <- ifelse(res_df$log2FoldChange > 1 & res_df$padj < 0.05, "Up",
                       ifelse(res_df$log2FoldChange < -1 & res_df$padj < 0.05, "Down", "NotSig"))

ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=group)) +
  geom_point(alpha=0.8, size=1.2) +
  scale_color_manual(values=c("Down"="#4DBBD5FF","NotSig"="gray80","Up"="#E64B35FF")) +
  theme_bw(base_size=12) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="gray40") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="gray40") +
  labs(x="log2 Fold Change", y="-log10 adjusted P-value", color="Regulation") +
  theme(legend.position="top")
ggsave("A_DESeq2_volcano.pdf", width=6, height=5)

library(ggpubr)
tpm <- read.csv("TPM_data.csv")
genes <- c("ACE","ANPEP","MMP2","PTGS2")
plist <- list()
for (g in genes) {
  p <- ggviolin(tpm, x="Group", y=g, fill="Group",
                palette=c("#E64B35","#4DBBD5"),
                add="boxplot", add.params=list(fill="white")) +
    stat_compare_means(method="wilcox.test", label="p.format", size=3) +
    labs(title=g, y="TPM (log2+1)") + theme_bw(base_size=12)
  plist[[g]] <- p
}
pdf("B_TPM_violin.pdf", width=8, height=4)
ggarrange(plist$ACE, plist$ANPEP, plist$MMP2, plist$PTGS2, ncol=4)
dev.off()

library(caret)
library(gbm)
set.seed(123)

data <- read.csv("model_training.csv")
data$Group <- factor(data$Group, levels=c("Control","MAFLD"))
ctrl <- trainControl(method="cv", number=5, classProbs=TRUE, summaryFunction=twoClassSummary)
fit_gbm <- train(Group ~ ., data=data, method="gbm", trControl=ctrl,
                 metric="ROC", verbose=FALSE)
imp <- varImp(fit_gbm)$importance
imp$Feature <- rownames(imp)
imp <- imp[order(-imp$Overall), ][1:8, ]

ggplot(imp, aes(x=reorder(Feature, Overall), y=Overall)) +
  geom_bar(stat="identity", fill="#4DBBD5") +
  coord_flip() + theme_bw(base_size=12) +
  labs(x=NULL, y="Feature Importance (GBM)") +
  theme(panel.grid.minor=element_blank())
ggsave("C_GBM_importance.pdf", width=6, height=4)

library(pROC)
library(DecisionCurve)

prob <- predict(fit_gbm, data, type="prob")[,2]
roc_gbm <- roc(response=data$Group, predictor=prob, levels=c("Control","MAFLD"))
pdf("D_ROC_GBM.pdf", width=5, height=5)
plot(roc_gbm, col="#E64B35", lwd=2, legacy.axes=TRUE)
dev.off()

dca_data <- data.frame(Group=as.numeric(data$Group=="MAFLD"), pred=prob)
dca_gbm <- decision_curve(Group ~ pred, data=dca_data, family=binomial(link="logit"),
                          thresholds=seq(0,1,by=0.01), bootstraps=100)
pdf("E_DCA_curve.pdf", width=6, height=5)
plot_decision_curve(list(dca_gbm), curve.names=c("GBM"), col=c("#E64B35"))
dev.off()

roc_anpep <- roc(response=data$Group, predictor=data$ANPEP, levels=c("Control","MAFLD"))
pdf("F_ROC_ANPEP.pdf", width=5, height=5)
plot(roc_anpep, col="#E64B35", lwd=2, legacy.axes=TRUE)
dev.off()
auc(roc_anpep)

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(ggplot2)

exposure_dat <- fread("GCST90468313_exposure.txt")
exposure_dat <- format_data(
  exposure_dat,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  samplesize_col = "n"
)
exposure_dat <- exposure_dat %>% filter(pval.exposure < 5e-8)
exposure_dat <- clump_data(exposure_dat, clump_r2 = 0.01, clump_kb = 10000)

outcome_dat <- fread("GCST90275050_outcome.txt")
outcome_dat <- format_data(
  outcome_dat,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  samplesize_col = "n"
)

dat <- harmonise_data(exposure_dat, outcome_dat)

mr_results <- mr(dat, method_list = c(
  "mr_ivw",
  "mr_egger_regression",
  "mr_weighted_median",
  "mr_simple_mode",
  "mr_weighted_mode"
))
write.csv(mr_results, "MR_results.csv", row.names = FALSE)

heterogeneity <- mr_heterogeneity(dat)
pleiotropy <- mr_pleiotropy_test(dat)
leaveoneout <- mr_leaveoneout(dat)
write.csv(heterogeneity, "MR_heterogeneity.csv")
write.csv(pleiotropy, "MR_pleiotropy.csv")

pdf("A_forest.pdf", width=6, height=3)
mr_forest_plot(mr_singlesnp(dat))
dev.off()

pdf("B_scatter.pdf", width=5, height=5)
mr_scatter_plot(mr_results, dat)
dev.off()

pdf("C_leaveoneout.pdf", width=6, height=6)
mr_leaveoneout_plot(leaveoneout)
dev.off()

pdf("D_funnel.pdf", width=6, height=6)
mr_funnel_plot(mr_singlesnp(dat))
dev.off()

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(viridis)
library(RColorBrewer)

dirs <- list.dirs("data/", full.names = TRUE, recursive = FALSE)
objects <- lapply(dirs, function(x) {
  obj <- Read10X(data.dir = x)
  obj <- CreateSeuratObject(counts = obj, project = basename(x), min.cells = 3, min.features = 200)
  return(obj)
})
sce <- merge(objects[[1]], y = objects[-1], add.cell.ids = basename(dirs), project = "MAFLD")
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
sce <- subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 3000)
sce <- ScaleData(sce)
sce <- RunPCA(sce, npcs = 50)
sce <- RunHarmony(sce, group.by.vars = "orig.ident")
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:30)
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:30)
sce <- FindClusters(sce, resolution = 0.5)

DimPlot(sce, reduction = "umap", group.by = "orig.ident")
ggsave("Harmony_batchcheck.pdf", width=6, height=5)
DimPlot(sce, reduction = "umap", label = TRUE, group.by = "seurat_clusters")
ggsave("Harmony_cluster.pdf", width=6, height=5)

markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(sce, features = top10$gene, group.by = "seurat_clusters") + scale_fill_viridis()
ggsave("Harmony_marker_heatmap.pdf", width=8, height=6)

FeaturePlot(sce, features = c("CD3D","NCAM1","S100A9","CD79A","IGHG1","ANPEP"), reduction = "umap", cols = c("lightblue","red"))
ggsave("Harmony_featureplot.pdf", width=8, height=6)

sce$celltype <- Idents(sce)
prop <- table(sce$celltype, sce$orig.ident)
prop_df <- as.data.frame(prop)
ggplot(prop_df, aes(x=Var2, y=Freq, fill=Var1)) + geom_bar(stat="identity", position="fill") +
  scale_fill_brewer(palette="Set2") + theme_bw()
ggsave("Harmony_cell_composition.pdf", width=6, height=5)

freq <- prop.table(prop, margin=2)
freq_log2 <- log2(freq + 1e-6)
pheatmap(freq_log2, cluster_rows = TRUE, cluster_cols = TRUE, color = colorRampPalette(c("white","orange","red"))(100))
ggsave("Harmony_celltype_heatmap.pdf", width=5, height=4)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)
library(RColorBrewer)

sce <- readRDS("sce_harmony.rds")

Idents(sce) <- "celltype"
subset_obj <- subset(sce, idents = c("Myeloid"))
subset_obj <- NormalizeData(subset_obj)
subset_obj <- FindVariableFeatures(subset_obj, selection.method = "vst", nfeatures = 3000)
subset_obj <- ScaleData(subset_obj)
subset_obj <- RunPCA(subset_obj, npcs = 50)
subset_obj <- FindNeighbors(subset_obj, dims = 1:30)
subset_obj <- FindClusters(subset_obj, resolution = 0.5)
subset_obj <- RunUMAP(subset_obj, dims = 1:30)

DimPlot(subset_obj, reduction = "umap", label = TRUE, group.by = "seurat_clusters")
ggsave("A_subcluster_umap.pdf", width=6, height=5)

DimPlot(subset_obj, reduction = "umap", split.by = "orig.ident", group.by = "seurat_clusters")
ggsave("B_subcluster_umap_split.pdf", width=8, height=4)

markers <- FindAllMarkers(subset_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(markers, "Myeloid_subcluster_markers.csv", row.names = FALSE)

markers_select <- c("S100A8","S100A9","CD1C","FCER1A","CD5L","MARCO","CLEC4C","XCR1")
DotPlot(subset_obj, features = markers_select) + scale_color_gradient(low = "blue", high = "red") + theme_bw()
ggsave("C_subcluster_dotplot.pdf", width=8, height=4)

FeaturePlot(subset_obj, features = markers_select, cols = c("lightblue","red"))
ggsave("D_subcluster_featureplot.pdf", width=8, height=6)
library(Seurat)
library(monocle3)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(patchwork)

sce <- readRDS("Myeloid_subcluster.rds")

FeaturePlot(sce, features = "ANPEP", reduction = "umap", cols = c("black","purple"))
ggsave("A_ANPEP_umap.pdf", width=6, height=5)

VlnPlot(sce, features = "ANPEP", group.by = "celltype", pt.size = 0)
ggsave("B_ANPEP_violin.pdf", width=6, height=4)

cds <- as.cell_data_set(sce)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
ggsave("C_pseudotime.pdf", width=5, height=5)
plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster = FALSE)
ggsave("D_pseudotime_celltype.pdf", width=5, height=5)
plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name == "ANPEP",], min_expr = 0.1)
ggsave("E_ANPEP_pseudotime.pdf", width=5, height=4)

sce$ANPEP_group <- ifelse(sce@assays$RNA@data["ANPEP",] > median(sce@assays$RNA@data["ANPEP",]), "High", "Low")
DimPlot(sce, reduction="umap", group.by="ANPEP_group", cols=c("#999999","#E41A1C"))
ggsave("F_ANPEP_group_umap.pdf", width=6, height=5)

Idents(sce) <- "celltype"
DimPlot(sce, reduction="umap", group.by="ANPEP_group", split.by="celltype", cols=c("#999999","#E41A1C"))
ggsave("G_ANPEP_split1.pdf", width=6, height=5)
DimPlot(sce, reduction="umap", group.by="ANPEP_group", split.by="orig.ident", cols=c("#999999","#E41A1C"))
ggsave("H_ANPEP_split2.pdf", width=6, height=5)

Idents(sce) <- "ANPEP_group"
deg <- FindMarkers(sce, ident.1="High", ident.2="Low", logfc.threshold=0.25, min.pct=0.1)
deg <- na.omit(deg)
deg <- deg[deg$p_val_adj < 0.05, ]
write.csv(deg, "DEG_ANPEP_High_vs_Low.csv")

gene_list <- rownames(deg)
gene_id <- bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
kk <- enrichKEGG(gene = gene_id$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
dotplot(kk, showCategory=10)
ggsave("I_KEGG_ANPEP_Low.pdf", width=6, height=5)

up_genes <- rownames(deg[deg$avg_log2FC > 0, ])
up_id <- bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
kk_up <- enrichKEGG(gene = up_id$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
dotplot(kk_up, showCategory=10)
ggsave("J_KEGG_ANPEP_High.pdf", width=6, height=5)
