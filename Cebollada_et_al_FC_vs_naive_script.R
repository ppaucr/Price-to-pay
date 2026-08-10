## Transcriptomics analysis for Paula Cebollada-Rica
## Interest in cell death
## Marta E. Camarena
## 22th April 2026

library(dplyr)
library(stringr)
library(tidyr)
library(pheatmap)
library(UpSetR)
library(org.Mm.eg.db)
library(GO.db)
library(clusterProfiler)
library(msigdbr)
library(ggplot2)
library(ggpubr)

## read tambe of expression, converted to csv
data=read.csv("full path to normalized expression table")
outdir="outdir_to_save_plots"

names(data) = gsub("^X","",names(data))
long_data = data %>% 
  pivot_longer(cols=-c("gene_id","gene_name", "gene_type"), names_to = "mouse", values_to = "CPM")
long_data$group = gsub("..*d","d",long_data$mouse)

## make sure there are no big differences per timepoint between both individuals
### correlations
corr_data <- long_data %>%
  dplyr::select(-mouse) %>%
  group_by(gene_id, gene_name, gene_type, group) %>%
  mutate(replicate = row_number()) %>%   # create replicate 1 and 2
  ungroup() %>%
  pivot_wider(
    names_from = replicate,
    values_from = CPM,
    names_prefix = "rep"
  )
replicate_cor <- corr_data %>%
  mutate(
    rep1 = log2(rep1 + 1),
    rep2 = log2(rep2 + 1)
  ) %>%
  group_by(group) %>%
  summarise(
    pearson_cor = cor(rep1, rep2, method = "pearson", use = "complete.obs"),
    spearman_cor = cor(rep1, rep2, method = "spearman", use = "complete.obs"),
    .groups = "drop"
  )
replicate_cor ## > 0.9

## compute differences
differences_df = long_data %>%
  dplyr::group_by(gene_id, gene_name, gene_type, group) %>%
  dplyr::mutate(
    log2FC_rep = abs(log2((CPM[1] + 1) / (CPM[2] + 1))),
    diff_CPM = abs(CPM[1] - CPM[2]),
    mean_logCPM = mean(log(CPM+1)))

# differences_df %>%
#   dplyr::select(gene_name, gene_id, group, log2FC_rep) %>%
#   unique() %>%
#   ggplot(aes(x=log2FC_rep)) +
#   scale_x_continuous(limits = c(0,5)) +
#   geom_vline(xintercept =1, color="red") +
#   labs(x="abs(rep1 log2CPM/rep2 log2CPM)") +
#   geom_density() +
#   labs(title="Absolute difference of expression per gene and replicate") +
#   facet_wrap(~group, nrow=2)
# ggsave(file.path(outdir,differences_per_gene.jpg"), height=2.5, width=5.85)

corr_data %>%
  mutate(
    rep1 = log2(rep1 + 1),
    rep2 = log2(rep2 + 1)
  ) %>%
  ggplot(aes(x=rep1, y=rep2)) +
  geom_point(size=.3) +
  labs(title="Correlation of gene expression per replicate",
       y="Replicate 2",
       x="Replicate 1") +
  facet_wrap(~ group, nrow=2) +
  stat_cor(method = "pearson", size=2, label.sep = "\n", digits=4)
ggsave(file.path(outdir,"Correlation_GeneExpression.jpg"), height=5.17, width=6.16)

# average value of logCPM
clean_average_data <- differences_df %>%
  group_by(gene_id) %>%
  # filter(all(abs(log2FC_rep) <= 1)) %>%   # gene kept only if ALL timepoints pass
  ungroup() %>%
  dplyr::select(gene_id, gene_name, gene_type, group, mean_logCPM) %>%
  distinct()

## compute log2FC vs d0
logFC_vs_d0 = clean_average_data %>%
  ungroup() %>%
  group_by(gene_id, gene_name, gene_type) %>%
  mutate(
    d0_expr = mean_logCPM[group == "d0"][1],
    log2FC_vs_d0 = mean_logCPM - d0_expr
  ) %>%
  ungroup()

## control 
d0 = logFC_vs_d0 %>% filter(group == "d0")
summary(d0$log2FC_vs_d0)

# ggplot(logFC_vs_d0, aes(x=log2FC_vs_d0)) +
#   geom_density() +
#   geom_vline(xintercept = 1) +
#   geom_vline(xintercept = -1) +
#   facet_wrap(~ group)
# dev.off()

############# log2FC > 1.5
log2FC_1point5_genes = logFC_vs_d0 %>%
  filter(log2FC_vs_d0 > 1.5) %>%
  subset(group != "d31") %>%
  pull(gene_name)
log2FC_1point5 = logFC_vs_d0 %>%
  filter(!(gene_name == "Gbp11" & gene_type != "protein_coding")) %>%
  dplyr::select(-c(mean_logCPM, d0_expr, gene_id, gene_type)) %>% unique() %>%
  filter(gene_name %in% log2FC_1point5_genes) %>%
  # distinct(gene_name, group, .keep_all = TRUE) %>%   # keep first duplicate
  pivot_wider(names_from = group, values_from = log2FC_vs_d0)

log2FC_1point5_matrix = as.matrix(log2FC_1point5)
rownames(log2FC_1point5_matrix) = log2FC_1point5$gene_name
log2FC_1point5_matrix = log2FC_1point5_matrix[,-1]
class(log2FC_1point5_matrix) <- "numeric"

renaming = log2FC_1point5_matrix
colnames(renaming) = c("Naive","d1 p.i.","d2 p.i.","d3 p.i.","d4 p.i.","d5 p.i.","d6 p.i.","d7 p.i.","d9 p.i.","d31 p.i.")

## heatmap of logFC
pdf(file.path(outdir,"heatmap_log2FC_1point5.main.pdf"), height=2.76, width=3.90)
ph = pheatmap(renaming, scale="none", 
              cluster_cols = F,
              treeheight_row = 0,
              show_rownames = F,
              border_color = NA,
              color = colorRampPalette(c("#F5F2D8","#58C9BC","#833D98"))(50),
              main="log2FC of those with log2FC>1.5",
              fontsize_row = 4,
              fontsize = 6,
              angle_col = 90
)
dev.off()
write.csv(renaming, file.path(outdir,"heatmap_log2FC_1point5.csv"))

## keep the rows order for showing the expression data
order_rows = ph$tree_row$order

################# complementary data #################
# ## display the CPM of those with log2FC > 1.5
# log2FC_1point5_CPMs = logFC_vs_d0 %>% 
#   filter(!(gene_name == "Gbp11" & gene_type != "protein_coding")) %>%
#   dplyr::select(-c(log2FC_vs_d0, d0_expr, gene_id, gene_type)) %>% unique() %>%
#   filter(gene_name %in% log2FC_1point5_genes) %>%
#   pivot_wider(names_from = group, values_from = mean_logCPM)
# 
# log2FC_1point5_CPMs_matrix = as.matrix(log2FC_1point5_CPMs)
# rownames(log2FC_1point5_CPMs_matrix) = log2FC_1point5_CPMs$gene_name
# log2FC_1point5_CPMs_matrix = log2FC_1point5_CPMs_matrix[,-1]
# class(log2FC_1point5_CPMs_matrix) <- "numeric"
# 
# ## order by the previous plot
# log2FC_1point5_CPMs_matrix <- log2FC_1point5_CPMs_matrix[order_rows,]
# pdf(file.path(outdir,"heatmap_log2FC_1point5.CPMs.pdf"), height=5.08, width=3.01)
# pheatmap(log2FC_1point5_CPMs_matrix, scale="none",
#          cluster_cols = F,
#          cluster_rows = F,
#          show_rownames = F,
#          border_color = NA,
#          color = colorRampPalette(c("#E2E6BD","#8E063B"))(50),
#          main="log2CPM of those with log2FC > 1.5",
#          fontsize_row = 4
# )
# dev.off()
####################################################################
## are these genes involved in pathways of cell death?
celldeath_paths <- AnnotationDbi::select(org.Mm.eg.db, keys=c("GO:0006915","GO:0097300","GO:0043068"), columns = c('SYMBOL','ENSEMBL'), keytype = "GOALL")
genes_4_enrichment <- merge(logFC_vs_d0 %>% select(gene_name, gene_id),log2FC_1point5 %>% select(gene_name))
genes_4_enrichment = merge(genes_4_enrichment, celldeath_paths %>% select(-EVIDENCEALL) %>% unique(), by.x="gene_id", by.y="ENSEMBL", all.x=T)
genes_4_enrichment = unique(genes_4_enrichment)
table(genes_4_enrichment$GOALL)
print(genes_4_enrichment %>% drop_na())

## paths of interest
pathName = data.frame("GOALL" = c("GO:0006915","GO:0097300","GO:0043068"),
                      "PATHWAY" = c("Apoptotic process", "Programmed necrotic cell death", "Positive regulation of programmed cell death"))
write.csv(merge(pathName,genes_4_enrichment %>% drop_na() %>% select(-c(ONTOLOGYALL,SYMBOL))), file.path(outdir,"genes_log2FC_1point5_celldeathpathways.csv"), row.names = F)

## enrichment of these genes - M5 subcollection GO: Gene Ontology
mouse_go_df <- msigdbr(species = "mouse", db_species = "MM", collection  = "M5", subcollection = "GO:BP")
term2gene_mapping <- mouse_go_df[, c("gs_name", "gene_symbol")]
# term2name_mapping <- mouse_go_df[, c("gs_name", "gene_symbol")]

enricher_results <- enricher(genes_4_enrichment$gene_name,
                             pAdjustMethod = "fdr",
                             universe = data$gene_name,
                             minGSSize = 20,
                             TERM2GENE = term2gene_mapping)
# write.csv(enricher_results, file.path(outdir,"enricher_log2FC_1point5.GO_BP_20GENES.csv"))

## plot enrichment
enricher_df <- as.data.frame(enricher_results)
enricher_df$GeneRatio_num <- as.numeric(sub("/\\d+", "", enricher_df$GeneRatio)) / as.numeric(sub("\\d+/", "", enricher_df$GeneRatio))
enricher_df$BgRatio_num <- as.numeric(sub("/\\d+", "", enricher_df$BgRatio)) / as.numeric(sub("\\d+/", "", enricher_df$BgRatio))
enricher_df$negLogFDR <- -log10(enricher_df$p.adjust)
enricher_df <- enricher_df %>%
  mutate(
    Description = Description %>%
      str_remove("^GOBP_") %>%      # remove prefix
      str_to_lower() %>%            # all lowercase
      str_replace_all("_", " ") %>% # replace underscores with spaces
      str_to_sentence()            # capitalize first letter
  ) 
## subset 1%
fdr1percent = enricher_df %>% subset(negLogFDR > 2)
# write.csv(fdr1percent, file.path(outdir,"enricherResults_log2FC_1point5.20GENESFDR001.csv"))

# we make a shorter version of a super long GO term
# Adaptive immune response via somatic recombination stands for regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains
enricher_df$Description = gsub("Adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains","Adaptive immune response via somatic recombination*", enricher_df$Description)

## initial plot
ggplot(enricher_df, aes(x = FoldEnrichment, y = reorder(Description, FoldEnrichment))) +
  geom_point(aes(size = Count, color = negLogFDR)) +
  scale_size_continuous(range = c(1,4), name = "Gene Count",  breaks = c(4,8,12,16)
  ) +
  scale_color_gradient(low = "red", high = "green", name = "-log10(FDR)") +
  guides(
    colour = guide_colourbar(order = 1, barwidth = 5, barheight = 0.5),
    size = guide_legend(order = 2, nrow=1)
  ) +
  labs(x = "Fold Enrichment", y = NULL)+
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8, color = "black"),
    axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
    axis.title.x = element_text(size = 10, color = "black", face = 'bold'),
    plot.title = element_text(hjust = 0.5, size = 10, color = "black", face = 'bold'),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    axis.text.y.right = element_text(hjust = 0),
    legend.text = element_text(size=10),
    legend.title = element_text(size = 10),
    legend.key.size = unit(0.4, 'cm'),
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  scale_y_discrete(position = "right") 
dev.off()


## Count > 15
## FDR1%
fdr1percent_count15 <- fdr1percent %>%
  filter(Count > 15)
fdr1percent_count15$Description = gsub("Adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains","Adaptive immune response via somatic recombination*", fdr1percent_count15$Description)
ggplot(fdr1percent_count15, aes(x = FoldEnrichment, y = reorder(Description, FoldEnrichment))) +
  geom_point(aes(size = Count, color = negLogFDR)) +
  scale_size_continuous(range = c(1,4), name = "Gene Count") +
  scale_color_gradient(low = "red", high = "green", name = "-log10(FDR)") +
  guides(
    colour = guide_colourbar(order = 1, barwidth = 5, barheight = 0.5),
    size = guide_legend(order = 2, nrow=1)
  ) +
  labs(x = "Fold Enrichment", y = NULL)+ #, title = "enricher (Gene Count > 15)") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8, color = "black"),
    axis.text.x = element_text(size = 11, color = "black", angle = 45, hjust = 1),
    axis.title.x = element_text(size = 10, color = "black", face = 'bold'),
    plot.title = element_text(hjust = 0.5, size = 10, color = "black", face = 'bold'),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    axis.text.y.right = element_text(hjust = 0),
    legend.text = element_text(size=10),
    legend.title = element_text(size = 10),
    legend.key.size = unit(0.4, 'cm'),
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  scale_y_discrete(
    position = "right",
    expand = expansion(mult = c(0, 0))
  )
ggsave(file.path(outdir,"enricher_plot_count15_FDR001.png"), dpi=500,height=10.47, width=8)