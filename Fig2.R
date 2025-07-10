library(tidyverse)
library(circlize)
library(ComplexHeatmap)





##### Figure 2B: Normal-to-leukaemia transcriptome comparison by logistic regression

# Load data
similarity_matrix = read.csv(file = paste0("Data/Fig2B_Logistic_Regression.csv"), row.names = 1, check.names = FALSE)
similarity_matrix = as.matrix(similarity_matrix)

# Ordering of rows and columns
celltypes_order = c('ZBTB16_neg', 'ZBTB16_pos', 'SCP')
cell_groupings_order = c('HSC_MPP', 'CYCLING_MPP', 'LMPP_MLP', 'DN(early)_T', 'DN(P)_T', 'DN(Q)_T', 'DP(P)_T', 'DP(Q)_T', 'ABT(ENTRY)', 'CD8+T', 'CD4+T', 'TREG', 'CYCLING_T', 'CD8AA', 'TYPE_1_INNATE_T', 'TYPE_3_INNATE_T', 'ILC2', 'ILC3', 'CYCLING_ILC', 'NK', 'CYCLING_NK', 'SCP')

# Extract cell groupings
cell_groupings = rownames(similarity_matrix) %>% sub(pattern = "::.*", replacement = "")

# Convert to discrete values
mtx = similarity_matrix
mtx_discrete = matrix(nrow = nrow(mtx), ncol = ncol(mtx))
mtx_discrete[mtx < 0.3] = "a"
mtx_discrete[mtx >= 0.3 & mtx < 0.4] = "b"
mtx_discrete[mtx >= 0.4 & mtx <= 0.6] = "c"
mtx_discrete[mtx > 0.6 & mtx <= 0.7] = "d"
mtx_discrete[mtx > 0.7] = "e"
rm(mtx)

# Add column names and convert mtx_discrete into a data frame
colnames(mtx_discrete) = colnames(similarity_matrix)
mtx_discrete = as.data.frame(mtx_discrete)

# Add cell_groupings
mtx_discrete$cell_groupings = cell_groupings

# Pivot longer
mtx_discrete = mtx_discrete %>%
  pivot_longer(-c(cell_groupings), names_to = "ref_celltypes", values_to = "similarity")

# Set order for similarity, ref_celltypes and cell_groupings
mtx_discrete = mtx_discrete %>%
  mutate(similarity = factor(similarity, levels = c("e", "d", "c", "b", "a")),
         ref_celltypes = factor(ref_celltypes, levels = celltypes_order),
         cell_groupings = factor(cell_groupings, levels = cell_groupings_order))

# Plot stacked barplot
p = ggplot(mtx_discrete, aes(y = 0, fill = similarity)) +
  geom_bar(stat = "count", position = "fill", color = "black", width = 0.1) +
  scale_x_reverse() +
  scale_fill_manual(values = c("a" = "#7F7F7F", "b" = "#B3B3B3", "c" = "#FFFFFF", "d" = "#77BDAD", "e" = "#309E88"),
                    labels = c("a"= "0.0 - <0.3", "b" = "0.3 - <0.4", "c" = "0.4 - 0.6", "d" = ">0.6 - 0.7", "e" = ">0.7 - 1.0")) +
  facet_grid(cell_groupings ~ ref_celltypes, switch = "x") +
  theme_void() +
  theme(strip.text.x = element_text(margin = margin(5, 5, 5, 5), angle = 90, vjust = 0.5, hjust = 1, size = 10),
        strip.text.y = element_text(margin = margin(5, 5, 5, 5), vjust = 0.5, hjust = 0, size = 10),
        legend.position = "bottom",
        plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"))
ggsave("Plots/Fig2B_Logisitc_Regression.pdf", p, height = 8, width = 8)





##### Fig2C: Heatmap of marker genes in T-ALL blasts

conv_genes = c("DNTT", "RAG1", "RAG2", "CD5", "CD27", "CD1A", "CD8A", "CD8B")
unconv_genes = c("GNG4", "ZNF683", "RORA", "KLRB1","EOMES", "KIT", "TBX21", "ZBTB16")

# Generate matrix
df = read.csv("Data/Fig2C_TALL_Marker_Genes.csv", check.names = FALSE)
df$category = paste0(df$timepoint, "_", df$response, "_", df$ZBTB16_status)
df = df %>%
  group_by(category) %>%
  summarise(across(all_of(c(conv_genes, unconv_genes)), mean)) %>%
  ungroup() %>%
  mutate(category = factor(category, levels = c('D0_Responsive_ZBTB16neg', 'D0_nonResponsive_ZBTB16neg', 'D28_nonResponsive_ZBTB16neg', 'D0_Responsive_ZBTB16pos', 'D0_nonResponsive_ZBTB16pos', 'D28_nonResponsive_ZBTB16pos'))) %>%
  arrange(category)
mat = as.matrix(df %>% select(-c('category')))
rownames(mat) = df$category
mat = apply(mat, 2, function(x) (x - min(x)) / (max(x) - min(x)))
mat = t(mat)

# Generate column annotation
column_df = df['category'] 
column_df[c('V1', 'V2', 'V3')] = str_split_fixed(string = column_df$category, pattern = "_", n = 3)
column_df$V1_V2 = paste0(column_df$V1, "_", column_df$V2)
column_ha = HeatmapAnnotation(
  which = "column",
  "Group" = column_df$V1_V2,
  col = list('Group' = c('D0_Responsive' = '#1D71B8', 'D0_nonResponsive' = '#F77F11', 'D28_nonResponsive' = '#D62728')),
  gp = gpar(col = "grey50"),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10),
  simple_anno_size = unit(3, "mm"),
  annotation_legend_param = list(border = "grey50")
)

# Generate row annotation
row_df = data.frame(genes = c(conv_genes, unconv_genes))
row_df$groups = "unsure"
row_df$groups[row_df$genes %in% conv_genes] = "conv"
row_df$groups[row_df$genes %in% unconv_genes] = "unconv"

# Plot heatmap
col_fun = colorRamp2(c(0, 1), c("#FFFFFF", "#000000"))
hm = Heatmap(
  mat,
  col = col_fun,
  name = "Scaled Mean Expression",
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_df$groups,
  row_gap = unit(2, "mm"),
  row_title = NULL,
  column_split = column_df$V3,
  column_gap = unit(2, "mm"),
  column_title = NULL,
  row_names_side = "left",
  column_names_side = "bottom",
  rect_gp = gpar(col = "grey50"),
  top_annotation = column_ha,
  bottom_annotation = column_ha,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(direction = "vertical", at = c(0, 1), border = "grey50")
)
pdf("Plots/Fig2C_TALL_Marker_Genes.pdf", width = 4, height = 6.2)
hm = draw(
  hm,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE,
)
dev.off()




