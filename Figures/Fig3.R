library(tidyverse)
library(patchwork)





##### Fig3B: ZBTB16 expression in Validation Cohort

# Load ZBTB16 expression data
df = read.csv("Source_Data/Fig3B_ZBTB16_Validation_Cohort.csv")

# MRD calculations for plotting
df$MRD_D28 = as.numeric(df$MRD_D28)
df$MRD_D28_adj = df$MRD_D28 / 100 + 1e-05

# Add MRD_D28_cat
df$MRD_D28_cat = NA
df$MRD_D28_cat[df$MRD_D28 == 0] = "0"
df$MRD_D28_cat[df$MRD_D28 > 0] = "0-1"
df$MRD_D28_cat[df$MRD_D28 >= 1] = "1-5"
df$MRD_D28_cat[df$MRD_D28 >= 5] = ">=5"

# Add Response_cat
df$Response_cat = "Responsive_Remission"
df$Response_cat[df$relapse == "Yes"] = "Responsive_Relapse"
df$Response_cat[df$response == "nonResponsive"] = "nonResponsive_Remission"
df$Response_cat[(df$response == "nonResponsive") & (df$relapse == "Yes")] = "nonResponsive_Relapse"
df$Response_cat[df$patient_ID %in% c("L074")] = "Induction_death"
df$Response_cat[df$patient_ID %in% c("T030")] = "Relapse_only"

# Define clusters with sample_ID
df = df %>% mutate(sample_cluster = paste0(sample_ID, "::", leiden))

# Calculate number of cells per cluster
sample_cluster_TABLE = df %>% count(sample_cluster)
match_index = match(x = df$sample_cluster, table = sample_cluster_TABLE$sample_cluster)
df$sample_cluster_COUNT = sample_cluster_TABLE$n[match_index]

# Calculate number of blasts per sample
sample_ID_TABLE = df %>% count(sample_ID)
match_index = match(x = df$sample_ID, table = sample_ID_TABLE$sample_ID)
df$sample_ID_COUNT = sample_ID_TABLE$n[match_index]

# Calculate fraction of blasts in each cluster
df$sample_cluster_FRACT = df$sample_cluster_COUNT / df$sample_ID_COUNT

# Calculate ZBTB16 expression per cluster
df_cluster = df[c("patient_ID", "sample_ID", "timepoint", "Response_cat", "sample_cluster", "sample_cluster_FRACT", "ZBTB16")] %>%
  group_by(sample_cluster) %>%
  mutate(exp = median(ZBTB16)) %>% 
  ungroup() %>% 
  distinct(sample_cluster, .keep_all = TRUE)

# Define clusters which are ZBTB16 positive and negative
df_cluster$exp_cat = NA
df_cluster$exp_cat[df_cluster$exp == 0] = "Zneg"
df_cluster$exp_cat[df_cluster$exp > 0] = "Zpos"



# Generate order for patient_ID (with L026 coming with all other Responsive_Remission patients)
clinical_df = df %>% distinct(patient_ID, .keep_all = TRUE)
clinical_df$Response_cat_v2 = clinical_df$Response_cat
clinical_df$Response_cat_v2[clinical_df$patient_ID %in% c("L026")] = "L026"
clinical_df$Response_cat_v2 = factor(clinical_df$Response_cat_v2, levels =  c("Induction_death", "L026", "Responsive_Remission", "Responsive_Relapse", "nonResponsive_Remission", "nonResponsive_Relapse", "Relapse_only"))
order_patient_ID = clinical_df %>%
  arrange(Response_cat_v2, MRD_D28, patient_ID) %>%
  pull(patient_ID)

# Set order for x-axis (patient_ID)
clinical_df$patient_ID = factor(clinical_df$patient_ID, levels = order_patient_ID)
df_cluster$patient_ID = factor(df_cluster$patient_ID, levels = order_patient_ID)

# Set order for horizontal facets (Response_cat)
clinical_df$Response_cat = factor(clinical_df$Response_cat, levels = c("Induction_death", "Responsive_Remission", "Responsive_Relapse", "nonResponsive_Remission", "nonResponsive_Relapse", "Relapse_only"))
df_cluster$Response_cat = factor(df_cluster$Response_cat, levels = c("Induction_death", "Responsive_Remission", "Responsive_Relapse", "nonResponsive_Remission", "nonResponsive_Relapse", "Relapse_only"))



# Plot 1: Induction Day 28 MRD
p1 = ggplot(clinical_df, aes(x = patient_ID, y = MRD_D28_adj, col = MRD_D28_cat)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_point(size = 1) +
  facet_grid(cols = vars(Response_cat), scales = "free_x", space = "free_x", as.table = FALSE) +
  scale_y_continuous(
    trans = "log10", 
    breaks = c(1e-05, 1e-04, 1e-03, 1e-02, 1e-01, 1), 
    labels = c('0%', '0.01%', '0.1%', '1%', '10%', '100%'), 
    limits = c(1e-05, 1), 
    expand = c(0.1, 0)
  ) +
  scale_color_manual(
    values = c("0" = "#1D71B8", "0-1" = "#20A7DB", "1-5" = "#A0D9EF", ">=5" = "#F77F11"),
    breaks = c("0", "0-1", "1-5", ">=5"),
    labels = c("0%", "0-1%", "1-5%", ">=5%")
  ) +
  ylab("Day 28 MRD (%)") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

# Plot2: ZBTB16 expression
p2 = ggplot() +
  #geom_point(data = df_cluster, mapping = aes(x = patient_ID, y = exp, size = sample_cluster_FRACT, color = exp_cat), shape = 19, show.legend = FALSE) +
  #scale_size_area(max_size = 3) +
  geom_point(data = df_cluster, mapping = aes(x = patient_ID, y = exp, color = exp_cat), shape = 19, size = 1.7, show.legend = FALSE) +
  facet_grid(cols = vars(Response_cat), rows = vars(timepoint), scales = "free_x", space = "free_x", as.table = FALSE) +
  scale_color_manual(values = c("Zneg" = "grey70", "Zpos" = "grey30")) +
  scale_y_continuous(breaks = 0:3, labels = 0:3) +
  coord_cartesian(
    ylim = c(-0.2, 3)
  ) +
  ylab("ZBTB16") +
  theme_bw() +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

# Combine plots
p = p2 + p1 + plot_layout(nrow = 2, ncol = 1, heights = c(2, 0.7), guides = "keep") &
  theme()
ggsave("Plots/Fig3B_ZBTB16_Validation_Cohort.pdf", p, width = 6.8, height = 5.3)


