library(tidyverse)
library(patchwork)
library(ggbeeswarm)
library(ggrastr)





##### Fig1F: ZBTB16 expression in discovery cohort

df = read.csv("Source_Data/Fig1FG_ZBTB16_Discovery_Cohort.csv")

df = df %>% 
  mutate(sample_cluster = paste0(sample_ID, "::", leiden)) %>%
  mutate(timepoint_response = paste0(timepoint, "_", response))

df$timepoint_response = factor(df$timepoint_response, levels = c("D0_Responsive", "D0_nonResponsive", "D28_nonResponsive"))

sample_cluster_TABLE = df %>% count(sample_cluster)
match_index = match(x = df$sample_cluster, table = sample_cluster_TABLE$sample_cluster)
df$sample_cluster_COUNT = sample_cluster_TABLE$n[match_index]

sample_ID_TABLE = df %>% count(sample_ID)
match_index = match(x = df$sample_ID, table = sample_ID_TABLE$sample_ID)
df$sample_ID_COUNT = sample_ID_TABLE$n[match_index]

df$sample_cluster_FRACT = df$sample_cluster_COUNT / df$sample_ID_COUNT
  
df_sample = df %>%
  group_by(sample_ID) %>%
  mutate(exp = median(.data[["ZBTB16"]])) %>%
  mutate(mean = mean(.data[["ZBTB16"]])) %>% 
  ungroup() %>% 
  distinct(sample_ID, .keep_all = TRUE) %>%
  arrange(exp, mean)

df_cluster = df %>%
  group_by(sample_cluster) %>%
  mutate(exp = median(.data[["ZBTB16"]])) %>% 
  ungroup() %>% 
  distinct(sample_cluster, .keep_all = TRUE)

df_binary = df[c("sample_ID", "timepoint_response", "ZBTB16")] %>%
  mutate(is_positive = .data[["ZBTB16"]] > 0)

df_sample$sample_ID = factor(df_sample$sample_ID, levels = df_sample$sample_ID)
df$sample_ID = factor(df$sample_ID, levels = df_sample$sample_ID)
df_cluster$sample_ID = factor(df_cluster$sample_ID, levels = df_sample$sample_ID)
df_binary$sample_ID = factor(df_binary$sample_ID, levels = df_sample$sample_ID)

p1 = ggplot() + 
  geom_point(data = df_cluster, mapping = aes(x = sample_ID, y = exp, fill = timepoint_response, size = sample_cluster_FRACT), 
             stroke = NA, shape = 21, show.legend = FALSE) + 
  scale_size_area(max_size = 4) + 
  scale_fill_manual(values = c("D0_Responsive" = "#1D71B8", "D0_nonResponsive" = "#F77F11", "D28_nonResponsive" = "#D62728")) +
  facet_grid(cols = vars(timepoint_response), scales = "free_x", space = "free_x") + 
  ylab("Cluster expression") + 
  ggtitle("ZBTB16") + 
  theme_bw() + 
  theme(
    axis.text.x = element_blank()
  )

p2 = ggplot(data = df, mapping = aes(x = sample_ID, y = .data[["ZBTB16"]], col = timepoint_response)) +
  rasterise(geom_quasirandom(size = 0.2, col = "#D3D3D3", show.legend = FALSE), dpi = 1200) +
  geom_boxplot(alpha = 0, outliers = FALSE) +
  scale_colour_manual(values = c("D0_Responsive" = "#1D71B8", "D0_nonResponsive" = "#F77F11", "D28_nonResponsive" = "#D62728")) +
  facet_grid(cols = vars(timepoint_response), scales = "free_x", space = "free_x") +
  ylab("Blast expression") +
  theme_bw() +
  theme(
    strip.text.x = element_blank(), 
    axis.text.x = element_blank()
  )

p3 = ggplot(df_binary, aes(x = sample_ID, fill = is_positive)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = c("#FFFFFF", "#000000"), labels = c("Not expressed", "Expressed"), name = "") + 
  facet_grid(cols = vars(timepoint_response), scales = "free_x", space = "free_x") + 
  scale_y_continuous(breaks = c(0, 1), labels = c(0, 1), limits = c(0, 1), expand = c(0, 0)) +
  ylab("Positive\nblast fraction") + 
  theme_bw() +
  theme(
    strip.text.x = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9)
  )

p = p1 + p2 + p3 + plot_layout(nrow = 3, ncol = 1, heights = c(2, 1.3, 1), guides = "collect") & 
  theme(
    strip.background = element_blank(),
    legend.position = "none",
    legend.title = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("Plots/Fig1F_ZBTB16_Discovery_Cohort.pdf", p, width = 7, height = 6)





##### Fig1G: ZBTB16 Day 0 vs Day 28

list_of_sample_ID = c(
  'P058_TALL_D0_BM_GEX', 'P058_TALL_D28_BM_GEX', 
  'P028_TALL_D0_BM_GEX', 'P028_TALL_D28_BM_GEX', 
  'P016_TALL_D0_BM_GEX', 'P016_TALL_D28_BM_GEX', 
  'P029_TALL_D0_BM_GEX', 'P029_TALL_D28_BM_GEX', 
  'P021_TALL_D0_BM_GEX', 'P021_TALL_D28_BM_GEX', 
  'P084_TALL_D0_BM_GEX', 'P084_TALL_D28_BM_GEX', 
  'P030_TALL_D0_BM_GEX', 'P030_TALL_D28_BM_GEX', 
  'P018_TALL_D0_BM_GEX', 'P018_TALL_D28_BM_GEX'
)

df = read.csv("Source_Data/Fig1FG_ZBTB16_Discovery_Cohort.csv")
df = df[df$sample_ID %in% list_of_sample_ID, ]

df_sample = df %>%
  group_by(sample_ID) %>%
  mutate(is_positive = .data[["ZBTB16"]] > 0) %>%
  mutate(positive_fract = mean(is_positive)) %>%
  filter(is_positive) %>%
  mutate(median_positive = median(.data[["ZBTB16"]])) %>%
  ungroup() %>%
  distinct(sample_ID, .keep_all = TRUE)
df_sample$sample_ID = factor(df_sample$sample_ID, levels = list_of_sample_ID)
df_sample = df_sample %>% arrange(sample_ID)
df_sample = df_sample[c('sample_ID', 'timepoint', 'positive_fract', 'median_positive')]

p1 = ggplot(df_sample, aes(x = sample_ID, y = median_positive, fill = timepoint)) +
  geom_col() +
  scale_fill_manual(values = c("D0" = "#F77F11", "D28" = "#D62728")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  ylab("Median") +
  theme_bw() +
  theme(
    axis.text.x = element_blank()
  )

p2 = ggplot(df_sample, aes(x = sample_ID, y = positive_fract, fill = timepoint)) +
  geom_col() +
  scale_fill_manual(values = c("D0" = "#F77F11", "D28" = "#D62728")) +
  scale_y_continuous(breaks = c(0, 1), labels = c(0, 1), limits = c(0, 1), expand = c(0, 0)) +
  ylab("Fraction") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

p = p1 + p2 + plot_layout(nrow = 2, ncol = 1, heights = c(1, 1), guides = "collect") &
  theme(
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("Plots/Fig1G_ZBTB16_D0_D28.pdf", p, width = 4, height = 4)







