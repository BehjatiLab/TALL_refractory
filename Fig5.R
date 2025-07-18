library(tidyverse)
library(ggrastr)





##### Fig5C: P058 WGS BAF

list_of_sample_ID = c('P058_D0_WGS_BAF', 'P058_D28_WGS_BAF')
list_of_chrom = c('chr9', 'chr17')
list_of_breakpoint = c(36706115, 57793447)

for (i in seq_along(list_of_sample_ID)) {
  
  for (j in seq_along(list_of_chrom)) {
    
    sample_ID = list_of_sample_ID[i]
    chrom = list_of_chrom[j]
    breakpoint = list_of_breakpoint[j]
    
    df = read.csv(paste0("Source_Data/Fig5C_", sample_ID, ".csv"))
    df = df[df$CHR == chrom, ]
    
    p = ggplot(df, aes(x = POS, y = BAF)) +
      rasterise(geom_point(aes(color = MajorAllele), size = 0.1), dpi = 1200) +
      geom_hline(yintercept = 0.5, linetype = "dashed") +
      geom_vline(xintercept = breakpoint, linetype = "dashed") +
      scale_x_continuous(breaks = c(breakpoint), labels = round(c(breakpoint)/1e6, digits = 1)) +
      scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0.0", "0.5", "1.0"), limits = c(0, 1)) +
      scale_color_manual(values = c("TRUE" = "#F77F11", "FALSE" = "#000000"), na.value = "#D3D3D3") +
      ylab("BAF") +
      xlab("Genomic position (Mb)") +
      theme_bw() +
      theme(
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
    ggsave(paste0("Plots/Fig5C_", sample_ID, "_", chrom, ".pdf"), p, width = 4, height = 2.6)
  }
}






##### Fig5D: P058 scRNA BAF

list_of_sample_ID = c('P058_D0_scRNA_BAF', 'P058_D28_scRNA_BAF')
list_of_chrom = c('chr9', 'chr17')
list_of_breakpoint = c(36706115, 57793447)

for (i in seq_along(list_of_sample_ID)) {
  
  for (j in seq_along(list_of_chrom)) {
    
    sample_ID = list_of_sample_ID[i]
    chrom = list_of_chrom[j]
    breakpoint = list_of_breakpoint[j]
    
    df = read.csv(paste0("Source_Data/Fig5D_", sample_ID, ".csv"))
    df = df[df$CHR == chrom, ]
    
    p = ggplot(df, aes(x = POS, y = BAF)) +
      geom_point(aes(color = MajorAllele), size = 1) +
      geom_hline(yintercept = 0.5, linetype = "dashed") +
      geom_vline(xintercept = breakpoint, linetype = "dashed") +
      scale_x_continuous(breaks = c(breakpoint), labels = round(c(breakpoint)/1e6, digits = 1)) +
      scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0.0", "0.5", "1.0"), limits = c(0, 1)) +
      scale_color_manual(values = c("TRUE" = "#F77F11", "FALSE" = "#000000"), na.value = "#D3D3D3") +
      ylab("BAF") +
      xlab("Genomic position (Mb)") +
      theme_bw() +
      theme(
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
    ggsave(paste0("Plots/Fig5D_", sample_ID, "_", chrom, ".pdf"), p, width = 4, height = 2.6)
  }
}




