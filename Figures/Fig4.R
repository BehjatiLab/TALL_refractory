library(tidyverse)
library(ggrepel)
library(ggbeeswarm)
library(survival)
library(survminer)
library(pammtools)
library(patchwork)





##### Fig4A: Princess Maxima Center (PMC) cohort - plotting against day 28 MRD categories

df = read.csv("Source_Data/Fig4A_PMC.csv")
df$D28_MRD_category = factor(df$D28_MRD_category, levels = c("0", "0-1", "1-5", ">=5"))

# Calculate number of individuals in each group
print("Fig4A_PMC")
print(table(df$D28_MRD_category, useNA = "always"))

# Calculate p-values
print("Fig4A_PMC_MRD_ZBTB16")
wilcox.test(x = df$ZBTB16_logCPM[df$D28_MRD_category %in% c(">=5")], y = df$ZBTB16_logCPM[df$D28_MRD_category %in% c("0", "0-1", "1-5")], alternative = "greater")
print("Fig4A_PMC_MRD_Module")
wilcox.test(x = df$Module_score[df$D28_MRD_category %in% c(">=5")], y = df$Module_score[df$D28_MRD_category %in% c("0", "0-1", "1-5")], alternative = "greater")

# Boxplot of ZBTB16
p = ggplot(df, aes(x = D28_MRD_category, y = ZBTB16_logCPM)) +
  geom_boxplot(aes(fill = D28_MRD_category), col ="#000000", outlier.shape = NA, linewidth = 0.7) +
  scale_fill_manual(values = c("0" = "#1D71B8", "0-1" = "#20A7DB", "1-5" = "#A0D9EF", ">=5" = "#F77F11")) +
  geom_beeswarm(col ="#000000", alpha = 0.5, cex = 2.7, size = 2.5) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(),
    axis.ticks.length.x = unit(2, "mm"),
    axis.ticks.length.y = unit(2, "mm"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave("Plots/Fig4A_PMC_MRD_ZBTB16.pdf", p, width = 3.2, height = 3.5)

# Boxplot of module score
p = ggplot(df, aes(x = D28_MRD_category, y = Module_score)) +
  geom_boxplot(aes(fill = D28_MRD_category), col ="#000000", outlier.shape = NA, linewidth = 0.7) +
  scale_fill_manual(values = c("0" = "#1D71B8", "0-1" = "#20A7DB", "1-5" = "#A0D9EF", ">=5" = "#F77F11")) +
  geom_beeswarm(col ="#000000", alpha = 0.5, cex = 2.7, size = 2.5) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(),
    axis.ticks.length.x = unit(2, "mm"),
    axis.ticks.length.y = unit(2, "mm"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave("Plots/Fig4A_PMC_MRD_Module.pdf", p, width = 3.2, height = 3.5)





##### Fig4B: COG AALL0434 cohort - plotting against day 28 MRD categories

df = read.csv("Source_Data/Fig4BC_COG.csv")
df$D28_MRD_category = factor(df$D28_MRD_category, levels = c("0", "0-1", "1-5", ">=5"))

# Calculate number of individuals in each group
print("Fig4B_COG")
print(table(df$D28_MRD_category, useNA = "always"))

# Calculate p-values
print("Fig4B_COG_MRD_ZBTB16")
wilcox.test(x = df$ZBTB16_logCPM[df$D28_MRD_category %in% c(">=5")], y = df$ZBTB16_logCPM[df$D28_MRD_category %in% c("0", "0-1", "1-5")], alternative = "greater")
print("Fig4B_COG_MRD_Module")
wilcox.test(x = df$Module_score[df$D28_MRD_category %in% c(">=5")], y = df$Module_score[df$D28_MRD_category %in% c("0", "0-1", "1-5")], alternative = "greater")

# Boxplot of ZBTB16
p = ggplot(df %>% drop_na(D28_MRD_category), aes(x = D28_MRD_category, y = ZBTB16_logCPM)) +
  geom_boxplot(aes(fill = D28_MRD_category), col ="#000000", outlier.shape = NA, linewidth = 0.7) +
  scale_fill_manual(values = c("0" = "#1D71B8", "0-1" = "#20A7DB", "1-5" = "#A0D9EF", ">=5" = "#F77F11")) +
  geom_beeswarm(col ="#000000", alpha = 0.5, cex = 0.6, size = 1.8) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(),
    axis.ticks.length.x = unit(2, "mm"),
    axis.ticks.length.y = unit(2, "mm"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave("Plots/Fig4B_COG_MRD_ZBTB16.pdf", p, width = 3.2, height = 3.5)

# Boxplot of module score
p = ggplot(df %>% drop_na(D28_MRD_category), aes(x = D28_MRD_category, y = Module_score)) +
  geom_boxplot(aes(fill = D28_MRD_category), col ="#000000", outlier.shape = NA, linewidth = 0.7) +
  scale_fill_manual(values = c("0" = "#1D71B8", "0-1" = "#20A7DB", "1-5" = "#A0D9EF", ">=5" = "#F77F11")) +
  geom_beeswarm(col ="#000000", alpha = 0.5, cex = 0.6, size = 1.8) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(),
    axis.ticks.length.x = unit(2, "mm"),
    axis.ticks.length.y = unit(2, "mm"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave("Plots/Fig4B_COG_MRD_Module.pdf", p, width = 3.3, height = 3.5)





##### Fig4C: COG AALL0434 cohort - survival curves

df = read.csv("Source_Data/Fig4BC_COG.csv")

# Define OS_status as TRUE/FALSE/NA
df$OS_status = as.logical(df$OS_status)

# Define EFS_status as TRUE/FALSE/NA
df$EFS_status = as.logical(df$EFS_status)

# Stratify by ZBTB16_logCPM
df$ZBTB16_strata = NA
df$ZBTB16_strata[df$ZBTB16_logCPM > quantile(df$ZBTB16_logCPM, 0.6667, na.rm = TRUE)] = "ZBTB16_high"
df$ZBTB16_strata[df$ZBTB16_logCPM < quantile(df$ZBTB16_logCPM, 0.3333, na.rm = TRUE)] = "ZBTB16_low"

# Stratify by Module_score
df$Module_strata = NA
df$Module_strata[df$Module_score > quantile(df$Module_score, 0.6667, na.rm = TRUE)] = "Module_high"
df$Module_strata[df$Module_score < quantile(df$Module_score, 0.3333, na.rm = TRUE)] = "Module_low"

# Calculate number of individuals in each group
print("Fig4C_COG_Survival_strata")
print(table(df$ZBTB16_strata))
print(table(df$Module_strata))

# Calculate number of individuals in CoxPH models
print("Fig4D_Cox_ETP")
print(table(df$ZBTB16_strata, df$Flow_subtype))



### OS by ZBTB16_strata

# Fit survival
sfit = survfit(formula = Surv(OS, OS_status) ~ ZBTB16_strata, data = df)

# Extract data from from Survival object
strata = sfit$strata
strata = rep(gsub(pattern = ".*=", replacement = "", x = names(strata)), times = strata)
plotdf = data.frame(
  'Time' = sfit$time,
  'Prob' = sfit$surv,
  'Event' = sfit$n.censor,
  'Lower' = sfit$lower,
  'Upper' = sfit$upper,
  'Group' = strata
)

# Survival plot
p = ggplot(plotdf, aes(x = Time, y = Prob, colour = Group, fill = Group)) +
  geom_stepribbon(aes(ymin = Lower, ymax = Upper), colour = NA, alpha = 0.25) +
  geom_step() +
  geom_point(data = plotdf[plotdf$Event != 0, ], shape = 108, size = 2.5) +
  scale_x_continuous(limits = c(0, 4100), breaks = c(0, 2000, 4000), expand = expansion(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = expansion(0, 0)) +
  scale_colour_manual(values = c("#595959", "#309E88")) +
  scale_fill_manual(values = c("#595959", "#309E88")) +
  ggtitle(surv_pvalue(sfit)$pval) +
  theme_bw() +
  theme(
    axis.ticks.length.x = unit(2, "mm"), 
    axis.ticks.length.y = unit(2, "mm"), 
    legend.position = "bottom", 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
  )
ggsave("Plots/Fig4C_COG_Survival_OS_ZBTB16_strata.pdf", p, width = 4, height = 4.2)



### EFS by ZBTB16_strata

# Fit survival
sfit = survfit(formula = Surv(EFS, EFS_status) ~ ZBTB16_strata, data = df)

# Extract data from from Survival object
strata = sfit$strata
strata = rep(gsub(pattern = ".*=", replacement = "", x = names(strata)), times = strata)
plotdf = data.frame(
  'Time' = sfit$time,
  'Prob' = sfit$surv,
  'Event' = sfit$n.censor,
  'Lower' = sfit$lower,
  'Upper' = sfit$upper,
  'Group' = strata
)

# Survival plot
p = ggplot(plotdf, aes(x = Time, y = Prob, colour = Group, fill = Group)) +
  geom_stepribbon(aes(ymin = Lower, ymax = Upper), colour = NA, alpha = 0.25) +
  geom_step() +
  geom_point(data = plotdf[plotdf$Event != 0, ], shape = 108, size = 2.5) +
  scale_x_continuous(limits = c(0, 4100), breaks = c(0, 2000, 4000), expand = expansion(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = expansion(0, 0)) +
  scale_colour_manual(values = c("#595959", "#309E88")) +
  scale_fill_manual(values = c("#595959", "#309E88")) +
  ggtitle(surv_pvalue(sfit)$pval) +
  theme_bw() +
  theme(
    axis.ticks.length.x = unit(2, "mm"), 
    axis.ticks.length.y = unit(2, "mm"), 
    legend.position = "bottom", 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
  )
ggsave("Plots/Fig4C_COG_Survival_EFS_ZBTB16_strata.pdf", p, width = 4, height = 4.2)



### OS by Module_strata

# Fit survival
sfit = survfit(formula = Surv(OS, OS_status) ~ Module_strata, data = df)

# Extract data from from Survival object
strata = sfit$strata
strata = rep(gsub(pattern = ".*=", replacement = "", x = names(strata)), times = strata)
plotdf = data.frame(
  'Time' = sfit$time,
  'Prob' = sfit$surv,
  'Event' = sfit$n.censor,
  'Lower' = sfit$lower,
  'Upper' = sfit$upper,
  'Group' = strata
)

# Survival plot
p = ggplot(plotdf, aes(x = Time, y = Prob, colour = Group, fill = Group)) +
  geom_stepribbon(aes(ymin = Lower, ymax = Upper), colour = NA, alpha = 0.25) +
  geom_step() +
  geom_point(data = plotdf[plotdf$Event != 0, ], shape = 108, size = 2.5) +
  scale_x_continuous(limits = c(0, 4100), breaks = c(0, 2000, 4000), expand = expansion(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = expansion(0, 0)) +
  scale_colour_manual(values = c("#595959", "#309E88")) +
  scale_fill_manual(values = c("#595959", "#309E88")) +
  ggtitle(surv_pvalue(sfit)$pval) +
  theme_bw() +
  theme(
    axis.ticks.length.x = unit(2, "mm"), 
    axis.ticks.length.y = unit(2, "mm"), 
    legend.position = "bottom", 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
  )
ggsave("Plots/Fig4C_COG_Survival_OS_Module_strata.pdf", p, width = 4, height = 4.2)



### EFS by Module_strata

# Fit survival
sfit = survfit(formula = Surv(EFS, EFS_status) ~ Module_strata, data = df)

# Extract data from from Survival object
strata = sfit$strata
strata = rep(gsub(pattern = ".*=", replacement = "", x = names(strata)), times = strata)
plotdf = data.frame(
  'Time' = sfit$time,
  'Prob' = sfit$surv,
  'Event' = sfit$n.censor,
  'Lower' = sfit$lower,
  'Upper' = sfit$upper,
  'Group' = strata
)

# Survival plot
p = ggplot(plotdf, aes(x = Time, y = Prob, colour = Group, fill = Group)) +
  geom_stepribbon(aes(ymin = Lower, ymax = Upper), colour = NA, alpha = 0.25) +
  geom_step() +
  geom_point(data = plotdf[plotdf$Event != 0, ], shape = 108, size = 2.5) +
  scale_x_continuous(limits = c(0, 4100), breaks = c(0, 2000, 4000), expand = expansion(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = expansion(0, 0)) +
  scale_colour_manual(values = c("#595959", "#309E88")) +
  scale_fill_manual(values = c("#595959", "#309E88")) +
  ggtitle(surv_pvalue(sfit)$pval) +
  theme_bw() +
  theme(
    axis.ticks.length.x = unit(2, "mm"), 
    axis.ticks.length.y = unit(2, "mm"), 
    legend.position = "bottom", 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
  )
ggsave("Plots/Fig4C_COG_Survival_EFS_Module_strata.pdf", p, width = 4, height = 4.2)





##### Fig4D: COG AALL0434 cohort - Cox-PH analysis of ZBTB16 (tertiles) vs ETP

df = read.csv("Source_Data/Fig4D_Cox_ETP.csv")
df = df[df$model %in% c("ZBTB16_tertiles_ETP_exclude", "ZBTB16_tertiles_ETP_include"), ]
df$model = factor(df$model, levels = c("ZBTB16_tertiles_ETP_exclude", "ZBTB16_tertiles_ETP_include"))
df$variable = factor(df$variable, levels = c("ETP", "ZBTB16"))
df$pval_label = round(-log10(df$pval), digits = 1)
df$pval_sig = "No"
df$pval_sig[df$pval < 0.05] = "Yes"

p1 = ggplot(df[df$outcome == "OS", ], aes(y = variable)) + 
  geom_text(aes(label = pval_label, col = pval_sig), x = 0.2, hjust = 0, vjust = 0.5) +
  facet_grid(rows = vars(model)) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("No" = "black", "Yes" = "red")) +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    
    strip.background = element_blank(), 
    strip.text.y = element_blank(),
    legend.position = "none"
  )

p2 = ggplot(df[df$outcome == "OS", ], aes(y = variable)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(aes(x = hr), shape = 15) +
  geom_segment(aes(x = lower_ci_95, xend = upper_ci_95)) +
  facet_grid(rows = vars(model)) +
  scale_x_continuous(breaks = 0:5, labels = 0:5, limits = c(0, 5), expand = expansion(c(0, 0))) +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    
    axis.title.x = element_blank(),
    axis.ticks.length.x = unit(2, "mm"), 
    
    strip.background = element_blank(), 
    strip.text.y = element_text(size = 6, angle = 0, hjust = 0, vjust = 0.5)
  )

p3 = ggplot(df[df$outcome == "EFS", ], aes(y = variable)) + 
  geom_text(aes(label = pval_label, col = pval_sig), x = 0.2, hjust = 0, vjust = 0.5) +
  facet_grid(rows = vars(model)) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("No" = "black", "Yes" = "red")) +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    
    strip.background = element_blank(), 
    strip.text.y = element_blank(),
    legend.position = "none"
  )

p4 = ggplot(df[df$outcome == "EFS", ], aes(y = variable)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(aes(x = hr), shape = 15) +
  geom_segment(aes(x = lower_ci_95, xend = upper_ci_95)) +
  facet_grid(rows = vars(model)) +
  scale_x_continuous(breaks = 0:4, labels = 0:4, limits = c(0, 4), expand = expansion(c(0, 0))) +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    
    axis.title.x = element_blank(),
    axis.ticks.length.x = unit(2, "mm"), 
    
    strip.background = element_blank(), 
    strip.text.y = element_text(size = 6, angle = 0, hjust = 0, vjust = 0.5)
  )

p = p1 + p2 + p3 + p4 + plot_layout(nrow = 1, ncol = 4, widths = c(1, 5, 1, 5), guides = "keep")
ggsave("Plots/Fig4D_Cox_ETP_tertiles.pdf", p, width = 8, height = 1.5)





##### Fig4E: COG AALL0434 cohort - Cox-PH analysis of ZBTB16 vs BMP

df = read.csv("Source_Data/Fig4E_Cox_BMP.csv")
df$model = factor(df$model, levels = c("ZBTB16_module_BMP_119", "ZBTB16_module_BMP_17", "ZBTB16_module_BMP_9"))
df$variable = factor(df$variable, levels = c("BMP", "ZBTB16"))
df$pval_label = round(-log10(df$pval), digits = 1)
df$pval_sig = "No"
df$pval_sig[df$pval < 0.05] = "Yes"

p1 = ggplot(df[df$outcome == "OS", ], aes(y = variable)) + 
  geom_text(aes(label = pval_label, col = pval_sig), x = 0.2, hjust = 0, vjust = 0.5) +
  facet_grid(rows = vars(model)) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("No" = "black", "Yes" = "red")) +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    
    strip.background = element_blank(), 
    strip.text.y = element_blank(),
    legend.position = "none"
  )

p2 = ggplot(df[df$outcome == "OS", ], aes(y = variable)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(aes(x = hr), shape = 15) +
  geom_segment(aes(x = lower_ci_95, xend = upper_ci_95)) +
  facet_grid(rows = vars(model)) +
  scale_x_continuous(breaks = 0:2, labels = 0:2, limits = c(0, 2), expand = expansion(c(0, 0))) +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    
    axis.title.x = element_blank(),
    axis.ticks.length.x = unit(2, "mm"), 
    
    strip.background = element_blank(), 
    strip.text.y = element_text(size = 6, angle = 0, hjust = 0, vjust = 0.5)
  )

p3 = ggplot(df[df$outcome == "EFS", ], aes(y = variable)) + 
  geom_text(aes(label = pval_label, col = pval_sig), x = 0.2, hjust = 0, vjust = 0.5) +
  facet_grid(rows = vars(model)) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("No" = "black", "Yes" = "red")) +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    
    strip.background = element_blank(), 
    strip.text.y = element_blank(),
    legend.position = "none"
  )

p4 = ggplot(df[df$outcome == "EFS", ], aes(y = variable)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(aes(x = hr), shape = 15) +
  geom_segment(aes(x = lower_ci_95, xend = upper_ci_95)) +
  facet_grid(rows = vars(model)) +
  scale_x_continuous(breaks = 0:2, labels = 0:2, limits = c(0, 2), expand = expansion(c(0, 0))) +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    
    axis.title.x = element_blank(),
    axis.ticks.length.x = unit(2, "mm"), 
    
    strip.background = element_blank(), 
    strip.text.y = element_text(size = 6, angle = 0, hjust = 0, vjust = 0.5)
  )

p = p1 + p2 + p3 + p4 + plot_layout(nrow = 1, ncol = 4, widths = c(1.5, 4.5, 1.5, 4.5), guides = "keep")
ggsave("Plots/Fig4E_Cox_BMP.pdf", p, width = 7.5, height = 2)


