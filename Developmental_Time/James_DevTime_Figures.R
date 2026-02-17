# Plotting developmental times and pupation times
# Script by Saudat Alishayeva, tweaked by James Tan

# Last Updated: 23/9/25

################################
##### Packages and Setup #######
################################

rm(list = ls())

library(ggplot2)
library(gplots)
library(ggpubr)
library(dplyr)
library(svglite)


#pupation dataset from Dila:
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\Developmental_Time")
pupae_compl_obs_long<-readRDS(file = 'Average_pupation_day_NEX_3006.RDS')
pupae_compl_obs_long$Cond <- ifelse(pupae_compl_obs_long$New_Combined_Cond_Name == 'Control Diet','Ctrl','HS')
palette = c("Ctrl" = "#C6B49E", "HS" = "#DF9F65")
logOddstest <- paste0('GLM log-odds = -1.62','\n','p < 0.001')
Pupation_Age_Plot <- ggplot(pupae_compl_obs_long, aes(x = pupae_compl_obs, y = count, color = Cond)) +
  geom_smooth(method = "loess", se = TRUE, size = 0.6, aes(fill = Cond), alpha = 0.2) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette, guide = "none") +
  labs(
    x = "Days after egg lay",
    y = "Pupae count",
    color = "Diet"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 0)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.position = c(0.85, 0.85),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 7),
    legend.background = element_blank()
  ) + 
  annotate("text", y = 0, x = 4.75, label = logOddstest, size = 3.5, hjust = 0,vjust = 0.5)
Pupation_Age_Plot
ggsave(filename = 'PupationAgePlot.svg',Pupation_Age_Plot,width = 2.5, height = 3,dpi=300)

# Save the cleaned eclosion dataframe 
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\Developmental_Time")
pupation_age_clean <- pupae_compl_obs_long[,c('REP_Group','pupae_compl_obs','Cond','count')]
colnames(pupation_age_clean) <- c('Replicate_ID','Days_after_egg_lay','Diet','Count')
write.csv(file = 'Pupation_age_clean.csv',pupation_age_clean)

# Survival dataset - 2025 experiment
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\Developmental_Time")
survival <- readRDS('220925_survival.RDS')
survival$condition = stringr::str_split_i(survival$Replica_ID, "_", 2)
survival$condition = gsub('CTRL', 'Ctrl', survival$condition)
survival$Survival_eclosion = gsub(',', '.', survival$Survival)
survival$Survival_eclosion = as.numeric(survival$Survival_eclosion)
survival_sbs = na.omit(survival) # One replicate is faulty
survival_sbs$Eggs = as.numeric(survival_sbs$Eggs)

# Plot
# Calculate statistics
wilcox_test <- wilcox.test(Survival_eclosion ~ condition, data = survival_sbs)
p_value <- format.pval(wilcox_test$p.value, digits = 2)
medians <- survival_sbs %>% 
  group_by(condition) %>% 
  summarise(median_val = median(Survival_eclosion))
median_diff <- round(medians$median_val[2] - medians$median_val[1], 2)  # Rounded to 2 digits
median_diff # 13%
palette = c("Ctrl" = "#C6B49E", "HS" = "#DF9F65")

EggAdultSurvivalPlot <- ggplot(data = survival_sbs, aes(x = condition, y = Survival_eclosion, fill = condition)) +
  geom_boxplot(alpha = 0.8,color = "darkgrey") +  # Remove legend
  geom_jitter(width = 0.2, size = 2, color = "darkgrey") +
  scale_fill_manual(values = palette) +
  labs(
    x = "Diet",
    y = "Egg-to-adult survival %"
  ) +
  theme_classic(base_size = 14) +
  ylim(50, max(survival_sbs$Survival_eclosion) * 1.1)+
  # Significance bar (p-value)
  geom_signif(
    comparisons = list(c("Ctrl", "HS")),
    annotations = paste0("p = ", p_value),
    y_position = max(survival_sbs$Survival_eclosion) * 1.02,  # Bar position
    tip_length = 0.01,
    textsize = 4,
    vjust = -0.5
  ) +
  theme(legend.position = "none")

EggAdultSurvivalPlot

# Save plot
ggsave(filename = 'EggAdultSurvivalPlot.svg',EggAdultSurvivalPlot,width = 2, height = 3,dpi=300)


# Effect of number of eggs on survival in each condition
survival_sbs_Ctrl <- subset(survival_sbs,condition=='Ctrl')
Ctrl_corr <- cor.test(survival_sbs_Ctrl$Egg,survival_sbs_Ctrl$Survival_eclosion)
cor = round(Ctrl_corr$estimate,3)
pval = round(Ctrl_corr$p.value,3)

Ctrl_corr_plot <- ggplot(survival_sbs_Ctrl,aes(y=Survival_eclosion,x=Eggs))+
  geom_point(color='darkgrey',alpha=0.8)+
  ylab('Egg-to-adult survival %')+
  xlab('Number of eggs per replicate (Ctrl)')+
  annotate("text", x = 90, y = 90, label = paste0("Pearson r = ",cor," p-value = ",pval), size = 3.5, vjust=1)+
  geom_smooth(method = "lm", col = "blue",linetype='dashed')+
  theme_classic()
Ctrl_corr_plot

# Repeat for HS
survival_sbs_HS <- subset(survival_sbs,condition=='HS')
HS_corr <- cor.test(survival_sbs_HS$Egg,survival_sbs_HS$Survival_eclosion)
cor = round(HS_corr$estimate,3)
pval = round(HS_corr$p.value,3)

HS_corr_plot <- ggplot(survival_sbs_HS,aes(y=Survival_eclosion,x=Eggs))+
  geom_point(color='darkgrey',alpha=0.8)+
  ylab('Egg-to-adult survival %')+
  xlab('Number of eggs per replicate (HS)')+
  annotate("text", x = 110, y = 100, label = paste0("Pearson r = ",cor," p-value = ",pval), size = 3.5, vjust=1)+
  geom_smooth(method = "lm", col = "blue",linetype='dashed')+
  theme_classic()
HS_corr_plot

# Save the two plots
ggsave(filename = 'Ctrl_corr_plot.svg',Ctrl_corr_plot,width = 3, height = 3,dpi=300)
ggsave(filename = 'HS_corr_plot.svg',HS_corr_plot,width = 3, height = 3,dpi=300)


# Cleaned dataframe already available

