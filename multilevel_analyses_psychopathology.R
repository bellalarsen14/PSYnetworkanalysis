############################################################
############################################################
############### Network Analysis, XXY and XYY ##############
############################################################
############################################################

### SET UP

library(magrittr)
library(dplyr)
library(ggplot2)
library(gplots)
library(forcats)
library(superheat)
library(RColorBrewer)
library(stringr)
library(reshape2)
library(psych)
library(igraph)
library(scales)
library(moments)
library(modelr)
library(purrr)
library(broom)
library(plotrix)
library(tidyr)
library(factoextra)
library(viridis)
library(ggpubr)
library(ggrepel)
library(cluster)
library(blockmodels)
library(Hmisc)
library(tibble)
library(purrrlyr)
library(skimr)
library(cocor)
library(devtools)
library(superheat)
library(networkD3)
library(ggsankey)
library(scmap)
library(deming)
library(SimplyAgree)
library(forcats)
library(poolr)

### Load in correlation matrices - two 53x53 matrices (XXY/KS and XYY) - Pearson Correlations between scales
matrix_xxy <- read.csv("matrix_xxy.csv")
matrix_xxy <- read.csv("matrix_xyy.csv")

### Load in effect size dataframes
cocor_xxy_df_eff_nodal <- read.csv("cocor_xxy_df_eff_nodal.csv")
cocor_xyy_df_eff_nodal <- read.csv("cocor_xyy_df_eff_nodal.csv")

### Load in permuted p-values
divergence_p_perm_long <- read.csv("divergence_p_perm_long.csv")
ag_p_diff_p_perm_long <- read.csv("ag_p_diff_p_perm_long.csv")

### Transform correlation coefficients using Fisher's Z
fisher_z_transform <- function(r) {
  fisherz(r)
} 

diag(matrix_xxy) <- NA
diag(matrix_xyy) <- NA

fisher_matrix_xxy <- matrix(sapply(matrix_xxy, fisher_z_transform), nrow = ncol(matrix_xxy))

fisher_matrix_xyy <- matrix(sapply(matrix_xyy, fisher_z_transform), nrow = ncol(matrix_xyy))

# change the row and column names to scale names
rownames(fisher_matrix_xxy) <-colnames(matrix_xxy)
colnames(fisher_matrix_xxy) <- colnames(matrix_xxy)

rownames(fisher_matrix_xyy) <-colnames(matrix_xyy)
colnames(fisher_matrix_xyy) <- colnames(matrix_xyy)
#===============================================================================
### Comparing the strength of each scale’s connectivity with all others

## Scale-level global coupling scores 

# 1: calculate nodal degree (row mean) for all scales
nodal_xxy <- rowMeans(fisher_matrix_xxy, na.rm = T) # the diagonal is NA
nodal_xyy <- rowMeans(fisher_matrix_xyy, na.rm = T)

# create dataframe of both disorder row means
df_nodal_xxy <- as.data.frame(nodal_xxy) %>% 
  mutate(group = "XXY")
df_nodal_xyy <- as.data.frame(nodal_xyy) %>% 
  mutate(group = "XYY")

df_nodal_xxy <- tibble::rownames_to_column(df_nodal_xxy, "scale_name") %>% 
  select(-group)
df_nodal_xyy <- tibble::rownames_to_column(df_nodal_xyy, "scale_name") %>% 
  select(-scale_name, -group)

df_nodal_total <- bind_cols(df_nodal_xxy, df_nodal_xyy) 

# 2: calculate the correlation between scale-level global coupling scores in XXY/KS and XYY
cor.test(df_nodal_total$nodal_xxy, df_nodal_total$nodal_xyy, method="pearson") # r = 0.78

# 3: Deming regression of nodal degree in XXY vs nodal degree in XYY
deming_packageSimplyAgree <- dem_reg(x="nodal_xxy", y="nodal_xyy", data = df_nodal_total) # slope: p = .661; intercept: p = .875
deming_slope <- deming_packageSimplyAgree$model$coef[2]
deming_intercept <- deming_packageSimplyAgree$model$coef[1]

# 3A: Extract fitted values
fitted_values <- deming_intercept + deming_slope*nodal_xxy

# 3B: Calculate residuals
residuals_nd <- nodal_xyy - fitted_values

# 3C: Display the residuals as a dataframe
residuals_nd <- as.data.frame(residuals_nd)
View(residuals_nd)

## Create Fig. 2A
# 1: color coding 
# 1A: define the scales
scales <- df_nodal_total %>% 
  mutate(scale_title = case_when(str_detect(scale_name, "CON") ~ "CON: Conners-3",
                                 str_detect(scale_name, "CBCL") ~ "CBCL: Child Behavior Checklist",
                                 str_detect(scale_name, "SRS") ~ "SRS: Social Responsiveness Scale",
                                 str_detect(scale_name, "SDQ") ~ "SDQ: Strengths and Difficulties Questionnaire",
                                 str_detect(scale_name, "SHRP") ~ "C-SHARP: Children’s Scale of Hostility and Aggression: Reactive/Proactive",
                                 str_detect(scale_name, "DCD") ~ "DCDQ: Developmental Coordination Disorder Questionnaire",
                                 str_detect(scale_name, "APSD") ~ "APSD: Antisocial Process Screening Device",
                                 str_detect(scale_name, "ARI") ~ "ARI: Affective Reactivity Index",
                                 str_detect(scale_name, "SCQ") ~ "SCQ: Social Communication Questionnaire"
  )) %>% 
  select(scale_name, scale_title)

scale_names <- scales %>% 
  pull(scale_name)
scale_titles <- scales %>% 
  distinct(scale_title) %>% 
  pull(scale_title)

df_nodal_total_color <- df_nodal_total %>% 
  mutate(scale_title = case_when(str_detect(scale_name, "CON") ~ "CON: Conners-3",
                                 str_detect(scale_name, "CBCL") ~ "CBCL: Child Behavior Checklist",
                                 str_detect(scale_name, "SRS") ~ "SRS: Social Responsiveness Scale",
                                 str_detect(scale_name, "SDQ") ~ "SDQ: Strengths and Difficulties Questionnaire",
                                 str_detect(scale_name, "SHRP") ~ "C-SHARP: Children’s Scale of Hostility and Aggression: Reactive/Proactive",
                                 str_detect(scale_name, "DCD") ~ "DCDQ: Developmental Coordination Disorder Questionnaire",
                                 str_detect(scale_name, "APSD") ~ "APSD: Antisocial Process Screening Device",
                                 str_detect(scale_name, "ARI") ~ "ARI: Affective Reactivity Index",
                                 str_detect(scale_name, "SCQ") ~ "SCQ: Social Communication Questionnaire"
  )) 

instrument_colors <- structure(c("#1f78b4", "#33a02c", "#984ea3", "#e31a1c",
                                 "#ff7f00", "#fb9a99", "#35978f", "#a65628",
                                 "black"), 
                               names = scale_titles)
# 2: plot Fig. 2A

figure_2a <- df_nodal_total_color %>% 
  ggplot(., aes(y=nodal_xyy, x=nodal_xxy, col=scale_title, label = scale_name))+
  geom_point(alpha=1,size=3)+
  geom_abline(slope=1.047404, intercept=0.008037, color="blue", linewidth=1)+
  geom_abline(intercept=0, slope = 1, lty=2, alpha=0.5, linewidth=1)+
  scale_color_manual(values = instrument_colors)+
  geom_text_repel(max.overlaps = 5)+
  theme(axis.text.x = element_text(size=14))+
  theme(axis.text.y = element_text(size=14))+
  theme(axis.title.x = element_text(size=14))+
  theme(axis.title.y = element_text(size=14))+
  theme(plot.subtitle = element_text(size=14))+
  theme(plot.title = element_text(size=16))+
  theme(legend.position="none")+
  labs(x = "Nodal degree: XXY", y = "Nodal degree: XYY", title = "Deming Regression of Nodal Degree",
       subtitle = "p = .661, r = 0.78")
figure_2a
#===============================================================================
## Comparing the profile of each scale’s connectivity with all others

# 1: create a 53x1 matrix of the correlation between edges in XXY/KS and XYY
f <- function(r1,r2){
  cor(r1,r2, use = "na.or.complete")
}

edge_corr <- map2_dfr(as.data.frame(fisher_matrix_xxy), as.data.frame(fisher_matrix_xyy), f)

edge_corr <- edge_corr %>% 
  pivot_longer(cols = 1:53, names_to = "scale_name", values_to = "correlation")

edge_corr_divergence <- edge_corr %>% 
  mutate(divergence = 1-correlation)

edge_corr_divergence %>% 
  summarise(min(divergence), max(divergence))

# 2: merge in permuted p-values
divergence_p_perm_long <- divergence_p_perm_long %>% 
  select(-observed)
edge_corr_all <- full_join(edge_corr_divergence, divergence_p_perm_long, by = "scale_name") %>% 
  mutate(p_perm_div = p_perm,
         p_perm.bonf_div = p_perm.bonf) %>% 
  select(-p_perm,-p_perm.bonf)

# 2A: create new variable to denote statistical significance
edge_corr_all <- edge_corr_all %>% 
  mutate(significance_div = case_when(p_perm_div < 0.05 ~ "Nominally Significant (p < .05)",
                                      p_perm.bonf_div < 0.05 ~ "Significant After Bonferroni Correction (p < 0.0009)",
                                      TRUE ~ "Nonsignificant"))

edge_corr_all <- edge_corr_all %>% 
  mutate(scale_title = case_when(str_detect(scale_name, "CON") ~ "CON: Conners-3",
                                 str_detect(scale_name, "CBCL") ~ "CBCL: Child Behavior Checklist",
                                 str_detect(scale_name, "SRS") ~ "SRS: Social Responsiveness Scale",
                                 str_detect(scale_name, "SDQ") ~ "SDQ: Strengths and Difficulties Questionnaire",
                                 str_detect(scale_name, "SHRP") ~ "C-SHARP: Children’s Scale of Hostility and Aggression: Reactive/Proactive",
                                 str_detect(scale_name, "DCD") ~ "DCDQ: Developmental Coordination Disorder Questionnaire",
                                 str_detect(scale_name, "APSD") ~ "APSD: Antisocial Process Screening Device",
                                 str_detect(scale_name, "ARI") ~ "ARI: Affective Reactivity Index",
                                 str_detect(scale_name, "SCQ") ~ "SCQ: Social Communication Questionnaire"
  )) 

## Create Fig. 2B

# 1: create legend
legend_edge <- ggplot(edge_corr_all, aes(x=reorder(scale_name,-desc(correlation)), y = 0.1)) + 
  geom_point(aes(color = scale_title), shape = 15, size = 3, show.legend = F) + 
  theme_classic()+
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm")) +
  labs(color = "Scale Name") +
  scale_color_manual(values = instrument_colors)

# 2: create plot
edgewise_div <- ggplot(edge_corr_all, aes(x = reorder(scale_name,-desc(divergence)), y = divergence, color = significance_div))+
  geom_point(size=2.5)+ 
  scale_color_manual(values = c("Nominally Significant (p < .05)" = "turquoise4",
                                "Nonsignificant" = "#999999",
                                "Significant After Bonferroni Correction (p < 0.0009)" = "#d8b365")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14,
                                   face = c("plain", "plain", "plain", "plain",
                                            "plain", "plain", "plain", "plain",
                                            "plain", "plain", "plain", "plain",
                                            "plain", "plain", "plain", "plain",
                                            "plain", "plain", "plain", "plain",
                                            "plain", "plain", "plain", "plain",
                                            "plain", "plain", "plain", "plain",
                                            "plain", "plain", "plain", "plain",
                                            "plain", "plain", "plain", "bold",
                                            "plain", "plain", "plain", "bold",
                                            "plain", "plain", "plain", "plain",
                                            "plain", "plain", "plain", "plain",
                                            "plain", "bold", "plain", "bold",
                                            "plain")),
        axis.title=element_text(size=14),
        axis.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14),
        legend.position="right",
        plot.title = element_text(size=16))+
  labs(title = "Fingerprint divergence score",
       x = "Scale Name", y = "Divergence (1 - Pearson's r)",
       color = "Statistical significance")

# 3: combine legend and plot
edgewise_div_leg <- edgewise_div + annotation_custom(ggplotGrob(legend_edge),
                                                     ymin = 0.035, ymax = 0.055, 
                                                     xmin = 0.2, xmax = 53.6) 
edgewise_div_leg

## Create Fig. 2C

# 1: pull out ag.p_SHRP, using raw correlations, not fisher's Z adjusted
diag(matrix_xxy) <- NA
diag(matrix_xyy) <- NA

both_fish_ag_p <- as.data.frame(matrix_xxy) %>% 
  select(ag.p_SHRP)
both_fish_ag_p_y <- as.data.frame(matrix_xyy) %>% 
  select(ag.p_SHRP)
both_fish_ag_p <- bind_cols(both_fish_ag_p,both_fish_ag_p_y)
both_fish_ag_p <- rownames_to_column(both_fish_ag_p)
colnames(both_fish_ag_p) <- c("scale_name","scale_xxy","scale_xyy")
both_fish_ag_p <- both_fish_ag_p %>% 
  filter(scale_name != "ag.p_SHRP") %>% 
  mutate(divergence_xxy = 1-scale_xxy,
         divergence_xyy = 1-scale_xyy) %>% 
  mutate(scale_title = case_when(str_detect(scale_name, "CON") ~ "CON: Conners-3",
                                 str_detect(scale_name, "CBCL") ~ "CBCL: Child Behavior Checklist",
                                 str_detect(scale_name, "SRS") ~ "SRS: Social Responsiveness Scale",
                                 str_detect(scale_name, "SDQ") ~ "SDQ: Strengths and Difficulties Questionnaire",
                                 str_detect(scale_name, "SHRP") ~ "C-SHARP: Children’s Scale of Hostility and Aggression: Reactive/Proactive",
                                 str_detect(scale_name, "DCD") ~ "DCDQ: Developmental Coordination Disorder Questionnaire",
                                 str_detect(scale_name, "APSD") ~ "APSD: Antisocial Process Screening Device",
                                 str_detect(scale_name, "ARI") ~ "ARI: Affective Reactivity Index",
                                 str_detect(scale_name, "SCQ") ~ "SCQ: Social Communication Questionnaire"
  ))

#deming 
deming_packageSimplyAgree_ag_p <- dem_reg(x="scale_xxy", y="scale_xyy", data = both_fish_ag_p) 

# use the permuted p-values to determine which scales are significantly more coupled to ag.p_SHRP in one group than the other 
ag_p_diff_p_perm_long %>% 
  filter(p_perm < 0.05)

# generate left panel of Fig. 2C
scatterplot_ag_p_SHRP <- both_fish_ag_p %>% 
  ggplot(., aes(y=scale_xyy, x=scale_xxy))+
  geom_point(aes(color = scale_title), alpha=1, size=3)+
  geom_text(data = subset(both_fish_ag_p, scale_name == "tot_APSD" | 
                            scale_name == "narc_APSD" | scale_name == "cond_SDQ" 
                          | scale_name == "soc_CBCL"| scale_name == "ax.dep_CBCL"
                          | scale_name == "ax.dep_CBCL"| scale_name == "int_CBCL"
                          | scale_name == "tot_CBCL"| scale_name == "som_CBCL"
                          | scale_name == "def_CON" | scale_name == "bul_SHRP"
                          | scale_name == "hos_SHRP" | scale_name == "rule_CBCL"
                          | scale_name == "tot_ARI"| scale_name == "adhd.att_CON"
                          | scale_name == "com_SCQ"| scale_name == "exec_CON"
                          | scale_name == "peer_CON" | scale_name == "peer_SDQ"
                          | scale_name == "ext_CBCL" | scale_name == "att_CBCL"
                          | scale_name == "att_CON" | scale_name == "wt.dep_CBCL"), 
            aes(label = scale_name), vjust = -1)+ 
  theme(axis.text.x = element_text(size=15))+
  theme(axis.text.y = element_text(size=15))+
  theme(legend.position = "none") +
  geom_abline(slope=2.2591, intercept=-0.3195, color="blue", linewidth=1)+
  geom_abline(intercept=0, slope = 1, lty=2, alpha=0.5, linewidth=1)+
  scale_color_manual(values = instrument_colors)+
  labs(x = "XXY/KS: Correlation with ag.p_SHRP",
       y = "XYY: Correlation with ag.p_SHRP",
       title = "ag.p_SHRP",
       subtitle = "Deming regression slope = 2.26; p = .300")
scatterplot_ag_p_SHRP

# 2: pull out soc_CBCL, using raw correlations, not fisher's Z adjusted

both_fish_soc_CBCL <- as.data.frame(matrix_xxy) %>% 
  select(soc_CBCL)
both_fish_soc_CBCL_y <- as.data.frame(matrix_xyy) %>% 
  select(soc_CBCL)
both_fish_soc_CBCL <- bind_cols(both_fish_soc_CBCL,both_fish_soc_CBCL_y)
both_fish_soc_CBCL <- rownames_to_column(both_fish_soc_CBCL)
colnames(both_fish_soc_CBCL) <- c("scale_name","scale_xxy","scale_xyy")
both_fish_soc_CBCL <- both_fish_soc_CBCL %>% 
  filter(scale_name != "soc_CBCL") %>% 
  mutate(divergence_xxy = 1-scale_xxy,
         divergence_xyy = 1-scale_xyy) %>% 
  mutate(scale_title = case_when(str_detect(scale_name, "CON") ~ "CON: Conners-3",
                                 str_detect(scale_name, "CBCL") ~ "CBCL: Child Behavior Checklist",
                                 str_detect(scale_name, "SRS") ~ "SRS: Social Responsiveness Scale",
                                 str_detect(scale_name, "SDQ") ~ "SDQ: Strengths and Difficulties Questionnaire",
                                 str_detect(scale_name, "SHRP") ~ "C-SHARP: Children’s Scale of Hostility and Aggression: Reactive/Proactive",
                                 str_detect(scale_name, "DCD") ~ "DCDQ: Developmental Coordination Disorder Questionnaire",
                                 str_detect(scale_name, "APSD") ~ "APSD: Antisocial Process Screening Device",
                                 str_detect(scale_name, "ARI") ~ "ARI: Affective Reactivity Index",
                                 str_detect(scale_name, "SCQ") ~ "SCQ: Social Communication Questionnaire"
  ))


#deming 
deming_packageSimplyAgree_soc <- dem_reg(x="scale_xxy", y="scale_xyy", data = both_fish_soc_CBCL) 

# use the permuted p-values to determine which scales are significantly more coupled to soc_CBCL in one group than the other 
soc_CBCL_p_perm_long %>% 
  filter(p_perm < 0.05)

# generate right panel of Fig. 2C
scatterplot_soc_CBCL <- both_fish_soc_CBCL %>% 
  ggplot(., aes(y=scale_xyy, x=scale_xxy))+
  geom_point(aes(color = scale_title), alpha=1, size=3)+
  geom_text(data = subset(both_fish_soc_CBCL, scale_name == "ag.p_SHRP" | 
                            scale_name == "imp_APSD" | scale_name == "att_CON" 
                          | scale_name == "adhd.att_CON"| scale_name == "cond_CON"
                          | scale_name == "rep_SCQ"| scale_name == "hyp.imp_SDQ"
                          | scale_name == "ax.dep_CBCL"| scale_name == "tot_SCQ"
                          | scale_name == "wt.dep_CBCL"
                          | scale_name == "def_CON" | scale_name == "ax.dep_CBCL"
                          | scale_name == "bul_SHRP"| scale_name == "peer_SDQ"
                          | scale_name == "mot_SRS"| scale_name == "tot_CBCL"
                          | scale_name == "tot_SDQ" | scale_name == "learn_CON"
                          | scale_name == "narc_APSD" | scale_name == "exec_CON"), 
            aes(label = scale_name), vjust = -1)+ 
  theme(axis.text.x = element_text(size=15))+
  theme(axis.text.y = element_text(size=15))+
  theme(legend.position = "none") +
  geom_abline(slope=1.00330, intercept=-0.01058, color="blue", linewidth=1)+
  geom_abline(intercept=0, slope = 1, lty=2, alpha=0.5, linewidth=1)+
  scale_color_manual(values = instrument_colors)+
  labs(x = "XXY/KS: Correlation with soc_CBCL",
       y = "XYY: Correlation with soc_CBCL",
       title = "Correlation of all scales with soc_CBCL",
       subtitle = "Deming regression slope = 1.00; p = .993")
scatterplot_soc_CBCL

## Create Fig. 2D

# 1: relate edge-wise correlation with differences in effect size (compared to XY) between groups
delta_eff_xxy <- cocor_xxy_df_eff_nodal %>% 
  mutate(effect_size_xxy = effect_size) %>% 
  select(scale_name, effect_size_xxy)
delta_eff_xyy <- cocor_xyy_df_eff_nodal %>% 
  mutate(effect_size_xyy = effect_size) %>% 
  select(scale_name, effect_size_xyy)

delta_eff <- full_join(delta_eff_xxy,delta_eff_xyy, by = "scale_name")
delta_eff <- delta_eff %>% 
  mutate(delta_eff = effect_size_xxy-effect_size_xyy) %>% 
  select(scale_name, delta_eff)

edge_corr_eff <- full_join(edge_corr_divergence, delta_eff, by = c("scale_name"))

edge_corr_eff <- edge_corr_eff %>% 
  mutate(abs_delta_eff = abs(delta_eff)) %>% 
  mutate(scale_title = case_when(str_detect(scale_name, "CON") ~ "CON: Conners-3",
                                 str_detect(scale_name, "CBCL") ~ "CBCL: Child Behavior Checklist",
                                 str_detect(scale_name, "SRS") ~ "SRS: Social Responsiveness Scale",
                                 str_detect(scale_name, "SDQ") ~ "SDQ: Strengths and Difficulties Questionnaire",
                                 str_detect(scale_name, "SHRP") ~ "C-SHARP: Children’s Scale of Hostility and Aggression: Reactive/Proactive",
                                 str_detect(scale_name, "DCD") ~ "DCDQ: Developmental Coordination Disorder Questionnaire",
                                 str_detect(scale_name, "APSD") ~ "APSD: Antisocial Process Screening Device",
                                 str_detect(scale_name, "ARI") ~ "ARI: Affective Reactivity Index",
                                 str_detect(scale_name, "SCQ") ~ "SCQ: Social Communication Questionnaire"
  )) 

# 2: correlate the absolute difference in effect size per group 
edge_eff_corr <- map2_dfr(as.data.frame(edge_corr_eff$abs_delta_eff), as.data.frame(edge_corr_eff$divergence), f)
edge_eff_corr <- pull(edge_eff_corr)

# 3: plot

scatterplot_eff_edge_corr <- edge_corr_eff %>% 
  ggplot(., aes(y=divergence, x=abs_delta_eff, label = scale_name))+
  geom_point(aes(color = scale_title), alpha=1)+
  theme(axis.text.x = element_text(size=15))+
  theme(axis.text.y = element_text(size=15))+
  theme(legend.position = "none")+
  geom_text_repel(max.overlaps = 4)+
  geom_smooth(method = "lm",se = FALSE) + 
  scale_color_manual(values = instrument_colors) +
  xlab("Absolute difference in effect size |XXY-XYY|")+
  ylab("Divergence (1 - Correlation between XXY and XYY)") +
  labs(title = "Divergence is not related to absolute effect size",
       subtitle = "r = -0.01")
scatterplot_eff_edge_corr

#===============================================================================
## Comparing the full distribution of pairwise coupling between all scales 

# 1: reduce matrices to half
upper_fisher_XXY <- as.data.frame(fisher_matrix_xxy[upper.tri(fisher_matrix_xxy, diag = FALSE)])
names(upper_fisher_XXY)[1] ="edge_strength"
upper_fisher_XXY <- upper_fisher_XXY %>% 
  mutate(group = "XXY")

upper_fisher_XYY <- as.data.frame(fisher_matrix_xyy[upper.tri(fisher_matrix_xyy, diag = FALSE)])
names(upper_fisher_XYY)[1] ="edge_strength"
upper_fisher_XYY <- upper_fisher_XYY %>% 
  mutate(group = "XYY")

upper_fisher <- full_join(upper_fisher_XXY, upper_fisher_XYY)

# 2: calculate mean and median edge strength per group
upper_fisher %>% 
  group_by(group) %>% 
  summarise(mean(edge_strength))

upper_fisher %>% 
  group_by(group) %>% 
  summarise(median(edge_strength))

# 3: directly test the mean delta Fisher's Z-adjusted pairwise correlation
mean_z_XXY <- upper_fisher %>% 
  filter(group == "XXY") %>% 
  summarise(mean(edge_strength)) %>% 
  pull('mean(edge_strength)')
mean_z_XYY <- upper_fisher %>% 
  filter(group == "XYY") %>% 
  summarise(mean(edge_strength)) %>% 
  pull('mean(edge_strength)')
abs_delta_z <- abs(mean_z_XXY - mean_z_XYY)

# 4: test for Kurtosis using package "moments" and function "kurtosis"
kurtosis_xxy <- kurtosis(upper_fisher_XXY$edge_strength, na.rm = FALSE)
kurtosis_xyy <- kurtosis(upper_fisher_XYY$edge_strength, na.rm = FALSE)

# compare the distribution of upper triangles of XXY and XYY 
XXYuppertriangle <- upper_fisher_XXY %>% 
  pull(edge_strength)
XYYuppertriangle <- upper_fisher_XYY %>% 
  pull(edge_strength)

res=ks.test(XXYuppertriangle, XYYuppertriangle)
res$statistic # D - max distance of cumulative distribution function (CDF) of XXYuppertriangle vs. XYYuppertriangle

ks_estimate_obs <- res$statistic

## Create Fig. 3

edge_strength_density <- ggplot(data=upper_fisher, aes(x=edge_strength, group=group, colour=group)) +
  geom_density(adjust=1.5, alpha=.5) +
  labs(title = "Edge strength in XXY vs XYY", x = "Edge Strength (Fisher's Z-transformed)",
       y = "Density", fill = "Group",
       subtitle = "Median(XXY) = 0.383, median(XYY) = 0.398, absolute delta z = 0.03, permuted p-value = .682") +
  scale_colour_manual(values = c("#ef8a62", "#67a9cf"))+
  geom_vline(xintercept=0.43, lty=2, size=0.5, colour = "#ef8a62")+
  geom_vline(xintercept=0.46, lty=2, size=0.5, colour = "#67a9cf")+
  annotate("text", x=0.16, y=-0.2, label= "mean(XXY) = 0.43", size=5, colour ="#ef8a62") +
  annotate("text", x=0.73, y=-0.2, label= "mean(XYY) = 0.46", size=5, colour = "#67a9cf") +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14),
        plot.title = element_text(size=16),
        plot.subtitle = element_text(size=14))
edge_strength_density
