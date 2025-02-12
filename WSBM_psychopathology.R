############################################################
############################################################
##################### WSBM, XXY and XYY ####################
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
matrix_xyy <- read.csv("matrix_xyy.csv")

### Load in cluster naming files
cluster_naming <- read.csv("cluster_naming.csv")

# 1: convert dataframes to matrices
rownames(matrix_xxy) <- matrix_xxy[,1]
matrix_xxy <- matrix_xxy %>% 
  select(-1)
matrix_xxy <- as.matrix(matrix_xxy)

rownames(matrix_xyy) <- matrix_xyy[,1]
matrix_xyy <- matrix_xyy %>% 
  select(-1)
matrix_xyy <- as.matrix(matrix_xyy)

# 2: Fisher's Z transform
fisher_z_transform <- function(r) {
  fisherz(r)
} 

diag(matrix_xxy) <- NA
diag(matrix_xyy) <- NA

fisher_matrix_xxy <- matrix(sapply(matrix_xxy, fisher_z_transform), nrow = ncol(matrix_xxy))

fisher_matrix_xyy <- matrix(sapply(matrix_xyy, fisher_z_transform), nrow = ncol(matrix_xyy))

#change the row and column names to scale names
rownames(fisher_matrix_xxy) <-colnames(matrix_xxy)
colnames(fisher_matrix_xxy) <- colnames(matrix_xxy)

rownames(fisher_matrix_xyy) <-colnames(matrix_xyy)
colnames(fisher_matrix_xyy) <- colnames(matrix_xyy)

# Create delta matrix by subtracting XXY and XYY fisher's Z corrected matrices
fisher_delta_matrix <- fisher_matrix_xxy - fisher_matrix_xyy

#===============================================================================
## Identifying sets of scales with coordinated differences in coupling between groups

# Cluster on the delta matrix using Weighted Stochastic Block Modeling (WSBM)
# 1: run WSBM on observed adjacency matrix: Delta
set.seed('123')

diag(fisher_delta_matrix) <- 1
delta_sbm_g_op <- BM_gaussian("SBM_sym", fisher_delta_matrix)

# 2: create Fig. S1
delta_sbm_g_op$estimate()

# 3: find optimum number of clusters 
which.max(delta_sbm_g_op$ICL) # answer = 6

# 4: make cluster assignments
delta_cluster_memberships <- apply(delta_sbm_g_op$memberships[[6]]$Z, 1, function(x) which.max(x))
data.frame(scale=colnames(fisher_delta_matrix), sbm_cluster=delta_cluster_memberships) %>% arrange(sbm_cluster)

# 5: determine block order by mean correlation with other blocks
delta_block.reorder <- c(1:6)[order(colMeans(delta_sbm_g_op$model_parameters[[6]]$mu))]
str(delta_sbm_g_op)
delta_sbm_g_op$model_parameters[[6]]$mu

# 6: make new membership vector with blocks re-ordered by their mean correlation with other blocks
delta_cluster_reordered <- as.factor(delta_cluster_memberships)
levels(delta_cluster_reordered) <- rank(colMeans(delta_sbm_g_op$model_parameters[[6]]$m)) #column mean of pairwise correlation between blocks
delta_cluster_reordered <- as.character(delta_cluster_reordered) %>% as.numeric

# 7: make a dataframe of the clusters and their memberships
delta_cluster_labels <- data.frame(scale = colnames(fisher_delta_matrix), 
                                   sbm_cluster = delta_cluster_memberships, 
                                   delta_sbm_cluster_ordered = factor(delta_cluster_memberships,levels = order(colMeans(delta_sbm_g_op$model_parameters[[6]]$mu)))) %>%
  arrange(delta_sbm_cluster_ordered)

# 8: label the blocks
delta_cluster_labels$names <- factor(delta_cluster_labels$delta_sbm_cluster_ordered)

# 9: name clusters based on clinician input (provided in csv file "cluster_naming")

delta_cluster_naming <- cluster_naming %>% 
  select(scale, Delta)
delta_cluster_labels <- full_join(delta_cluster_labels, delta_cluster_naming, by = "scale")

# 10: ensure that the clusters are a factor ordered by mean connectivity to other clusters
delta_cluster_labels$Delta <- factor(delta_cluster_labels$Delta, levels = c("Externalizing/Dissociality","Neurodevelopmental differences-1 (attn/EF/social)", 
                                                                                      "Neurodevelopmental differences-2 (motor/social)", "ASD related", "Externalizing","Internalizing"))

arrange(delta_cluster_labels, delta_cluster_labels$names)

# 11: prepare data for WSBM heatmap

delta_cluster_reordered <- as.numeric(as.character(delta_cluster_reordered))
delta_cluster_order <- as.numeric(as.character(delta_cluster_labels$sbm_cluster))

# 11A: flip the order of the columns and the scales within the columns
delta_cluster_order_cols <- rev(as.numeric(as.character(delta_cluster_labels$sbm_cluster)))

delta_cluster_labels_cols <- delta_cluster_labels %>% 
  arrange(desc(Delta)) %>% 
  group_by(Delta) %>%
  mutate(order_scale = row_number()) %>% 
  ungroup()

delta_cluster_labels_cols <- delta_cluster_labels_cols %>% 
  arrange(desc(Delta), desc(order_scale))

print(delta_cluster_labels_cols)

custom_cols_delta <- delta_cluster_labels_cols$scale
custom_rows_delta <- delta_cluster_labels$scale

# 11B: Get column indices
col_indices_delta <- match(custom_cols_delta, colnames(fisher_delta_matrix))
row_indices_delta <- match(custom_rows_delta, rownames(fisher_delta_matrix))

diag(fisher_delta_matrix) <- 0

# 12: find range for heat map
min(fisher_delta_matrix)
max(fisher_delta_matrix)

# 13: change the order of the rows and columns based on above indices
fisher_delta_matrix_ordered <- fisher_delta_matrix[row_indices_delta,col_indices_delta]

## Create Fig. 4A

superheat(fisher_delta_matrix_ordered, membership.rows = delta_cluster_labels$Delta, 
          membership.cols = delta_cluster_labels_cols$Delta, 
          title = 'WSBM clustering delta', 
          bottom.label = "variable", bottom.label.text.size = 7, 
          left.label = "none",
          bottom.label.text.angle = 90, 
          bottom.label.col = "white",
          bottom.label.text.col = c("black","black","black","black",
                                    
                                    "black","black","black","black",
                                    
                                    "#e31a1c","#e31a1c","#e31a1c","#e31a1c",
                                    "#33a02c", "#1f78b4",
                                    
                                    "black","#a65628","#a65628","#a65628",
                                    "#a65628","#a65628","#a65628","#35978f",
                                    "#fb9a99","#fb9a99","#fb9a99","#ff7f00",
                                    "#e31a1c",
                                    
                                    "#a65628","#fb9a99","#ff7f00","#ff7f00",
                                    "#ff7f00",
                                    
                                    
                                    "black","#35978f","#35978f","#e31a1c",
                                    "#e31a1c","#e31a1c",
                                    
                                    "black","#35978f","#35978f","#35978f",
                                    "#e31a1c","#e31a1c","#e31a1c", "#984ea3",
                                    "#33a02c","#33a02c","#33a02c","#33a02c",
                                    "#1f78b4","#1f78b4","#1f78b4" 
                                    
          ),
          heat.lim = c(-0.6590713,0.5548819),
          heat.pal = c("#67a9cf", "#f7f7f7", "#ef8a62"),
          legend.height = 0.2,
          legend.width = 3,
          legend.text.size = 16) # assign colors to values

# note: cluster names and greyed out boxes were added outside of R during figure creation
#===============================================================================

## Create Fig. 4B

# 1: make an adjacency matrix with the clusters from the SBM
delta_sbm_g_op_adj <- delta_sbm_g_op$model_parameters[[6]]$mu

# 2: order the columns and rows by their mean correlation with other blocks
delta_adj_ordered <- delta_sbm_g_op_adj[delta_block.reorder, delta_block.reorder]

colnames(delta_adj_ordered) <- levels(delta_cluster_labels$Delta)
rownames(delta_adj_ordered) <- levels(delta_cluster_labels$Delta)

# 3: use the graph_from_adjacency_matrix function from the package igraph to make the network structure

delta_ordered_grph <- graph_from_adjacency_matrix(delta_adj_ordered, mode="upper", weighted=T, diag=F, add.colnames=T)

V(delta_ordered_grph)$label <- colnames(delta_adj_ordered)
V(delta_ordered_grph)$label.dist <- 1.5
V(delta_ordered_grph)$label.cex <- plotrix::rescale(colMeans(delta_adj_ordered), c(1, 3))
V(delta_ordered_grph)$label.color <- "black"
V(delta_ordered_grph)$size <- plotrix::rescale(colMeans(delta_adj_ordered), c(7, 17))

E(delta_ordered_grph)$width <- round(abs(get.edge.attribute(delta_ordered_grph)$weight)*20)
E(delta_ordered_grph)$weight <- round(abs(get.edge.attribute(delta_ordered_grph)$weight)*20)

E(delta_ordered_grph)$width <- plotrix::rescale(get.edge.attribute(delta_ordered_grph)$weight, c(5,20))
E(delta_ordered_grph)$weight <- plotrix::rescale(get.edge.attribute(delta_ordered_grph)$weight, c(5,20))


plot(delta_ordered_grph,
     vertex.frame.color = c("#67a9cf", "black", "black",
                            "black", "black","black"),
     vertex.frame.width = c(8,1,1,1,1,1),
     edge.lty=c("dashed","dashed","dashed","dashed","dashed",
                "solid", "dashed", "dashed", "dashed","dashed","dashed",
                "dashed","solid","solid","solid"),
     edge.color = c("grey","grey","grey","grey","grey","#67a9cf","grey","grey","grey",
                    "grey", "grey","grey", "#ef8a62", "#ef8a62","#ef8a62"),
     vertex.color = c("grey", "grey", "grey",
                      "grey", "grey","grey"))

#===============================================================================

## Calculate blockwise averages 

# 1: create function 
calc_blockwise_avgedge <- function(fisher_matrix,cluster_labels) 
{ 
  cluster_assignment <- sort(unique(cluster_labels$sbm_cluster))
  blockwise_avgedge <- matrix(rep(NaN,max(cluster_assignment)*max(cluster_assignment)),
                              ncol = max(cluster_assignment))
  rownames(blockwise_avgedge) <- cluster_assignment
  colnames(blockwise_avgedge) <- cluster_assignment
  for (m in 1:max(cluster_assignment)){
    for (n in 1:m){
      scalename1 <- cluster_labels$scale[cluster_labels$sbm_cluster==cluster_assignment[m]]
      scalename2 <- cluster_labels$scale[cluster_labels$sbm_cluster==cluster_assignment[n]]
      blockvalue <- fisher_matrix[scalename1,scalename2]
      if (m==n){
        blockwise_avgedge[m,n] <- mean(blockvalue[upper.tri(blockvalue, diag = FALSE)])}
      else {blockwise_avgedge[m,n] <- mean(blockvalue)}
    }
  }
  return(data.frame(ClusterAssignment=cluster_assignment,BlockwiseAvgedge=blockwise_avgedge))
}

# 2: apply blockwise average function to delta
delta_blockwise_avgedge_obs <- calc_blockwise_avgedge(fisher_delta_matrix,delta_cluster_labels) %>% 
  select(-1) %>% 
  as.matrix()

colnames(delta_blockwise_avgedge_obs) <- c("Neurodevelopmental differences-1 (attn/EF/social)",
                                           "Neurodevelopmental differences-2 (motor/social)",
                                           "Externalizing", "Internalizing","Externalizing/Dissociality",
                                           "ASD related")
rownames(delta_blockwise_avgedge_obs) <- c("Neurodevelopmental differences-1 (attn/EF/social)",
                                           "Neurodevelopmental differences-2 (motor/social)",
                                           "Externalizing", "Internalizing","Externalizing/Dissociality",
                                           "ASD related")

#===============================================================================

## Create Fig. S2, Panel A

# Run WSBM on observed adj matrix: XXY
set.seed('123') 
diag(matrix_xxy) <- 1
xxy_sbm_g_op <- BM_gaussian("SBM_sym", matrix_xxy) 
xxy_sbm_g_op$estimate()

# 1: find optimum number of clusters 
which.max(xxy_sbm_g_op$ICL) # answer = 7
xxy_sbm_g_op$model_parameters[[7]]$mu

# 2: make cluster assignments
xxy_cluster_memberships <- apply(xxy_sbm_g_op$memberships[[7]]$Z, 1, function(x) which.max(x)) 
data.frame(scale=colnames(matrix_xxy), sbm_cluster=xxy_cluster_memberships) %>% arrange(sbm_cluster) #view as dataframe

# 3: determine block order by mean correlation with other blocks (column mean of pairwise correlation between blocks)
xxy_block.reorder <- c(1:7)[order(colMeans(xxy_sbm_g_op$model_parameters[[7]]$mu))] #mu is the block average correlation

# 4: make new membership vector with blocks re-ordered by their mean correlation with other blocks
xxy_cluster_reordered <- as.factor(xxy_cluster_memberships)
levels(xxy_cluster_reordered) <- xxy_block.reorder # reorder the factor variable based on the correlation

# 5: make df of the clusters and their memberships
xxy_cluster_labels <- data.frame(scale = colnames(c.t.h.XXY[,-c(1:3, 8:11, 29:31, 36:41, 70:71)]), 
                                 sbm_cluster = xxy_cluster_memberships, 
                                 xxy_sbm_cluster_ordered = factor(xxy_cluster_memberships,levels = order(colMeans(xxy_sbm_g_op$model_parameters[[7]]$mu)))) %>%
  arrange(xxy_sbm_cluster_ordered) # 3 6 2 1 4 5 7

# 6: label the blocks
xxy_cluster_labels$names <- factor(xxy_cluster_labels$xxy_sbm_cluster_ordered)

# 7: name clusters
XXY_cluster_naming <- cluster_naming %>% 
  select(scale, XXY)
xxy_cluster_labels <- full_join(xxy_cluster_labels, XXY_cluster_naming, by = "scale")

xxy_cluster_labels$XXY <- factor(xxy_cluster_labels$XXY, levels = c("Dissociality-1","Dissociality-2",
                                                                              "ASD related-2", "Internalizing",
                                                                              "ADHD-related","Externalizing",
                                                                              "ASD related-1"
))

arrange(xxy_cluster_labels, xxy_cluster_labels$names)

xxy_cluster_reordered <- as.numeric(as.character(xxy_cluster_reordered))
xxy_cluster_order <- as.numeric(as.character(xxy_cluster_labels$sbm_cluster))

# 8: find range for heat map
min(matrix_xxy)
max(matrix_xxy)

# 9: flip the order of the columns
xxy_cluster_order_cols <- rev(as.numeric(as.character(xxy_cluster_labels$sbm_cluster)))

xxy_cluster_order_cols <- xxy_cluster_labels %>% 
  arrange(desc(XXY)) %>% 
  group_by(XXY) %>%
  mutate(order_scale = row_number()) %>% 
  ungroup()

xxy_cluster_order_cols <- xxy_cluster_order_cols %>% 
  arrange(desc(XXY), desc(order_scale))

print(xxy_cluster_order_cols)

custom_cols_xxy <- xxy_cluster_order_cols$scale
custom_rows_xxy <- xxy_cluster_labels$scale

# 10: Get column indices
col_indices_xxy <- match(custom_cols_xxy, colnames(matrix_xxy))
row_indices_xxy <- match(custom_rows_xxy, rownames(matrix_xxy))

diag(matrix_xxy) <- 1

# 11: order the matrix
matrix_xxy_ordered <- matrix_xxy[row_indices_xxy,col_indices_xxy]

## Plot the XXY heatmap
superheat(matrix_xxy_ordered, membership.rows = xxy_cluster_labels$XXY, 
          membership.cols = xxy_cluster_order_cols$XXY, # remove left.label = "variable" to show cluster names
          title = 'WSBM clustering XXY', 
          bottom.label = "variable", bottom.label.text.size = 7, 
          left.label = "none",
          bottom.label.text.angle = 90, 
          bottom.label.text.col = c("black","black","black","#a65628",
                                    "#a65628","#a65628",
                                    
                                    "#a65628","#a65628",
                                    
                                    "#35978f","#e31a1c","black","black",
                                    "black","#35978f","#e31a1c","#e31a1c",
                                    "#e31a1c",
                                    
                                    "#984ea3","#33a02c","#33a02c","#33a02c","#33a02c",
                                    "black","#35978f",
                                    
                                    "#e31a1c","#e31a1c","#e31a1c","#e31a1c",
                                    "#e31a1c","black","black",
                                    
                                    "black","black","#a65628","#35978f",
                                    "#35978f","#a65628","#fb9a99","#fb9a99",
                                    "#fb9a99","#fb9a99","#ff7f00","#ff7f00",
                                    
                                    
                                    "#ff7f00","#ff7f00","#35978f","#33a02c",
                                    "#e31a1c","#e31a1c","#1f78b4","#1f78b4","#1f78b4",
                                    "#1f78b4"
          ),
          bottom.label.col = "white",
          heat.lim = c(-0.2024269,1),
          heat.pal.values= c(0,1), 
          heat.pal = c("white", "#ef8a62"),# assign colors to values
          legend.breaks = c(-0.2, 0, 0.5, 1),
          legend.height = 0.2,
          legend.width = 3,
          legend.text.size = 16) 
# note: cluster names were added outside of R during figure creation
#===============================================================================

## Create Fig. S2, Panel B

# Run WSBM on observed adj matrix: XYY
set.seed('123')
diag(matrix_xyy) <- 1
xyy_sbm_g_op <- BM_gaussian("SBM_sym", matrix_xyy)
xyy_sbm_g_op$estimate()

# 1: find optimum number of clusters 
which.max(xyy_sbm_g_op$ICL) # answer = 9

# 2: make cluster assignments
xyy_cluster_memberships <- apply(xyy_sbm_g_op$memberships[[9]]$Z, 1, function(x) which.max(x))
data.frame(scale=colnames(matrix_xyy), sbm_cluster=xyy_cluster_memberships) %>% arrange(sbm_cluster)

# 3: determine block order by mean correlation with other blocks
xyy_block.reorder <- c(1:9)[order(colMeans(xyy_sbm_g_op$model_parameters[[9]]$mu))]
str(xyy_sbm_g_op)

# 4: make new membership vector with blocks re-ordered by their mean correlation with other blocks
xyy_cluster_reordered <- as.factor(xyy_cluster_memberships)
levels(xyy_cluster_reordered) <- xyy_block.reorder 

# 5: make dataframe of the clusters and their memberships
xyy_cluster_labels <- data.frame(scale = colnames(c.t.h.XYY[,-c(1:3, 8:11, 29:31, 36:41, 70:71)]), 
                                 sbm_cluster = xyy_cluster_memberships, 
                                 xyy_sbm_cluster_ordered = factor(xyy_cluster_memberships,levels = order(colMeans(xyy_sbm_g_op$model_parameters[[9]]$mu)))) %>%
  arrange(xyy_sbm_cluster_ordered) # 1 2 4 8 6 9 7 3 5

# 6: label the blocks
xyy_cluster_labels$names <- factor(xyy_cluster_labels$xyy_sbm_cluster_ordered)

# 7: name clusters

XYY_cluster_naming <- cluster_naming %>% 
  select(scale, XYY)
xyy_cluster_labels <- full_join(xyy_cluster_labels, XYY_cluster_naming, by = "scale")

xyy_cluster_labels$XYY <- factor(xyy_cluster_labels$XYY, levels = c("ASD related-1","Dissociality", "Internalizing",
                                                                                    "ASD related-3", "Hyperactivity/Impulsivity", "Externalizing",
                                                                                    "ASD related-2", "Inattention", "Total"))

arrange(xyy_cluster_labels, xyy_cluster_labels$names)

xyy_cluster_reordered <- as.numeric(as.character(xyy_cluster_reordered))
xyy_cluster_order <- as.numeric(as.character(xyy_cluster_labels$sbm_cluster))

# 8: flip the order of the columns
xyy_cluster_order_cols <- rev(as.numeric(as.character(xyy_cluster_labels$sbm_cluster)))

xyy_cluster_order_cols <- xyy_cluster_labels %>% 
  arrange(desc(XXY_cluster_naming)) %>% 
  group_by(XXY_cluster_naming) %>%
  mutate(order_scale = row_number()) %>% 
  ungroup()

xyy_cluster_order_cols <- xyy_cluster_order_cols %>% 
  arrange(desc(XXY_cluster_naming), desc(order_scale))

print(xyy_cluster_order_cols)

custom_cols_xyy <- xyy_cluster_order_cols$scale
custom_rows_xyy <- xyy_cluster_labels$scale

# 9: get column indices
col_indices_xyy <- match(custom_cols_xyy, colnames(matrix_xyy))
row_indices_xyy <- match(custom_rows_xyy, rownames(matrix_xyy))

diag(matrix_xyy) <- 1

# 11: order the matrix
matrix_xyy_ordered <- matrix_xyy[row_indices_xyy,col_indices_xyy]

# 12: find range for heat map
min(matrix_xyy)

## Plot the XYY heatmap
superheat(matrix_xyy_ordered, membership.rows = xyy_cluster_labels$XXY_cluster_naming, 
          membership.cols = xyy_cluster_order_cols$XXY_cluster_naming, 
          title = 'WSBM clustering XYY', 
          bottom.label = "variable", bottom.label.text.size = 7, 
          bottom.label.col = "white",
          bottom.label.text.col = c("black","#35978f","#e31a1c","#33a02c",
                                    
                                    
                                    "black","black","#35978f","#e31a1c",
                                    "#e31a1c",
                                    
                                    "#a65628","#a65628","#a65628",
                                    "#a65628","#a65628","#a65628",
                                    
                                    
                                    "black","black","black","#35978f","#e31a1c",
                                    "#e31a1c","#e31a1c","#984ea3","#33a02c",
                                    "#33a02c","#33a02c","#1f78b4",
                                    
                                    
                                    "#e31a1c","#e31a1c","#1f78b4",
                                    
                                    "black","#35978f","#35978f","#ff7f00",
                                    "#ff7f00","#ff7f00","#e31a1c",
                                    
                                    "black","black","black","#35978f","#e31a1c",
                                    
                                    "#33a02c","#1f78b4","#1f78b4",
                                    
                                    "black","#a65628","#fb9a99","#fb9a99",
                                    "#fb9a99","#fb9a99","#ff7f00","#e31a1c"
                                    
          ),
          left.label = "none",
          bottom.label.text.angle = 90, 
          heat.lim = c(-0.2280801,1),
          heat.pal.values= c(0,1), 
          legend.breaks = c(-0.2, 0, 0.5, 1),
          heat.pal = c("white", "#67a9cf"),
          legend.height = 0.2,
          legend.width = 3,
          legend.text.size = 16) # assign colors to values
#===============================================================================

## Create Fig. S2, Panel C

## Create the Sankey diagram

# Sankey diagrams ordered by mean cluster connectivity (column means)
# XXY order: 3 6 2 1 4 5 7
# XYY order: 1 2 4 8 7 9 6 3 5

# 1: manually create the links between clusters and scales across groups
links <- data.frame(
  source=c("Dissociality-1_xxy","Dissociality-1_xxy", "Dissociality-1_xxy", "Dissociality-1_xxy", 
           "Dissociality-1_xxy", "Dissociality-1_xxy",
           
           "Dissociality-2_xxy","Dissociality-2_xxy",
           
           "ASD related-2_xxy","ASD related-2_xxy","ASD related-2_xxy","ASD related-2_xxy",
           "ASD related-2_xxy","ASD related-2_xxy","ASD related-2_xxy","ASD related-2_xxy","ASD related-2_xxy",
           
           "Internalizing_xxy","Internalizing_xxy","Internalizing_xxy","Internalizing_xxy","Internalizing_xxy",
           "Internalizing_xxy","Internalizing_xxy",
           
           "ADHD-related_xxy","ADHD-related_xxy","ADHD-related_xxy","ADHD-related_xxy",
           "ADHD-related_xxy","ADHD-related_xxy","ADHD-related_xxy",
           
           "Externalizing_xxy","Externalizing_xxy","Externalizing_xxy","Externalizing_xxy",
           "Externalizing_xxy","Externalizing_xxy","Externalizing_xxy","Externalizing_xxy",
           "Externalizing_xxy","Externalizing_xxy","Externalizing_xxy","Externalizing_xxy",
           
           "ASD related-1_xxy","ASD related-1_xxy","ASD related-1_xxy","ASD related-1_xxy",
           "ASD related-1_xxy","ASD related-1_xxy","ASD related-1_xxy","ASD related-1_xxy",
           "ASD related-1_xxy","ASD related-1_xxy",
           
           
           "learn_CON", "cont_DCD", "tot_SCQ", "soc_SCQ", "com_SCQ", "rep_SCQ", "mot_SRS", "wt.dep_CBCL",
           
           "cu_APSD", "narc_APSD", "ag.p_SHRP",
           
           "exec_CON","emo_SDQ", "ax.dep_CBCL","som_CBCL","int_CBCL",
           
           "peer_CON","fine_DCD", "coord_DCD","tot_DCD","peer_SDQ","soc_SDQ","tht_CBCL",
           
           "imp_APSD","hyp.imp_CON", "adhd.hyp_CON",
           
           "tot_APSD","ag.v_SHRP", "bul_SHRP","ag.c_SHRP","tot_ARI","def_CON",
           "cond_CON","odd_CON", "cond_SDQ","rule_CBCL","ag_CBCL","ext_CBCL",
           
           "tot_SRS","awr_SRS", "cog_SRS","com_SRS","man_SRS", "asd_SRS",
           
           "att_CON","adhd.att_CON", "hyp.imp_SDQ","soc_CBCL","att_CBCL",
           
           "hos_SHRP","tot_CON", "tot_SDQ","tot_CBCL"
           
  ), 
  target=c("tot_APSD","cu_APSD","narc_APSD","imp_APSD","learn_CON","exec_CON",
           
           "ag.p_SHRP","soc_SDQ",
           
           "cont_DCD","fine_DCD","coord_DCD","tot_DCD","tot_SCQ",
           "soc_SCQ", "com_SCQ", "rep_SCQ", "awr_SRS",
           
           "emo_SDQ","peer_SDQ","mot_SRS","ax.dep_CBCL","wt.dep_CBCL",
           "som_CBCL","int_CBCL",
           
           "att_CON", "hyp.imp_CON","tot_CON","adhd.att_CON","adhd.hyp_CON",
           "hyp.imp_SDQ","att_CBCL",
           
           "ag.v_SHRP","bul_SHRP","ag.c_SHRP","hos_SHRP","tot_ARI","def_CON",
           "cond_CON","odd_CON","cond_SDQ","rule_CBCL","ag_CBCL","ext_CBCL",
           
           "peer_CON", "tot_SDQ", "tot_SRS", "cog_SRS","com_SRS", "man_SRS", "asd_SRS", 
           "soc_CBCL","tht_CBCL", "tot_CBCL",
           
           
           "ASD related-1_xyy","ASD related-1_xyy", "ASD related-1_xyy", "ASD related-1_xyy","ASD related-1_xyy","ASD related-1_xyy", "ASD related-1_xyy", "ASD related-1_xyy",
           "Dissociality_xyy","Dissociality_xyy", "Dissociality_xyy",
           "Internalizing_xyy","Internalizing_xyy", "Internalizing_xyy","Internalizing_xyy","Internalizing_xyy",
           "ASD related-3_xyy","ASD related-3_xyy", "ASD related-3_xyy","ASD related-3_xyy","ASD related-3_xyy","ASD related-3_xyy","ASD related-3_xyy",
           "Hyperactivity/Impulsivity_xyy","Hyperactivity/Impulsivity_xyy", "Hyperactivity/Impulsivity_xyy",
           "Externalizing_xyy","Externalizing_xyy", "Externalizing_xyy","Externalizing_xyy","Externalizing_xyy","Externalizing_xyy",
           "Externalizing_xyy","Externalizing_xyy", "Externalizing_xyy","Externalizing_xyy","Externalizing_xyy","Externalizing_xyy",
           "ASD related-2_xyy","ASD related-2_xyy", "ASD related-2_xyy","ASD related-2_xyy","ASD related-2_xyy", "ASD related-2_xyy",
           "Inattention_xyy","Inattention_xyy", "Inattention_xyy","Inattention_xyy","Inattention_xyy",
           "Total_xyy","Total_xyy", "Total_xyy","Total_xyy"
  ), 
  value=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
          1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
          1,1,1,
          1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
          1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
          1,1,1)
)

# 2: From these flows, we need to create a node data frame: it lists every entity involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# 3: With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# 4: Visualize the Network

p  <- sankeyNetwork(Links = links, Nodes = nodes,
                    Source = "IDsource", Target = "IDtarget",
                    Value = "value", NodeID = "name", 
                    sinksRight=FALSE,
                    iterations = 0,
                    fontSize = 14,
                    nodeWidth = 30) #iterations = 0 preserves the order of the nodes
p

# note: Sankey diagram was reflected vertically and color coding was edited outside of R
