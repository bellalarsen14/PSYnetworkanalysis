############################################################
############################################################
############### Permutation, XXY and XYY ###################
############################################################
############################################################

### SET UP

library(magrittr)
library(dplyr)
library(ggplot2)
library(gplots)
library(purrr)
library(psych)
library(tidyr)
library(stringr)
library(tibble)
library(purrrlyr)
library(reshape2)
library(skimr)
library(cocor)
library(devtools)
library(SimplyAgree)
library(moments)
library(Hmisc)
#library(DescTools)
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("scmap")
library(scmap)
library(tictoc)

### Load in WSBM assignment files (R data file)
load(file="data/WSBM_assignment_git.RData")

#===============================================================================

### Define functions to be used

# 1: Transform correlation coefficients using Fisher's Z
fisher_z_transform <- function(r) {
  fisherz(r)
} 

# 2: Calculate edge wise correlation
calculate_edge_corr <- function(r1,r2){ 
  cor(r1,r2, use = "na.or.complete")
}

# 3: Calculate block-wise avg edge based on WSBM assignments from the observed data
calc_blockwise_avgedge <- function(fisher_matrix,cluster_labels) 
{ 
  cluster_assignment <- sort(unique(cluster_labels$sbm_cluster))
  blockwise_avgedge <- matrix(rep(NaN,max(cluster_assignment)*max(cluster_assignment)),
                              ncol = max(cluster_assignment)) # create empty matrix to hold the average edge values
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

# 4: calculate nodal degree, edge correlation, upper fisher Z, blockwise avg edge, using the null data
calc_metrics <- function(zscore_xxy,zscore_xyy,xxy_cluster_labels,xyy_cluster_labels,delta_cluster_labels)
{
  # create pairwise correlation matrices from the two dataframes
  matrix_xxy <- zscore_xxy %>%
    rcorr(.,) %>%
    extract2(1)
  matrix_xyy <- zscore_xyy %>%
    rcorr(.,) %>%
    extract2(1)
  
  diag(matrix_xxy) <- NA
  diag(matrix_xyy) <- NA
  
  fisher_matrix_xxy <- matrix(sapply(matrix_xxy, fisher_z_transform), nrow = ncol(matrix_xxy)) 
  fisher_matrix_xyy <- matrix(sapply(matrix_xyy, fisher_z_transform), nrow = ncol(matrix_xyy)) 
  
  # change the row and column names to scale names
  row.names(fisher_matrix_xxy) <-colnames(matrix_xxy)
  colnames(fisher_matrix_xxy) <- colnames(matrix_xxy)
  
  row.names(fisher_matrix_xyy) <-colnames(matrix_xyy)
  colnames(fisher_matrix_xyy) <- colnames(matrix_xyy)
  
  ## extract nodal degree for all scales
  
  # calculate the row means for each row
  nodal_xxy <- rowMeans(fisher_matrix_xxy, na.rm = T)
  nodal_xyy <- rowMeans(fisher_matrix_xyy, na.rm = T)
  
  # create df of both disorder row means
  df_nodal_xxy <- as.data.frame(nodal_xxy) %>% 
    mutate(group = "XXY")
  df_nodal_xyy <- as.data.frame(nodal_xyy) %>% 
    mutate(group = "XYY")
  
  df_nodal_xxy <- tibble::rownames_to_column(df_nodal_xxy, "scale_name") %>% 
    select(-group)
  df_nodal_xyy <- tibble::rownames_to_column(df_nodal_xyy, "scale_name") %>% 
    select(-scale_name, -group)
  
  df_nodal_total <- bind_cols(df_nodal_xxy, df_nodal_xyy) 
  
  # 4.1: make a vector of the upper triangle of the Fisher's Z matrices, turn into dataframe for future use
  upper_fisher_XXY <- fisher_matrix_xxy[upper.tri(fisher_matrix_xxy, diag = FALSE)]
  
  upper_fisher_XYY <- fisher_matrix_xyy[upper.tri(fisher_matrix_xyy, diag = FALSE)]
  
  upper_fisher_XXY_df <- as.data.frame(upper_fisher_XXY)
  names(upper_fisher_XXY_df)[1] ="edge_strength"
  upper_fisher_XXY_df <- upper_fisher_XXY_df %>% 
    mutate(group = "XXY")
  
  upper_fisher_XYY_df <- as.data.frame(upper_fisher_XYY)
  names(upper_fisher_XYY_df)[1] ="edge_strength"
  upper_fisher_XYY_df <- upper_fisher_XYY_df %>% 
    mutate(group = "XYY")
  
  upper_fisher <- bind_rows(upper_fisher_XXY_df, upper_fisher_XYY_df)
  
  # 4.2: calculate the edge-wise correlation between XXY and XYY and fingerprint divergence score (1-r)
  
  edge_corr <- map2_dfr(as.data.frame(fisher_matrix_xxy), as.data.frame(fisher_matrix_xyy), calculate_edge_corr)
  edge_corr <- edge_corr %>% 
    pivot_longer(cols = 1:53, names_to = "scale_name", values_to = "correlation") %>% 
    mutate(divergence = 1-correlation)
  edge_corr_r <- setNames(edge_corr$correlation, edge_corr$scale_name)
  edge_corr_div <- setNames(edge_corr$divergence, edge_corr$scale_name)
  
  # 4.3: save the correlations per group for two scales
  diag(matrix_xxy) <- NA
  diag(matrix_xyy) <- NA
  
  ag_p_xxy <- as.data.frame(matrix_xxy) %>% 
    select(ag.p_SHRP) %>% 
    filter(ag.p_SHRP != "ag.p_SHRP")
  ag_p_xxy <- setNames(ag_p_xxy$ag.p_SHRP, rownames(ag_p_xxy))
  ag_p_xyy <- as.data.frame(matrix_xyy) %>% 
    select(ag.p_SHRP) %>% 
    filter(ag.p_SHRP != "ag.p_SHRP")
  ag_p_xyy <- setNames(ag_p_xyy$ag.p_SHRP, rownames(ag_p_xyy))
  
  
  soc_CBCL_xxy <- as.data.frame(matrix_xxy) %>% 
    select(soc_CBCL) %>% 
    filter(soc_CBCL != "soc_CBCL")
  soc_CBCL_xxy <- setNames(soc_CBCL_xxy$soc_CBCL, rownames(soc_CBCL_xxy))
  soc_CBCL_xyy <- as.data.frame(matrix_xyy) %>% 
    select(soc_CBCL) %>% 
    filter(soc_CBCL != "soc_CBCL")
  soc_CBCL_xyy <- setNames(soc_CBCL_xyy$soc_CBCL, rownames(soc_CBCL_xyy))
  
  # 4.4: create delta matrix
  fisher_delta_matrix <- fisher_matrix_xxy - fisher_matrix_xyy
  row.names(fisher_delta_matrix) <-colnames(matrix_xyy)
  colnames(fisher_delta_matrix) <- colnames(matrix_xyy)
  
  # 4.5: calculate block-wise avg edge based on WSBM assignments from the observed data
  delta_blockwise_avgedge <- calc_blockwise_avgedge(fisher_delta_matrix,delta_cluster_labels) %>% 
    select(-1) %>% 
    as.matrix()
  
  # 4.6: directly test the mean delta Z
  mean_z_XXY <- upper_fisher %>% 
    filter(group == "XXY") %>% 
    summarise(mean(edge_strength)) %>% 
    pull('mean(edge_strength)')
  mean_z_XYY <- upper_fisher %>% 
    filter(group == "XYY") %>% 
    summarise(mean(edge_strength)) %>% 
    pull('mean(edge_strength)')
  abs_delta_z <- abs(mean_z_XXY - mean_z_XYY)
  
  # 4.7: run a Two-sample Kolmogorov-Smirnov (KS) test
  ks_res = ks.test(upper_fisher_XXY, upper_fisher_XYY)
  ks_estimate <- ks_res$statistic
  
# 5: Return outputs from (#4)
  return(list(nodal_xxy=nodal_xxy,nodal_xyy=nodal_xyy, upper_fisher_XXY=upper_fisher_XXY,
              upper_fisher_XYY=upper_fisher_XYY,
              edge_corr_r=edge_corr_r,
              edge_corr_div=edge_corr_div,
              ag_p_xxy=ag_p_xxy,
              ag_p_xyy=ag_p_xyy,
              soc_CBCL_xxy=soc_CBCL_xxy,
              soc_CBCL_xyy=soc_CBCL_xyy,
              matrix_xxy=matrix_xxy,
              matrix_xyy=matrix_xyy,
              fisher_delta_matrix=fisher_delta_matrix,
              delta_blockwise_avgedge=delta_blockwise_avgedge,
              abs_delta_z=abs_delta_z,
              ks_estimate=ks_estimate
              )) # add here to return and then enter the calculation above for each
  
}

#===============================================================================
### Prepare permutation

## Read in data
# 1: z-score matrices for the observed z-scores per group
load(file="data/zscore_xxy.RData") 
load(file="data/zscore_xyy.RData")

# 2: observed nodal degree and edgewise correlation and divergence values
nodal_xxy_obs <- read.csv(file ="data/nodal_xxy_obs.csv")
nodal_xyy_obs <- read.csv(file ="data/nodal_xyy_obs.csv")
edge_corr_div_obs <- read.csv(file ="data/edge_corr_divergence_obs.csv")

nodal_xxy_obs <- nodal_xxy_obs %>% 
  column_to_rownames(var = "X") %>% 
  rename(nodal_degree = x)
nodal_xyy_obs <- nodal_xyy_obs %>% 
  column_to_rownames(var = "X") %>% 
  rename(nodal_degree = x)

# 3: observed correlation between all scales and two example scales
ag_p_xxy_obs <- read.csv(file ="data/ag_p_xxy_obs.csv")
colnames(ag_p_xxy_obs)[1] <- "scale_name"
colnames(ag_p_xxy_obs)[2] <- "corr_ag_p"
ag_p_xyy_obs <- read.csv(file ="data/ag_p_xyy_obs.csv")
colnames(ag_p_xyy_obs)[1] <- "scale_name"
colnames(ag_p_xyy_obs)[2] <- "corr_ag_p"
soc_CBCL_xxy_obs <- read.csv(file ="data/soc_CBCL_xxy_obs.csv")
colnames(soc_CBCL_xxy_obs)[1] <- "scale_name"
colnames(soc_CBCL_xxy_obs)[2] <- "corr_soc_CBCL"
soc_CBCL_xyy_obs <- read.csv(file ="data/soc_CBCL_xyy_obs.csv")
colnames(soc_CBCL_xyy_obs)[1] <- "scale_name"
colnames(soc_CBCL_xyy_obs)[2] <- "corr_soc_CBCL"

# 4: absolute value of the observed delta z-score
load(file="data/abs_delta_z_obs.RData")

# 5: observed blockwise average edge values for the delta matrix
load(file="data/delta_blockwise_avgedge_obs.RData")

# 6: load observed KS test results
load(file="data/ks_estimate_obs.RData")

#===============================================================================
### Set up continued
combined_zscore <- rbind(zscore_xxy,zscore_xyy)
group <-  c(rep('XXY',nrow(zscore_xxy)),rep('XYY',nrow(zscore_xyy)))

## 1: set the number of permutations
num_perm <- 10000

## 2: set seed for reproducibility
set.seed(03172020)

## 3: create null storage systems in which to save the permutation values

  # 3.1: nodal degree
nodal_xxy_null <- matrix(rep(NaN, 53*num_perm), ncol=53)
colnames(nodal_xxy_null) <- colnames(combined_zscore)
nodal_xyy_null <- matrix(rep(NaN, 53*num_perm), ncol=53)
colnames(nodal_xyy_null) <- colnames(combined_zscore)

  # 3.2: divergence scores
edge_corr_div_null <- matrix(rep(NaN, 53*num_perm), ncol=53)
colnames(edge_corr_div_null) <- colnames(combined_zscore)

  # 3.3: differences for 2 scales
colnames_ag_p <- combined_zscore[,-9]
ag_p_xxy_null <- matrix(rep(NaN, 52*num_perm), ncol=52)
colnames(ag_p_xxy_null) <- colnames(colnames_ag_p)
ag_p_xyy_null <- matrix(rep(NaN, 52*num_perm), ncol=52)
colnames(ag_p_xyy_null) <- colnames(colnames_ag_p)

colnames_soc_CBCL <- combined_zscore[,-46]
soc_CBCL_xxy_null <- matrix(rep(NaN, 52*num_perm), ncol=52)
colnames(soc_CBCL_xxy_null) <- colnames(colnames_soc_CBCL)
soc_CBCL_xyy_null <- matrix(rep(NaN, 52*num_perm), ncol=52)
colnames(soc_CBCL_xyy_null) <- colnames(colnames_soc_CBCL)

  # 3.4: upper fisher's triangle
upper_fisher_xxy_null <- matrix(rep(NaN, 1378*num_perm), ncol=1378) 
upper_fisher_xyy_null <- matrix(rep(NaN, 1378*num_perm), ncol=1378) 

  # 3.5: delta z
upper_deltaz_null <- rep(NaN, num_perm)

  # 3.6: blockwise avgeage edge from null xxy, xyy, delta matrices 
delta_blockwise_null <- array(numeric(),c(max(delta_cluster_labels$sbm_cluster),max(delta_cluster_labels$sbm_cluster),num_perm))

  # 3.7: xxy and xyy correlation matrices 
xxy_matrix_null <- array(numeric(),c(53,53,num_perm))
xyy_matrix_null <- array(numeric(),c(53,53,num_perm))
fisher_delta_matrix_null <- array(numeric(),c(53,53,num_perm))

  # 3.8: ks results
ks_estimate_null <- rep(NaN, num_perm)

#===============================================================================
### Run permutation

tic() # to calculate the permutation run time

## Loop - creating 10,000 null distributions for each test statistic, shuffling XXY/XYY group assignment
for (n in 1:num_perm)
{
  print(n)
  # A: shuffling the karyotype group label to make a new "group" variable
  randorder <- sample(length(group))
  group_null <- group[randorder]
  # B: calculate z-score
  zscore_xxy_null <- combined_zscore[which(group_null=='XXY'),]
  zscore_xyy_null <- combined_zscore[which(group_null=='XYY'),]
  # 3.1: calculate nodal degree (numbers correspond to above functions)
  metric_list <- calc_metrics(zscore_xxy_null,zscore_xyy_null,
                              xxy_cluster_labels,xyy_cluster_labels,delta_cluster_labels) 
  nodal_xxy_null[n,] <- metric_list$nodal_xxy
  nodal_xyy_null[n,] <- metric_list$nodal_xyy
  # 3.2: pull divergence score
  edge_corr_div_null[n,] <- metric_list$edge_corr_div
  # 3.3: pull differences for 2 scales
  ag_p_xxy_null[n,] <- metric_list$ag_p_xxy
  ag_p_xyy_null[n,] <- metric_list$ag_p_xyy
  soc_CBCL_xxy_null[n,] <- metric_list$soc_CBCL_xxy
  soc_CBCL_xyy_null[n,] <- metric_list$soc_CBCL_xyy
  # 3.4: pull upper fisher
  upper_fisher_xxy_null[n,] <- metric_list$upper_fisher_XXY
  upper_fisher_xyy_null[n,] <- metric_list$upper_fisher_XYY
  # 3.5: pull delta z
  upper_deltaz_null[n] <- metric_list$abs_delta_z
  # 3.6: pull blockwise average edge correlation (block is defined by the observed data)
  delta_blockwise_null[,,n] <- metric_list$delta_blockwise_avgedge
  # 3.7: pull xxy and xyy correlation matrices
  xxy_matrix_null[,,n] <- metric_list$matrix_xxy
  xyy_matrix_null[,,n] <- metric_list$matrix_xyy
  fisher_delta_matrix_null[,,n] <- metric_list$fisher_delta_matrix
  # 3.8: pull ks test estimate
  ks_estimate_null[n] <- metric_list$ks_estimate
}
toc() # approx. 193 seconds (~3min)

#===============================================================================
## Calculate p_perm based the null distribution (numbers correspond to above, 3.1, etc.)

# 3.1 nodal degree
  # XXY
nodal_degree_xxy_p_perm <- matrix(rep(NaN,ncol(nodal_xxy_null)),ncol=ncol(nodal_xxy_null))
colnames(nodal_degree_xxy_p_perm)<- colnames(nodal_xxy_null)
for (m in 1:ncol(nodal_xxy_null))
{
  nodal_degree_xxy_p_perm[1,m] <- length(which(abs(nodal_xxy_null[,m])>abs(nodal_xxy_obs$nodal_degree[m])))/num_perm
}

nodal_degree_xxy_p_perm_long <- data.frame(nodal_degree_xxy_p_perm) %>% 
  pivot_longer(cols = 1:53, names_to = "scale_name", values_to = "p_perm") %>% 
  mutate(p_perm.bonf=p.adjust(p_perm, method="bonferroni")) 

nodal_degree_xxy_p_perm_long$observed <- nodal_xxy_obs$nodal_degree

  # XYY
nodal_degree_xyy_p_perm <- matrix(rep(NaN,ncol(nodal_xyy_null)),ncol=ncol(nodal_xyy_null))
colnames(nodal_degree_xyy_p_perm)<- colnames(nodal_xyy_null)
for (m in 1:ncol(nodal_xyy_null))
{
  nodal_degree_xyy_p_perm[1,m] <- length(which(abs(nodal_xyy_null[,m])>abs(nodal_xyy_obs$nodal_degree[m])))/num_perm
}

nodal_degree_xyy_p_perm_long <- data.frame(nodal_degree_xyy_p_perm) %>% 
  pivot_longer(cols = 1:53, names_to = "scale_name", values_to = "p_perm") %>% 
  mutate(p_perm.bonf=p.adjust(p_perm, method="bonferroni"))

nodal_degree_xyy_p_perm_long$observed <- nodal_xyy_obs$nodal_degree

# 3.2: fingerprint divergence score
divergence_p_perm <- matrix(rep(NaN,ncol(edge_corr_div_null)),ncol=ncol(edge_corr_div_null))
colnames(divergence_p_perm)<- colnames(edge_corr_div_null)
for (m in 1:ncol(edge_corr_div_null))
{
  divergence_p_perm[1,m] <- length(which(abs(edge_corr_div_null[,m])>abs(edge_corr_div_obs$divergence[m])))/num_perm
} 

divergence_p_perm_long <- data.frame(divergence_p_perm) %>% 
  pivot_longer(cols = 1:53, names_to = "scale_name", values_to = "p_perm") %>% 
  mutate(p_perm.bonf=p.adjust(p_perm, method="bonferroni")) 

divergence_p_perm_long$observed <- edge_corr_div_obs$divergence # reproduces column "P" in Table S1

# 3.3: differences for 2 scales

#ag.p_SHRP
ag_p_diff_p_perm <- matrix(rep(NaN,ncol(ag_p_xxy_null)),ncol=ncol(ag_p_xxy_null))
colnames(ag_p_diff_p_perm)<- colnames(ag_p_xxy_null)
for (m in 1:ncol(ag_p_xxy_null))
{
  ag_p_diff_p_perm[1,m] <- length(which(abs(ag_p_xxy_null[,m]-ag_p_xyy_null[,m])
                                        >abs(ag_p_xxy_obs$corr_ag_p[m]-ag_p_xyy_obs$corr_ag_p[m])))/num_perm
}

ag_p_diff_p_perm_long <- data.frame(ag_p_diff_p_perm) %>% 
  pivot_longer(cols = 1:52, names_to = "scale_name", values_to = "p_perm") %>% 
  mutate(p_perm.bonf=p.adjust(p_perm, method="bonferroni"))

ag_p_diff_p_perm_long$observed <- (ag_p_xxy_obs$corr_ag_p)-(ag_p_xyy_obs$corr_ag_p)

#soc_CBCL
soc_CBCL_p_perm <- matrix(rep(NaN,ncol(soc_CBCL_xxy_null)),ncol=ncol(soc_CBCL_xxy_null))
colnames(soc_CBCL_p_perm)<- colnames(soc_CBCL_xxy_null)
for (m in 1:ncol(soc_CBCL_xxy_null))
{
  soc_CBCL_p_perm[1,m] <- length(which(abs(soc_CBCL_xxy_null[,m]-soc_CBCL_xyy_null[,m])
                                       >abs(soc_CBCL_xxy_obs$corr_soc_CBCL[m]-soc_CBCL_xyy_obs$corr_soc_CBCL[m])))/num_perm
}

soc_CBCL_p_perm_long <- data.frame(soc_CBCL_p_perm) %>% 
  pivot_longer(cols = 1:52, names_to = "scale_name", values_to = "p_perm") %>% 
  mutate(p_perm.bonf=p.adjust(p_perm, method="bonferroni"))

soc_CBCL_p_perm_long$observed <- (soc_CBCL_xxy_obs$corr_soc_CBCL)-(soc_CBCL_xyy_obs$corr_soc_CBCL)


# 3.5: compare observed vs. null absolute delta mean z (XXY-XYY)
delta_z_perm <- length(which(abs(upper_deltaz_null)>abs_delta_z))/num_perm

# 3.6: blockwise average edge correlation
delta_blockwise_p_perm <- matrix(rep(NaN,ncol(delta_blockwise_null)*nrow(delta_blockwise_avgedge_obs)),
                                 ncol=ncol(delta_blockwise_null))
colnames(delta_blockwise_p_perm)<- c("cluster_1", "cluster_2","cluster_3","cluster_4","cluster_5","cluster_6")
rownames(delta_blockwise_p_perm)<- c("cluster_1", "cluster_2","cluster_3","cluster_4","cluster_5","cluster_6")

for (m in 1:nrow(delta_blockwise_avgedge_obs))
{
  
  for (n in 1:m)
  {
    delta_blockwise_p_perm[m,n] <- length(which(abs(delta_blockwise_null[m,n,])>abs(delta_blockwise_avgedge_obs[m,n])))/num_perm
  }
}

colnames(delta_blockwise_p_perm)<- c("cluster_1", "cluster_2","cluster_3","cluster_4","cluster_5","cluster_6")
rownames(delta_blockwise_p_perm)<- c("cluster_1", "cluster_2","cluster_3","cluster_4","cluster_5","cluster_6")

delta_blockwise_p_perm_adjusted <- reshape2::melt(delta_blockwise_p_perm, na.rm = TRUE)

delta_blockwise_p_perm_adjusted <- delta_blockwise_p_perm_adjusted %>% 
  mutate(p_adjusted = p.adjust(value, method="bonferroni")) %>% 
  rename(p_value = value) %>% 
  mutate(Var1 = case_when(Var1 == "cluster_1" ~ "Neurodevelopmental differences-1 (attn/EF/social)",
                          Var1 == "cluster_2" ~ "Neurodevelopmental differences-2 (motor/social)", 
                          Var1 == "cluster_3" ~ "Externalizing",
                          Var1 == "cluster_4" ~ "Internalizing",
                          Var1 == "cluster_5" ~ "Externalizing/Dissociality",
                          Var1 == "cluster_6" ~ "ASD related")) %>% 
  mutate(Var2 = case_when(Var2 == "cluster_1" ~ "Neurodevelopmental differences-1 (attn/EF/social)",
                          Var2 == "cluster_2" ~ "Neurodevelopmental differences-2 (motor/social)", 
                          Var2 == "cluster_3" ~ "Externalizing",
                          Var2 == "cluster_4" ~ "Internalizing",
                          Var2 == "cluster_5" ~ "Externalizing/Dissociality",
                          Var2 == "cluster_6" ~ "ASD related"))

# 3.7: KS test permuted p-value
ks_estimate_p_perm <- length(which(ks_estimate_null>ks_estimate_obs))/num_perm

