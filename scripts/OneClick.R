# example for scib-functions
rm(list = ls())

source("utils/scib/scib_utils.R")

# input: raw output from scib
scib_bm_table <- read.csv("data/13_scib_metrics_concat2.csv", row.names = 1)
scib_bm_table <- t(scib_bm_table) %>% as.data.frame()

# calculate scores
scib_scores <- scib_calculate_scores(scib_bm_table)

# add information manually
scib_summary_tab <- scib_scores %>% 
  dplyr::select(Batch_Correction, Bio_conservation, Overall_Score) %>% 
  mutate(Methods = rownames(.)) %>% 
  mutate(Features = c("HVG","HVG","FULL", rep("HVG", 8))) %>% 
  mutate(Scaling = c("unscaled", "unscaled", "unscaled", "scaled", "unscaled", "scaled", 
                     "unscaled", "scaled", "unscaled", "unscaled", "scaled")) %>% 
  mutate(Output = c("Embedding", "Embedding", "Features", "Graph", "Graph","Features", 
                    "Features", "Embedding", "Embedding", "Embedding", "Embedding"))

# get the plot directly
plot <- scib_OneClick(scib_summary_tab)
