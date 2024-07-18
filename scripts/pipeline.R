# set your work dir
set.wd("your/path/")

# import functions
source("utils/internal_functions.R")
source("utils/scib_calculate_scores.R")
source("utils/scib_score2tab.R")
source("utils/scib_tab2plot.R")
source("utils/scib_NicePlot.R")

# step 1: calculate scores from scib output
scib_bm_table <- read.csv("data/scib_metrics_output.csv", row.names = 1)
scib_bm_table <- t(scib_bm_table) %>% as.data.frame()
 
# preview
#                                  CCA    Harmony       scVI     scANVI     Combat  Scanorama    FastMNN MIRA_feature MIRA_topic      BBKNN   CellHint
# NMI_cluster/label         0.68113424 0.40519090 0.70671790 0.78339236 0.67791204 0.07367235 0.70544284   0.64829393  0.6195890 0.63326321 0.63792718
# ARI_cluster/label         0.42592147 0.18460916 0.45410402 0.60476635 0.44332553 0.01662426 0.45040404   0.38200270  0.3622267 0.38102480 0.37931134
# ASW_label                 0.52110099 0.39779555 0.53132162 0.57205420 0.50618038 0.45541007 0.54055252   0.47874877  0.4217216 0.55051497 0.55051497
# ASW_label/batch           0.76473195 0.74994410 0.81212279 0.81441418 0.87448696 0.76810873 0.81290315   0.63291911  0.4328971 0.83497916 0.83497916
# PCR_batch                 0.59717657 0.98060867 0.72960228 0.65156498 0.97685082 0.36158519 0.58487625   0.31726962  0.2983176 0.39905948 0.39905948
# cell_cycle_conservation   0.52890935 0.62995832 0.64323095 0.60428295 0.45395620 0.33703927 0.48226669   0.41737277  0.3747284 0.72444419 0.72444419
# isolated_label_F1         0.07960199 0.09638554 0.07729469 0.07881773 0.18750000 0.02777778 0.08080808   0.06698565  0.1750000 0.08121827 0.08121827
# isolated_label_silhouette 0.68413717 0.71319951 0.75965893 0.75051224 0.66156498 0.47277215 0.75466657   0.61397110  0.6961596 0.69666968 0.69666968
# graph_conn                0.91086152 0.69119826 0.94227532 0.94778690 0.90399530 0.08817273 0.92683494   0.73696937  0.5444755 0.99941979 0.99987096
# kBET                      0.83627548 0.69261077 0.79482423 0.79398228 0.66163020 0.03976692 0.82879844   0.75732297  0.6348482 0.74590008 0.74406771
# iLISI                     0.11064468 0.13437375 0.09931888 0.09130411 0.08422109 0.08968606 0.09891523   0.10286301  0.1113102 0.50229592 0.44128234
# cLISI                     0.98150708 0.94894612 0.98188344 0.99171863 0.97944120 0.87262356 0.98292546   0.97101663  0.9679921 0.81405685 0.90387700
# hvg_overlap                       NA         NA         NA         NA         NA         NA         NA           NA         NA         NA         NA
# trajectory                        NA         NA         NA         NA         NA         NA         NA           NA         NA         NA         NA


scib_scores <- scib_calculate_scores(scib_bm_table)

# preview
#              ... graph_conn       kBET Batch_Correction Bio_conservation Overall_Score
# scANVI       ... 0.94778690 0.79398228        0.6598105        0.6265064     0.6398280
# scVI         ... 0.94227532 0.79482423        0.6756287        0.5934588     0.6263268
# Combat       ... 0.90399530 0.66163020        0.7002369        0.5585543     0.6152273
# CellHint     ... 0.99987096 0.74406771        0.6838519        0.5677089     0.6141661
# BBKNN        ... 0.99941979 0.74590008        0.6963309        0.5544560     0.6112060
# FastMNN      ... 0.92683494 0.82879844        0.6504656        0.5710095     0.6027919
# CCA          ... 0.91086152 0.83627548        0.6439380        0.5574732     0.5920591
# Harmony      ... 0.69119826 0.69261077        0.6497471        0.4822979     0.5492776
# MIRA_feature ... 0.73696937 0.75732297        0.5094688        0.5111988     0.5105068
# MIRA_topic   ... 0.54447549 0.63484822        0.4043697        0.5167739     0.4718122
# Scanorama    ... 0.08817273 0.03976692        0.2694639        0.3222742     0.3011501


# Prepare these data manually
scib_summary_tab <- scib_scores %>% 
  dplyr::select(Batch_Correction, Bio_conservation, Overall_Score) %>% 
  mutate(Methods = rownames(.)) %>% 
  mutate(Features = c("HVG","HVG","FULL", rep("HVG", 8))) %>% 
  mutate(Scaling = c("unscaled", "unscaled", "unscaled", "scaled", "unscaled", "scaled", 
                     "unscaled", "scaled", "unscaled", "unscaled", "scaled")) %>% 
  mutate(Output = c("Embedding", "Embedding", "Features", "Graph", "Graph","Features", 
                    "Features", "Embedding", "Embedding", "Embedding", "Embedding"))

# preview
#              Batch_Correction Bio_conservation Overall_Score      Methods Features  Scaling    Output
# scANVI              0.6598105        0.6265064     0.6398280       scANVI      HVG unscaled Embedding
# scVI                0.6756287        0.5934588     0.6263268         scVI      HVG unscaled Embedding
# Combat              0.7002369        0.5585543     0.6152273       Combat     FULL unscaled  Features
# CellHint            0.6838519        0.5677089     0.6141661     CellHint      HVG   scaled     Graph
# BBKNN               0.6963309        0.5544560     0.6112060        BBKNN      HVG unscaled     Graph
# FastMNN             0.6504656        0.5710095     0.6027919      FastMNN      HVG   scaled  Features
# CCA                 0.6439380        0.5574732     0.5920591          CCA      HVG unscaled  Features
# Harmony             0.6497471        0.4822979     0.5492776      Harmony      HVG   scaled Embedding
# MIRA_feature        0.5094688        0.5111988     0.5105068 MIRA_feature      HVG unscaled Embedding
# MIRA_topic          0.4043697        0.5167739     0.4718122   MIRA_topic      HVG unscaled Embedding
# Scanorama           0.2694639        0.3222742     0.3011501    Scanorama      HVG   scaled Embedding


preprocessed_data <- scib_preprocessing_tab(scib_summary_tab)
# [1] "All required columns are present."
# [1] "Values in Output are corrected."
# [1] "Values in Scaling are corrected."
# [1] "Values in Features are corrected."
# [1] "Data preprocessing has been completed..."

# preview
#                   Methods Features  Scaling overall_rank bio_rank batch_rank    Output output_img Overall_Score Bio_conservation Batch_Correction
# scANVI             scANVI      HVG unscaled          1.0      1.0        0.6 Embedding      embed     0.6398280        0.6265064        0.6598105
# scVI                 scVI      HVG unscaled          0.9      0.9        0.7 Embedding      embed     0.6263268        0.5934588        0.6756287
# Combat             Combat     FULL unscaled          0.8      0.6        1.0  Features       gene     0.6152273        0.5585543        0.7002369
# CellHint         CellHint      HVG   scaled          0.7      0.7        0.8     Graph      graph     0.6141661        0.5677089        0.6838519
# BBKNN               BBKNN      HVG unscaled          0.6      0.4        0.9     Graph      graph     0.6112060        0.5544560        0.6963309
# FastMNN           FastMNN      HVG   scaled          0.5      0.8        0.5  Features       gene     0.6027919        0.5710095        0.6504656
# CCA                   CCA      HVG unscaled          0.4      0.5        0.3  Features       gene     0.5920591        0.5574732        0.6439380
# Harmony           Harmony      HVG   scaled          0.3      0.1        0.4 Embedding      embed     0.5492776        0.4822979        0.6497471
# MIRA_feature MIRA_feature      HVG unscaled          0.2      0.2        0.2 Embedding      embed     0.5105068        0.5111988        0.5094688
# MIRA_topic     MIRA_topic      HVG unscaled          0.1      0.3        0.1 Embedding      embed     0.4718122        0.5167739        0.4043697
# Scanorama       Scanorama      HVG   scaled          0.0      0.0        0.0 Embedding      embed     0.3011501        0.3222742        0.2694639



prepared_data_list <- scib_tab2plot(preprocessed_data, label_top3 = T, Image_path = "./images")
# [1] "Return a list object containing plot_data."
# Warning message:
# Setting row names on a tibble is deprecated.


plot <- scib_nice_plot(plot_data = prepared_data_list, label_top3 = T)
# The final plot is showed in README.









