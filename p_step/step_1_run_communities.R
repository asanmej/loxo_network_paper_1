# Author: Alejandro Santos Mejías
# Date last update: 2026-06-12
# Script name: step_1_run_communities.R
# Aims: to run preparation steps and community detection algorithm


rm(list = ls())
gc()

source("98_packages.R")

path_null_model <- paste0(path_output,"null_model")
null_model_files <- list.files(path_null_model, full.names = TRUE)

path_output <- paste0(path_output, "community_detection")
if (!dir.exists(path_output)) {
  dir.create(path_output, recursive = TRUE)
}

file_list <- list.files(path = path_data,
                       pattern = "(2012|2017|2022)_(H|M)_(60-69|70-79|80-89|90)\\.csv$",
                       full.names = TRUE)

seed <- 29101997
set.seed(seed)

for (file_current in file_list) {

  name_group <- tools::file_path_sans_ext(basename(file_current))
  path_folder_group <- paste0(path_output, "/",gsub("-", "", Sys.Date()))
  if (!dir.exists(path_folder_group)) {
    dir.create(path_folder_group, recursive = TRUE)
  }
  
  cat("Iniciando calculo para", name_group)
  df <- read_csv(file_current, show_col_types = FALSE)
  
  rm_disease <- c(
    "CHROMOSOMAL_ABNORMALITIES",
    "INFLAMMATORY_BOWEL_DISEASES",
    "MULTIPLE_SCLEROSIS",
    "SCHIZOPHRENIA_AND_DELUSIONAL_DISEASES"
  )
  
  if (grepl("_M_", name_group)) {
    rm_disease <- c(rm_disease, "PROSTATE_DISEASES")
  }
  
  #   Remove diseases
  df <- df %>% select(all_of(setdiff(names(df), rm_disease)))
  
  # Calculate regularized tetrachoric of original data
  tetra_res <- tetrachoric(df, smooth = TRUE)
  matriz_r <- tetra_res$rho
  gamma_value <- 0.25
  red_ebic <- EBICglasso(matriz_r, n = nrow(df), gamma = gamma_value)
  colnames(red_ebic) <- colnames(df)
  rownames(red_ebic) <- colnames(df)
  
  #   Read null model
  null_model_rho <- read.csv(grep(name_group, null_model_files, value = T))$x
  threshold <- quantile(abs(null_model_rho), 0.95, na.rm = TRUE)
  
  #   Get components of the network
  edges_filter <- as.data.frame(as.table(red_ebic)) %>%
    rename(from = Var1,
           to = Var2,
           correlation = Freq) %>%
    filter(as.character(from) < as.character(to)) %>%
    mutate(weight = abs(correlation))  %>% 
    filter(weight > threshold)
  
  nodes_unique <- unique(c(edges_filter$from, edges_filter$to))
  prevalence <- colMeans(df[nodes_unique], na.rm = TRUE)
  
  nodes <- data.frame(
    id = nodes_unique,
    label = nodes_unique,
    prevalencia_pct = prevalence * 100,
    size = 10 + (prevalence * 100)
  )
  
  # #   Run Leiden
  edges_pos <- edges_filter %>% filter(correlation > 0)
  leiden_only_pos(edges = edges_pos, nodes = nodes, nombre_base = name_group, ruta_subcarpeta_estrato = paste0(path_folder_group, "/Leiden/", name_group))

  # # Run Spinnglass with only positive
  spinnglass_pos(edges = edges_pos, nodes = nodes, nombre_base = name_group, ruta_subcarpeta_estrato = paste0(path_folder_group, "/Spinnglass_pos/", name_group))

  # Run Spinnglass with negatives
  spinnglass_neg(edges = edges_filter, nodes = nodes, nombre_base = name_group, ruta_subcarpeta_estrato = paste0(path_folder_group, "/Spinnglass_neg/", name_group))
  
  
  
  
}
