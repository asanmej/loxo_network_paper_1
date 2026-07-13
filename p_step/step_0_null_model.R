# Author: Alejandro Santos Mejías
# Date last update: 2026-06-11
# Script name: step_0_null_model.R
# Aims: To run null model for each binary dataset.


rm(list = ls())
gc()

source("98_packages.R")

if (!dir.exists(path_output)) {
  dir.create(path_output, recursive = TRUE)
}

file_list <- list.files(path = path_data,
                       pattern = "^(2012|2017|2022)_(H|M)_(60-69|70-79|80-89|90)\\.csv$",
                       full.names = TRUE)

seed <- 29101997
set.seed(seed)

for (file_current in file_list) {
  
  name_group <- tools::file_path_sans_ext(basename(file_current))
  path_folder_group <- paste0(path_output, "null_model/")
  if (!dir.exists(path_folder_group)) {
    dir.create(path_folder_group, recursive = TRUE)
  }
  
  cat("Iniciando calculo del modelo nulo para ", name_group)
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
  
  #   Run null model
  gamma_value <- 0.25
  progressr::with_progress(null_model <- modelo_nulo_ebicglasso(X = df, B = 200, gamma_val = gamma_value, seed = seed))
  
  # Result null model 
  summary(null_model)
  # Save null model
  null_model_filename <- paste0(path_folder_group, "null_model_ebicglasso_", name_group, ".csv")
  write.csv(null_model, null_model_filename, fileEncoding = "UTF-8", row.names = FALSE)
}
