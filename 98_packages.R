# Author: Alejandro Santos Mejías
# Date last update: 2026-06-11
# Script name: 98_package.R
# Aims: Load package needed for LOXO Network analysis


##### PACKAGES
pkg <- c(
  "psych",
  "qgraph",
  "Matrix",
  "future",
  "furrr",
  "progressr",
  "qgraph",
  "tibble",
  "igraph",
  "readr",
  "visNetwork",
  "htmlwidgets",
  "dplyr"
)

lapply(pkg, function (x){if(!require(x, character.only = T)){install.packages(x , character.only = T)}})
rm(pkg)

# PATHING
source("99_path.R")

# FUNCTIONS 
functions <- list.files(path = path_funct, pattern = "\\.R$", full.names = TRUE, ignore.case = TRUE)
sapply(functions, source)




