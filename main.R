# Author: Alejandro Santos Mejías
# Date last update: 2026-06-11
# Script name: main.R
# Aims: Main script to orchestrate all steps


rm(list = ls())
gc()

source("98_packages.R")

step_0 <- T

if(step_0){
  source("p_step/step_0_null_model.R")
}
