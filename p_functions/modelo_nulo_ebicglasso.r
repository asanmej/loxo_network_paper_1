library(psych)
library(qgraph)
library(Matrix)
library(future)
library(furrr)
library(progressr)

modelo_nulo_ebicglasso <- function(X, B = 500, gamma_val = 0.25, seed = 29101997) {
  
  set.seed(seed)
  N <- nrow(X)
  p <- ncol(X)

  
  # Configurar la estrategia de paralelización para Windows
  num_cores <- availableCores() - 1
  plan(multisession, workers = num_cores)
  
  # Crear un canal de señalización de progreso para las B iteraciones
  p_bar <- progressor(steps = B)
  
  # Ejecución en paralelo optimizada con reporte de avance
  null_correlations <- future_map(1:B, function(b) {
    # Indicar que esta iteración ha comenzado/terminado para actualizar la barra
    p_bar(message = sprintf("Procesando permutación %d de %d", b, B))
    
    # Permutación directa por columnas
    X_perm <- as.data.frame(X)
    for (i in 1:p) {
      X_perm[, i] <- sample(X_perm[, i])
    }
    
    # Cálculo de la matriz tetracórica
    tc_res <- suppressWarnings(psych::tetrachoric(X_perm, smooth = TRUE))
    matriz_r <- tc_res$rho
    
    # Aplicar EBICglasso
    red_ebic_null <- qgraph::EBICglasso(matriz_r, n = N, gamma = gamma_val, verbose = F)
    
    valores_aristas <- red_ebic_null[upper.tri(red_ebic_null)]
    return(valores_aristas)
  }, .options = furrr_options(seed = TRUE))
  
  # Cerrar las conexiones en segundo plano
  plan(sequential)
  
  return(unlist(null_correlations))
}